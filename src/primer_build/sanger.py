"""Configurable, local Sanger validation against explicit edited references.

Coordinates in JSON/CSV reports are 1-based; internal positions are 0-based.
Peak ratios are screening heuristics, never calibrated mixture fractions.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import statistics
import struct

from Bio import Align, SeqIO
from Bio.Seq import Seq


DEFAULTS = dict(base_quality=20, mutation_quality=30, min_aligned_bases=150,
                min_identity=0.85, trim_error_probability=0.05,
                secondary_peak_ratio=0.33, min_peak_signal=100,
                trace_radius=1, alignment_padding=150, coverage_quality=30,
                min_bidirectional_fraction=0.90, indel_context_bases=10)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def natural_key(value):
    return tuple((0,int(p)) if p.isdigit() else (1,p.casefold())
                 for p in re.split(r'(\d+)',str(value)) if p)


def complement(base):
    return str(Seq(base).complement())


def mott_bounds(qualities, cutoff=0.05):
    """Maximum positive-score interval, retaining original trace coordinates."""
    total, best, start, bounds = 0.0, 0.0, 0, (0, 0)
    for i, quality in enumerate(qualities):
        total += cutoff - 10 ** (-quality / 10)
        if total > best:
            best, bounds = total, (start, i + 1)
        if total < 0:
            total, start = 0.0, i + 1
    return bounds


def primer_site(reference, primer):
    primer = primer.upper()
    if not primer or set(primer) - set('ACGT'):
        raise ValueError('Primer must contain only A/C/G/T')
    hits = []
    for strand, query in [(1, primer), (-1, str(Seq(primer).reverse_complement()))]:
        doubled = reference + reference[:len(primer)-1]
        hits.extend((i, strand) for i in range(len(reference))
                    if doubled.startswith(query, i))
    if len(hits) != 1:
        raise ValueError(f'Primer has {len(hits)} exact sites, expected one')
    return hits[0][0], len(primer), hits[0][1]


def directed_window(reference, site, length):
    """Reference after primer 3' end in raw sequencing orientation."""
    start, width, strand = site
    positions = [((start + width + i) if strand == 1 else (start - 1 - i))
                 % len(reference) for i in range(min(length, len(reference)-width))]
    sequence = ''.join(reference[p] if strand == 1 else complement(reference[p])
                       for p in positions)
    return sequence, positions


def read_reference(path):
    suffix = path.suffix.lower()
    fmt = 'snapgene' if suffix == '.dna' else 'genbank' if suffix in {'.gb', '.gbk', '.genbank'} else 'fasta'
    return SeqIO.read(path, fmt)


def build_reference(spec, base):
    path = (base / spec['reference']).resolve()
    record = read_reference(path)
    parental = str(record.seq).upper()
    if set(parental) - set('ACGT'):
        raise ValueError(f'Reference contains ambiguous bases: {path}')
    expected = list(parental)
    seen = set()
    for edit in spec.get('edits', []):
        pos = edit['position1'] - 1
        if pos in seen or not 0 <= pos < len(parental):
            raise ValueError('Duplicate or out-of-bounds edit')
        if len(edit['ref']) != 1 or len(edit['alt']) != 1 or edit['alt'] not in 'ACGT':
            raise ValueError('This reference-edit schema supports single-base substitutions only')
        if parental[pos] != edit['ref']:
            raise ValueError(f'Parental allele mismatch at {path}:{pos+1}')
        seen.add(pos)
        expected[pos] = edit['alt']
    required = set()
    for region in spec['required_regions']:
        if not 1 <= region['start1'] <= region['end1'] <= len(parental):
            raise ValueError('Invalid required region (use separate regions across the origin)')
        required.update(range(region['start1'] - 1, region['end1']))
    if not required or not seen <= required:
        raise ValueError('Required regions must contain every intended edit')
    return dict(path=path, source_sha256=sha256(path), record=record, parental=parental,
                expected=''.join(expected), required=required, edits=spec.get('edits', []),
                spec=spec)


def trace_evidence(raw, index, called, parameters):
    order = raw.get('FWO_1', b'')
    order = order.decode() if isinstance(order, bytes) else order
    peaks = raw.get('PLOC2')
    if set(order) != set('ACGT') or len(order) != 4 or peaks is None:
        return None
    channels = {base: raw.get(f'DATA{9+i}') for i, base in enumerate(order)}
    if any(channel is None for channel in channels.values()) or index >= len(peaks):
        return None
    center, radius = int(peaks[index]), parameters['trace_radius']
    heights = {base: max(channel[max(0, center-radius):min(len(channel), center+radius+1)], default=0)
               for base, channel in channels.items()}
    if called not in heights:
        return dict(center=center, heights=heights, secondary_ratio=None, secondary_base=None)
    alt = max((b for b in 'ACGT' if b != called), key=lambda b: heights[b])
    ratio = heights[alt] / max(heights[called], 1)
    return dict(center=center, heights=heights, secondary_ratio=round(ratio, 4),
                secondary_base=alt, primary_signal=heights[called])


def valid_trace(raw, sequence_length):
    order=raw.get('FWO_1',b'')
    order=order.decode() if isinstance(order,bytes) else order
    peaks=raw.get('PLOC2')
    channels=[raw.get(f'DATA{i}') for i in range(9,13)]
    if (len(order)!=4 or set(order)!=set('ACGT') or peaks is None
            or len(peaks)!=sequence_length or not peaks
            or any(channel is None for channel in channels)):
        return False
    lengths={len(channel) for channel in channels}
    return (len(lengths)==1 and min(peaks)>=0 and max(peaks)<len(channels[0])
            and all(a<=b for a,b in zip(peaks,peaks[1:])))


def align_read(path, reference, primer, parameters):
    record = SeqIO.read(path, 'abi')
    sequence = str(record.seq).upper()
    qualities = record.letter_annotations.get('phred_quality', [])
    if len(qualities) != len(sequence):
        raise ValueError('AB1 has no matching Phred quality array')
    lo, hi = mott_bounds(qualities, parameters['trim_error_probability'])
    if hi-lo < parameters['min_aligned_bases']:
        raise ValueError('Read too short after quality trimming')
    site = primer_site(reference['expected'], primer)
    window, positions = directed_window(reference['expected'], site,
                                         hi+parameters['alignment_padding'])
    query = sequence[lo:hi]
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-4,
                                    open_gap_score=-8, extend_gap_score=-1)
    # Overlap alignment: unsequenced reference and residual query ends are free.
    aligner.left_deletion_score = aligner.right_deletion_score = 0
    aligner.left_insertion_score = aligner.right_insertion_score = 0
    alignment = aligner.align(window, query)[0]
    observations, gaps = {}, []
    matches, aligned = 0, 0
    raw = record.annotations.get('abif_raw', {})
    blocks = list(zip(alignment.coordinates.T, alignment.coordinates.T[1:]))
    for block_index, ((a, c), (b, d)) in enumerate(blocks):
        a, b, c, d = map(int, (a, b, c, d))
        if b > a and d > c:
            for offset in range(b-a):
                refpos, readindex = positions[a+offset], lo+c+offset
                rawbase = sequence[readindex]
                base = rawbase if site[2] == 1 else complement(rawbase)
                peak = trace_evidence(raw, readindex, rawbase, parameters)
                if peak and peak.get('secondary_base') and site[2] == -1:
                    peak['secondary_base'] = complement(peak['secondary_base'])
                    peak['heights'] = {complement(k): v for k, v in peak['heights'].items()}
                mixed = bool(peak and peak.get('primary_signal', 0) >= parameters['min_peak_signal']
                             and peak.get('secondary_ratio', 0) >= parameters['secondary_peak_ratio'])
                observations[refpos] = dict(base=base, quality=qualities[readindex],
                    read_index0=readindex, strand=site[2], peak=peak, mixed=mixed)
                matches += base == reference['expected'][refpos]
                aligned += 1
        elif block_index not in (0, len(blocks)-1):
            # Only internal gaps are biological candidates. End clipping is not a deletion.
            context=parameters['indel_context_bases']
            neighborhood = qualities[max(lo,lo+c-context):min(hi,lo+d+context)]
            quality = min(neighborhood, default=0)
            if b > a:
                affected = positions[a:b]
                gaps.append(dict(kind='DEL', positions1=[p+1 for p in affected],
                    position1=min(affected)+1, ref=''.join(reference['expected'][p] for p in sorted(affected)),
                    alt='', quality=quality, strand=site[2], read_index0=lo+c))
            elif d > c:
                boundary = (positions[a] if site[2] == 1 else positions[a]+1) % len(reference['expected'])
                allele = query[c:d] if site[2] == 1 else str(Seq(query[c:d]).reverse_complement())
                gaps.append(dict(kind='INS', positions1=[boundary+1], position1=boundary+1,
                    ref='', alt=allele, quality=quality, strand=site[2], read_index0=lo+c))
    identity = matches / max(aligned, 1)
    if aligned < parameters['min_aligned_bases'] or identity < parameters['min_identity']:
        raise ValueError(f'Insufficient alignment: {aligned} bases, identity={identity:.3f}')
    seq_path = path.with_suffix('.seq')
    sidecar_matches = None
    if seq_path.exists():
        sidecar = ''.join(line.strip() for line in seq_path.read_text().splitlines()
                          if not line.startswith('>')).upper().replace(' ', '')
        sidecar_matches = sidecar == sequence
    summary = dict(path=str(path), length=len(sequence), trim_start0=lo, trim_end0=hi,
                   aligned_bases=aligned, identity=round(identity,5), strand=site[2],
                   median_quality=statistics.median(qualities[lo:hi]),
                   trace_available=valid_trace(raw,len(sequence)),
                   seq_sidecar_matches=sidecar_matches)
    return dict(summary=summary, observations=observations, gaps=gaps)


def aggregate_clone(sample, construct, colony, reads, reference, parameters):
    """Never let a low-quality agreeing read override high-quality disagreement."""
    expected, required = reference['expected'], reference['required']
    pile = defaultdict(list)
    for read in reads:
        for pos, obs in read['observations'].items():
            pile[pos].append(dict(obs, read_name=read['name']))
    issues, evidence_rows, missing, single, trace_flags = [], [], [], [], []
    outside_findings=[]
    hq = parameters['base_quality']
    good_positions = set()
    for pos in sorted(required):
        calls = [x for x in pile[pos] if x['quality'] >= hq]
        clean = [x for x in calls if x['base'] == expected[pos]]
        directions = {x['strand'] for x in clean}
        covered = [x for x in clean if x['quality'] >= parameters['coverage_quality']]
        if not covered:
            missing.append(pos+1)
        else:
            good_positions.add(pos)
        if directions != {-1,1}:
            single.append(pos+1)
        bad = [x for x in calls if x['base'] != expected[pos]
               and x['quality'] >= parameters['coverage_quality']]
        mixed_by_allele = defaultdict(set)
        for x in calls:
            if x['mixed']:
                mixed_by_allele[x['peak']['secondary_base']].add(x['strand'])
                trace_flags.append(dict(position1=pos+1,read=x['read_name'],quality=x['quality'],
                                       peak=x['peak']))
        bidirectional_mixed = any(v == {-1,1} for v in mixed_by_allele.values())
        if bad or bidirectional_mixed:
            alleles = Counter(x['base'] for x in calls)
            alt_directions = {x['strand'] for x in bad if x['base'] != expected[pos]}
            hard = len(alleles)==1 and len(alt_directions)==2 and all(x['base'] in 'ACGT' for x in calls)
            issues.append(dict(kind='SNV_OR_MIXED', position1=pos+1,
                expected=expected[pos], calls=dict(alleles), confirmed_difference=hard,
                mixed_flag=bidirectional_mixed,
                details=[dict(read=x['read_name'],base=x['base'],quality=x['quality'],
                              mixed=x['mixed'],peak=x['peak']) for x in calls]))
        evidence_rows.append(dict(position1=pos+1,expected=expected[pos],
            confident_reads=len(calls),clean_strands=sorted(directions),
            calls=[dict(read=x['read_name'],**{k:x[k] for k in ['base','quality','read_index0','mixed','peak','strand']})
                   for x in pile[pos]]))
    gap_issues = [dict(gap,read=read['name'],confirmed_difference=False)
                  for read in reads for gap in read['gaps']
                  if gap['quality'] >= parameters['coverage_quality']
                  and any(p-1 in required for p in gap['positions1'])]
    # Two-direction support confirms a gap; one-direction gaps still require review.
    for gap in gap_issues:
        gap['confirmed_difference'] = any(other['kind']==gap['kind'] and other['ref']==gap['ref']
            and other['alt']==gap['alt'] and other['position1']==gap['position1']
            and other['strand'] != gap['strand'] for other in gap_issues)
    issues.extend(gap_issues)
    outside_gaps=[dict(gap,read=read['name'],confirmed_difference=False,
                      outside_region=True)
                  for read in reads for gap in read['gaps']
                  if gap['quality']>=parameters['coverage_quality']
                  and not any(p-1 in required for p in gap['positions1'])]
    for gap in outside_gaps:
        gap['confirmed_difference']=any(other['kind']==gap['kind'] and other['ref']==gap['ref']
            and other['alt']==gap['alt'] and other['position1']==gap['position1']
            and other['strand']!=gap['strand'] for other in outside_gaps)
        if gap['confirmed_difference']:
            issues.append(gap)
    outside_findings.extend(outside_gaps)
    for pos in sorted(set(pile)-required):
        calls=[x for x in pile[pos] if x['quality']>=parameters['coverage_quality']]
        bad=[x for x in calls if x['base']!=expected[pos]]
        if not bad:
            continue
        confirmed=any({x['strand'] for x in bad if x['base']==base}=={-1,1}
                      for base in 'ACGT' if base!=expected[pos])
        finding=dict(kind='OUTSIDE_REGION_SNV',position1=pos+1,expected=expected[pos],
                     confirmed_difference=confirmed,
                     calls=[dict(read=x['read_name'],base=x['base'],quality=x['quality'],strand=x['strand'],
                                 read_index0=x['read_index0']) for x in calls])
        outside_findings.append(finding)
        if confirmed:
            issues.append(finding)
    edit_results = []
    for edit in reference['edits']:
        pos=edit['position1']-1
        calls=pile[pos]
        clean=[x for x in calls if x['quality'] >= parameters['mutation_quality']
               and x['base']==edit['alt'] and not x['mixed']]
        edit_results.append(dict(**edit, quality_threshold=parameters['mutation_quality'],
             bidirectional_hq={x['strand'] for x in clean}=={-1,1},
             observations=[dict(read=x['read_name'],**{k:x[k] for k in ['base','quality','mixed','peak','strand','read_index0']})
                           for x in calls]))
    error_reads = [r['name'] for r in reads if r['summary'].get('error')]
    if any(x['confirmed_difference'] for x in issues):
        status='REJECT_DIFFERENCE'
    elif issues:
        status='REVIEW_TRACE_OR_DISCORDANCE'
    elif missing:
        status='INCOMPLETE_COVERAGE'
    elif (1-len(single)/len(required) < parameters['min_bidirectional_fraction']
          or error_reads or not all(e['bidirectional_hq'] for e in edit_results)):
        status='SINGLE_STRAND_OR_LOW_Q'
    elif not all(r['summary'].get('trace_available') and r['summary'].get('seq_sidecar_matches') is not False for r in reads):
        status='REVIEW_INPUT_OR_TRACE'
    elif not single:
        status='PASS_BIDIRECTIONAL'
    else:
        status='PASS_COMPLEMENTARY'
    region_status=status
    if status.startswith('PASS_') and outside_findings:
        status='REVIEW_OUTSIDE_REGION'
    target_qualities=[x['quality'] for p in required for x in pile[p]
                      if x['quality']>=hq and x['base']==expected[p]]
    ranking_q = min((max(x['quality'] for x in pile[p]
                        if x['quality']>=hq and x['base']==expected[p])
                     for p in required if p in good_positions), default=0)
    return dict(sample=sample,construct=str(construct),colony=str(colony),status=status,
        region_status=region_status,
        required_bases=len(required),covered_expected_bases=len(good_positions),
        bidirectional_clean_bases=len(required)-len(single),
        target_min_quality=ranking_q,target_median_quality=statistics.median(target_qualities) if target_qualities else 0,
        missing_positions1=missing,single_strand_positions1=single,
        issues=issues,trace_flags=trace_flags,outside_region_findings=outside_findings,
        edit_results=edit_results,evidence=evidence_rows,
        reads=[dict(name=r['name'],gaps=r['gaps'],**r['summary']) for r in reads])


def analyze(config_path):
    config_path=Path(config_path).resolve()
    config=json.loads(config_path.read_text())
    base=config_path.parent
    if config.get('schema_version') != 1:
        raise ValueError('Unsupported schema_version')
    parameters=dict(DEFAULTS,**config.get('parameters',{}))
    unknown=set(config.get('parameters',{}))-set(DEFAULTS)
    if unknown:
        raise ValueError(f'Unknown parameters: {unknown}')
    for key in ['base_quality','coverage_quality','mutation_quality']:
        if not isinstance(parameters[key],int) or not 0 <= parameters[key] <= 93:
            raise ValueError(f'Invalid Phred threshold: {key}')
    for key in ['min_identity','min_bidirectional_fraction']:
        if not 0 < parameters[key] <= 1:
            raise ValueError(f'{key} must be in (0,1]')
    if not 0 < parameters['trim_error_probability'] < 1:
        raise ValueError('Invalid trimming error probability')
    if not 0 < parameters['secondary_peak_ratio'] <= 1:
        raise ValueError('Secondary peak ratio must be in (0,1]')
    if (not isinstance(parameters['trace_radius'],int) or not 0 <= parameters['trace_radius'] <= 3
            or not isinstance(parameters['min_aligned_bases'],int) or parameters['min_aligned_bases'] <= 0
            or not isinstance(parameters['indel_context_bases'],int) or parameters['indel_context_bases']<2
            or parameters['min_peak_signal']<=0 or parameters['alignment_padding']<0):
        raise ValueError('Invalid trace radius or minimum alignment length')
    input_dir=(base/config['input_directory']).resolve()
    regex=re.compile(config['sample_regex'])
    if not {'sample','construct','colony','primer'}<=set(regex.groupindex):
        raise ValueError('sample_regex requires named groups: sample, construct, colony, primer')
    references={str(k):build_reference(v,base) for k,v in config['constructs'].items()}
    primers=config['primers']
    grouped=defaultdict(list)
    input_hashes={}
    skipped=[]
    identities={}
    validated_bindings=set()
    for path in sorted(p for p in input_dir.rglob('*') if p.is_file() and p.suffix.lower()=='.ab1'):
        match=regex.fullmatch(path.name)
        if not match:
            skipped.append(str(path));continue
        tags=match.groupdict()
        construct,sample,colony,primer=(tags[k] for k in ['construct','sample','colony','primer'])
        if construct not in references or primer not in primers:
            raise ValueError(f'Unknown construct or primer: {path.name}')
        if (construct,primer) not in validated_bindings:
            primer_site(references[construct]['expected'],primers[primer])
            validated_bindings.add((construct,primer))
        if sample in identities and identities[sample] != (construct,colony):
            raise ValueError(f'Sample identity conflict: {sample}')
        identities[sample]=(construct,colony)
        if any(x['name']==primer for x in grouped[sample]):
            raise ValueError(f'Duplicate sample/primer {sample}/{primer}; use explicit distinct primer IDs for replicate reads')
        input_hashes[str(path)]=sha256(path)
        if path.with_suffix('.seq').exists():
            input_hashes[str(path.with_suffix('.seq'))]=sha256(path.with_suffix('.seq'))
        try:
            read=align_read(path,references[construct],primers[primer],parameters)
        except (ValueError, KeyError, OSError, EOFError, struct.error) as exc:
            read=dict(summary=dict(path=str(path),error=str(exc)),observations={},gaps=[])
        read['name']=primer
        grouped[sample].append(read)
    if not grouped:
        raise ValueError('No AB1 files matched the configured filename regex')
    if skipped and not config.get('allow_unmatched_files',False):
        raise ValueError(f'Unmatched AB1 files require explicit mapping or allow_unmatched_files: {skipped}')
    clones=[]
    for sample, reads in grouped.items():
        construct,colony=identities[sample]
        clones.append(aggregate_clone(sample,construct,colony,reads,references[construct],parameters))
    clones.sort(key=lambda c:(natural_key(c['construct']),natural_key(c['colony'])))
    priorities={s:i for i,s in enumerate(['PASS_BIDIRECTIONAL','PASS_COMPLEMENTARY','SINGLE_STRAND_OR_LOW_Q',
        'INCOMPLETE_COVERAGE','REVIEW_TRACE_OR_DISCORDANCE','REVIEW_INPUT_OR_TRACE','REJECT_DIFFERENCE'])}
    priorities['REVIEW_OUTSIDE_REGION']=priorities['REVIEW_TRACE_OR_DISCORDANCE']
    for construct in references:
        rows=sorted((c for c in clones if c['construct']==construct),
                    key=lambda c:(priorities[c['status']],len(c['outside_region_findings']),
                                  len(c['trace_flags']),
                                  -c['bidirectional_clean_bases']/c['required_bases'],
                                  -c['target_min_quality'],-c['target_median_quality'],c['sample']))
        for rank,row in enumerate(rows,1):
            row['rank']=rank
            row['recommended']=rank==1 and row['status'].startswith('PASS_')
    missing_samples=sorted(set(config.get('expected_samples',[]))-set(grouped))
    return dict(config=config,parameters=parameters,references=references,clones=clones,
                input_hashes=input_hashes,missing_samples=missing_samples,skipped_files=skipped)
