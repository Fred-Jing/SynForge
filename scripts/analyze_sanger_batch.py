#!/usr/bin/env python3
"""Validate paired Sanger AB1 reads using an explicit batch JSON manifest."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import csv
import gzip
import html
import json
from pathlib import Path
import platform
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from Bio import SeqIO, __version__ as biopython_version
from primer_build.sanger import analyze, sha256, valid_trace, natural_key

LABELS={'PASS_BIDIRECTIONAL':'双向逐碱基通过',
        'PASS_COMPLEMENTARY':'双向互补覆盖通过',
        'REJECT_DIFFERENCE':'有明确非预期差异',
        'INCOMPLETE_COVERAGE':'高质量覆盖不足',
        'SINGLE_STRAND_OR_LOW_Q':'单向支持或突变位点质量不足',
        'REVIEW_TRACE_OR_DISCORDANCE':'峰形或双向不一致需复核',
        'REVIEW_INPUT_OR_TRACE':'输入或峰图需复核',
        'REVIEW_OUTSIDE_REGION':'目标区域通过，外侧差异需复核'}


def write_csv(path, rows, fields):
    with path.open('w',encoding='utf-8-sig',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def display_path(path):
    path=Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def trace_panel(path, read_index, radius=8):
    record=SeqIO.read(path,'abi');raw=record.annotations['abif_raw']
    peaks=raw.get('PLOC2');order=raw.get('FWO_1',b'')
    order=order.decode() if isinstance(order,bytes) else order
    if not valid_trace(raw,len(record)) or not 0<=read_index<len(record):
        return dict(unavailable=True)
    left=max(0,read_index-radius);right=min(len(record),read_index+radius+1)
    x0=max(0,int(peaks[left])-6);x1=min(len(raw['DATA9']),int(peaks[right-1])+7)
    return dict(sequence=str(record.seq)[left:right],
        qualities=record.letter_annotations['phred_quality'][left:right],
        positions=[int(p)-x0 for p in peaks[left:right]],
        signals={base:list(raw[f'DATA{9+i}'][x0:x1]) for i,base in enumerate(order)},
        selected_index=read_index-left,raw_read_start1=left+1)


def describe(issue):
    prefix='外侧 ' if issue.get('outside_region') else ''
    if issue['kind']=='DEL':
        return f"{prefix}{issue['position1']} 缺失 {issue['ref']}"
    if issue['kind']=='INS':
        return f"{prefix}{issue['position1']} 前插入 {issue['alt']}"
    if issue['kind']=='OUTSIDE_REGION_SNV':
        return f"{issue['position1']} 外侧预期 {issue['expected']}："+', '.join(
            f"{x['read']}={x['base']}(Q{x['quality']})" for x in issue['calls'])
    return f"{issue['position1']} 预期 {issue['expected']}，读段 {issue['calls']}"


def render_report_html(result, summary, output):
    payload=[]
    for clone in result['clones']:
        selected=set(e['position1'] for e in clone['edit_results'])
        selected.update(i['position1'] for i in clone['issues'])
        selected.update(i['position1'] for i in clone['outside_region_findings'])
        selected.update(i['position1'] for i in clone['trace_flags'])
        if not selected:
            selected.add(clone['evidence'][len(clone['evidence'])//2]['position1'])
        snapshots={}
        paths={r['name']:r['path'] for r in clone['reads']}
        evidence_by_pos={e['position1']:e for e in clone['evidence']}
        for pos in sorted(selected):
            panels=[]
            entry=evidence_by_pos.get(pos)
            if entry:
                for call in entry['calls']:
                    panels.append(dict(read=call['read'],base=call['base'],quality=call['quality'],
                        strand=call['strand'],**trace_panel(paths[call['read']],call['read_index0'])))
            # Deleted reference bases have no mapped query base: show gap boundary.
            for issue in clone['issues']+clone['outside_region_findings']:
                if issue['position1']==pos and issue['kind'] in {'DEL','INS'}:
                    if not any(p['read']==issue['read'] for p in panels):
                        panels.append(dict(read=issue['read'],base='gap',quality=issue['quality'],
                            strand=issue['strand'],**trace_panel(paths[issue['read']],issue['read_index0'])))
            if not panels:
                for finding in clone['outside_region_findings']:
                    if finding['position1']!=pos:continue
                    # Outside-region read index is stored by the core for trace review.
                    for call in finding.get('calls',[]):
                        if 'read_index0' in call:
                            panels.append(dict(read=call['read'],base=call['base'],quality=call['quality'],
                                strand=call['strand'],**trace_panel(paths[call['read']],call['read_index0'])))
            snapshots[str(pos)]=panels
        payload.append(dict(sample=clone['sample'],construct=clone['construct'],status=clone['status'],
            label=LABELS[clone['status']],recommended=clone['recommended'],
            usable=clone['status'].startswith('PASS_'),rank=clone['rank'],
            covered=clone['covered_expected_bases'],required=clone['required_bases'],
            bidirectional=clone['bidirectional_clean_bases'],qmin=clone['target_min_quality'],
            trace_flags=len(clone['trace_flags']),edit_positions=[e['position1'] for e in clone['edit_results']],
            issues=[describe(i) for i in clone['issues']],
            outside=[describe(i) for i in clone['outside_region_findings']],
            snapshots=snapshots,evidence=[dict(position1=e['position1'],expected=e['expected'],
                calls=[{k:c[k] for k in ['read','base','quality','mixed']} for c in e['calls']]) for e in clone['evidence']]))
    data=json.dumps(payload,ensure_ascii=False).replace('<','\\u003c')
    title=html.escape(result['config']['batch_name'])
    document='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>__TITLE__</title>
<style>body{font:15px system-ui,sans-serif;margin:25px;color:#263443;max-width:1500px}p{line-height:1.65}.muted{color:#667}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px;border:1px solid #dbe2e9;text-align:left}th{background:#edf2f7;position:sticky;top:0}.pass{background:#edf8ef}.review{background:#fff5de}.reject{background:#ffe9e9}button,select{padding:7px;margin:5px}tr{cursor:pointer}svg{width:100%;border:1px solid #ddd;background:white}h2{margin-top:30px}.base{font:13px monospace}.bad{background:#ffdfd9}.q{color:#789;font-size:11px}.controls{position:sticky;top:0;background:white;z-index:3;padding:6px}</style>
<h1>__TITLE__</h1><p>__SUMMARY__</p>
<p>筛选范围为配置中的必验区域。每个位点至少一个 Q30 匹配，≥90% 位点双向 Q20 一致；预期突变必须双向 Q30。互补覆盖通过不代表每个碱基都双向高质量，也不代表全质粒测序通过。单向副峰为人工复核线索，不是混样比例。</p>
<label><input type="checkbox" id="only"/>只显示全部达标样本（不限数量）</label><table><thead><tr><th>样本</th><th>判定</th><th>排名</th><th>目标覆盖</th><th>双向质量支持位点</th><th>每位点最佳 Q 的最低值</th><th>单向副峰线索数</th><th>差异/提示</th></tr></thead><tbody id="rows"></tbody></table>
<h2 id="detailtitle">峰图复核</h2><div class="controls"><select id="sample"></select><select id="pos"></select><label><input id="allbases" type="checkbox"/>显示整个验证区域</label></div>
<p id="note"></p><div id="traces"></div><h3>逐碱基证据（参考方向）</h3><div id="bases"></div>
<p class="muted">原始峰图保留测序读取方向；反向读段峰图上的 G 对应参考方向的 C。下方逐碱基表已统一到参考正链。低质量尾端、单向差异和缺失数据均保留在 evidence.json.gz / CSV 中。完整方法与限制见 report.md。</p>
<script>const DATA=__DATA__; const PRIMERS=__PRIMERS__; const esc=x=>String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const byId=id=>document.getElementById(id); const sample=byId('sample'), pos=byId('pos');
for(const d of DATA){const o=document.createElement('option');o.value=d.sample;o.textContent=d.sample+' · '+d.label;sample.appendChild(o)}
function cls(d){return d.status.startsWith('PASS_')?'pass':d.status.startsWith('REJECT')?'reject':'review'}
function rows(){byId('rows').innerHTML=DATA.filter(d=>!byId('only').checked||d.usable).map(d=>'<tr class="'+cls(d)+'" data-id="'+esc(d.sample)+'"><td>'+esc(d.sample)+(d.recommended?' ★':'')+'</td><td>'+d.label+'</td><td>'+d.rank+'</td><td>'+d.covered+'/'+d.required+'</td><td>'+d.bidirectional+'/'+d.required+'</td><td>'+d.qmin+'</td><td>'+d.trace_flags+'</td><td>'+esc([...d.issues,...d.outside].join('; '))+'</td></tr>').join('');for(const row of byId('rows').rows)row.onclick=()=>{sample.value=row.dataset.id;choose()}}
function trace(p){if(p.unavailable)return '<p>缺少原始通道</p>';const W=1100,H=210,top=40,bottom=175, N=Math.max(...Object.values(p.signals).map(x=>x.length)),max=Math.max(1,...Object.values(p.signals).flat()),xx=x=>20+x/Math.max(1,N-1)*(W-40),yy=y=>bottom-y/max*(bottom-top);const colors={A:'#159447',C:'#246dcc',G:'#222',T:'#d83b3b'};let s='<h3>'+esc(p.read)+'：参考正链 '+esc(p.base)+'，Q'+p.quality+'；原始读取方向 '+(p.strand===1?'正向':'反向')+'</h3><svg viewBox="0 0 '+W+' '+H+'">';s+='<rect x="'+(xx(p.positions[p.selected_index])-8)+'" y="25" width="16" height="155" fill="#fff0ba"/>';for(const [b,v] of Object.entries(p.signals))s+='<polyline fill="none" stroke="'+colors[b]+'" stroke-width="1.6" points="'+v.map((y,i)=>xx(i).toFixed(1)+','+yy(y).toFixed(1)).join(' ')+'"/>';for(let i=0;i<p.sequence.length;i++)s+='<text x="'+xx(p.positions[i])+'" y="18" text-anchor="middle" fill="'+(colors[p.sequence[i]]||'#777')+'">'+esc(p.sequence[i])+'</text><text x="'+xx(p.positions[i])+'" y="200" text-anchor="middle" font-size="11">'+p.qualities[i]+'</text>';return s+'</svg>'}
function detail(){const d=DATA.find(x=>x.sample===sample.value),p=pos.value;byId('detailtitle').textContent=d.sample+' 峰图与碱基复核';byId('note').textContent=d.label+'。参考坐标 '+p+'；黄色标记为所选位置或缺失边界。';byId('traces').innerHTML=(d.snapshots[p]||[]).map(trace).join('')||'<p>该位置没有可展示的峰图。</p>';const pp=Number(p),ev=d.evidence.filter(e=>byId('allbases').checked||Math.abs(e.position1-pp)<=12);byId('bases').innerHTML='<table><thead><tr><th>坐标</th><th>预期</th><th>T7</th><th>T7T</th></tr></thead><tbody>'+ev.map(e=>'<tr><td>'+e.position1+'</td><td>'+e.expected+'</td>'+['T7','T7T'].map(name=>{const c=e.calls.find(x=>x.read===name);return '<td class="'+(c&&c.quality>=30&&c.base!==e.expected?'bad':'')+'">'+(c?esc(c.base)+' <span class="q">Q'+c.quality+(c.mixed?' · 单向副峰':'')+'</span>':'未覆盖')+'</td>'}).join('')+'</tr>').join('')+'</tbody></table>'}
function choose(){const d=DATA.find(x=>x.sample===sample.value);pos.innerHTML=Object.keys(d.snapshots).map(p=>'<option value="'+p+'">'+p+(d.edit_positions.includes(Number(p))?' · 预期突变':' · 复核线索')+'</option>').join('');pos.value=String(d.edit_positions[0]||Object.keys(d.snapshots)[0]);detail()}
sample.onchange=choose;pos.onchange=detail;byId('allbases').onchange=detail;byId('only').onchange=rows;rows();choose();</script></html>'''
    document=document.replace('<th>T7</th><th>T7T</th>',"'+PRIMERS.map(name=>'<th>'+esc(name)+'</th>').join('')+'")
    document=document.replace("['T7','T7T'].map(name=>",'PRIMERS.map(name=>')
    document=document.replace('__PRIMERS__',json.dumps(list(result['config']['primers'])))
    parameters=result['parameters']
    document=document.replace('每个位点至少一个 Q30 匹配，≥90% 位点双向 Q20 一致；预期突变必须双向 Q30。',
        f"每个位点至少一个 Q{parameters['coverage_quality']} 匹配，≥{parameters['min_bidirectional_fraction']*100:g}% 位点双向 Q{parameters['base_quality']} 一致；预期突变必须双向 Q{parameters['mutation_quality']}。")
    document=document.replace('c.quality>=30',f"c.quality>={parameters['coverage_quality']}")
    document=document.replace('__TITLE__',title).replace('__SUMMARY__',html.escape(summary)).replace('__DATA__',data)
    output.write_text(document)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--overwrite',action='store_true',help='Replace generated reports in a managed result directory')
    args=parser.parse_args();config_path=args.config.resolve()
    output=(args.output or config_path.parent/'results').resolve()
    config=json.loads(config_path.read_text());config_digest=sha256(config_path)
    input_dir=(config_path.parent/config['input_directory']).resolve()
    if output==input_dir or input_dir in output.parents:
        parser.error('Output must be outside the original sequencing input directory')
    if output.exists() and any(output.iterdir()):
        manifest_path=output/'manifest.json'
        if not args.overwrite or not manifest_path.exists():
            parser.error('Output is not empty; use a new directory, or --overwrite for a managed result directory')
        previous=json.loads(manifest_path.read_text())
        if previous.get('workflow')!='sanger_batch_v1' or previous.get('config_path')!=str(config_path):
            parser.error('Output is not managed by this workflow/config')
    result=analyze(config_path);clones=result['clones'];counts=Counter(c['status'] for c in clones)
    parameters=result['parameters']
    output.mkdir(parents=True,exist_ok=True)
    with gzip.open(output/'evidence.json.gz','wt',encoding='utf-8') as handle:
        json.dump(clones,handle,ensure_ascii=False)
    fields=['sample','construct','colony','status','status_label','usable','recommended','rank',
            'required_bases','covered_expected_bases','bidirectional_clean_bases','bidirectional_fraction',
            'target_min_quality','target_median_quality','target_mutation_bidir_hq','trace_flag_count',
            'outside_region_flags','details']
    sample_rows=[]
    for c in clones:
        sample_rows.append(dict(c,status_label=LABELS[c['status']],usable=c['status'].startswith('PASS_'),
            target_mutation_bidir_hq=all(e['bidirectional_hq'] for e in c['edit_results']),
            trace_flag_count=len(c['trace_flags']),outside_region_flags=len(c['outside_region_findings']),
            bidirectional_fraction=round(c['bidirectional_clean_bases']/c['required_bases'],4),
            details='; '.join(dict.fromkeys(describe(i) for i in c['issues']+c['outside_region_findings']))))
    write_csv(output/'samples.csv',sample_rows,fields)
    usable_rows=[row for row in sample_rows if row['usable']]
    write_csv(output/'usable_samples.csv',usable_rows,fields)
    selection=['# 全部质量达标样本','',
        f"共 {len(usable_rows)} 个；不限制每个构建的数量。排名仅供参考，不作为额外排除条件。",'',
        '|构建|全部可用样本（按原编号）|数量|','|---|---|---:|']
    for construct in sorted(result['references'],key=natural_key):
        rows=sorted((r for r in usable_rows if r['construct']==construct),key=lambda r:natural_key(r['colony']))
        selection.append(f"|{construct}|{', '.join(r['sample'] for r in rows) or '无'}|{len(rows)}|")
    selection+=['','这些样本达到当前配置的目标区序列及质量要求，可作为后续转化/表达的候选；所有样本均保留原始名称。验证范围以 batch.json 的 required_regions 为准，不代表全质粒或表达性能已经验证。','',
                '未通过、待复核和缺失样本见 report.md；逐样本指标见 usable_samples.csv。','']
    (output/'usable_samples.md').write_text('\n'.join(selection))
    read_rows=[]
    for c in clones:
        read_rows.extend(dict(sample=c['sample'],construct=c['construct'],**r) for r in c['reads'])
    write_csv(output/'reads.csv',read_rows,['sample','construct','name','path','length','trim_start0',
        'trim_end0','aligned_bases','identity','strand','median_quality','trace_available','seq_sidecar_matches','error'])
    issues=[]
    for c in clones:
        for item in c['issues']+c['outside_region_findings']:
            issues.append(dict(sample=c['sample'],construct=c['construct'],status=c['status'],
                kind=item['kind'],position1=item['position1'],confirmed_difference=item['confirmed_difference'],
                description=describe(item)))
    issues=[dict(items) for items in dict.fromkeys(tuple(row.items()) for row in issues)]
    write_csv(output/'variants.csv',issues,['sample','construct','status','kind','position1','confirmed_difference','description'])
    summaries=[]
    for construct in result['references']:
        rows=sorted((c for c in clones if c['construct']==construct),key=lambda c:c['rank'])
        passed=[c['sample'] for c in rows if c['status'].startswith('PASS_')]
        summaries.append(dict(construct=construct,received=len(rows),passed=len(passed),
            recommended=next((c['sample'] for c in rows if c['recommended']),''),
            backups=', '.join(passed[1:]),review=', '.join(c['sample'] for c in rows if c['status'].startswith('REVIEW') or c['status'] in {'INCOMPLETE_COVERAGE','SINGLE_STRAND_OR_LOW_Q'}),
            rejected=', '.join(c['sample'] for c in rows if c['status'].startswith('REJECT'))))
    write_csv(output/'construct_summary.csv',summaries,list(summaries[0]))
    with (output/'expected_references.fasta').open('w') as handle:
        for key,ref in result['references'].items():
            handle.write(f'>{key} expected_full_reference_not_sequencing_consensus\n{ref["expected"]}\n')
    with (output/'observed_validated_regions.fasta').open('w') as handle:
        for c in clones:
            # Observed best-quality allele, not a reference-filled assembly.
            sequence=[]
            for e in c['evidence']:
                calls=sorted((x for x in e['calls'] if x['quality']>=result['parameters']['coverage_quality']),key=lambda x:-x['quality'])
                alleles={x['base'] for x in calls}
                sequence.append(calls[0]['base'] if calls and len(alleles)==1 else 'N')
            handle.write(f'>{c["sample"]} {c["status"]} required_positions_only_gaps_as_N_insertions_in_variants\n'+''.join(sequence)+'\n')
    codes=[ROOT/'src/primer_build/sanger.py',Path(__file__).resolve(),ROOT/'uv.lock']
    manifest=dict(workflow='sanger_batch_v1',config_path=str(config_path),config_sha256=config_digest,
        input_directory=str(input_dir),input_files=[dict(path=p,sha256=h) for p,h in sorted(result['input_hashes'].items())],
        references=[dict(construct=k,path=str(r['path']),source_sha256=r['source_sha256'],
            expected_sequence_sha256=__import__('hashlib').sha256(r['expected'].encode()).hexdigest(),
            edits=r['edits'],required_regions=r['spec']['required_regions']) for k,r in result['references'].items()],
        code=[dict(path=str(p),sha256=sha256(p)) for p in codes],
        software=dict(python=platform.python_version(),biopython=biopython_version),
        parameters=result['parameters'],created_at=datetime.now().astimezone().isoformat(),
        counts=dict(counts),samples_received=len(clones),reads_received=len(read_rows),
        missing_samples=result['missing_samples'],skipped_files=result['skipped_files'])
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    passed=sum(c['status'].startswith('PASS_') for c in clones)
    summary=f'{len(clones)} 个克隆，{len(read_rows)} 条 Sanger 读段；{passed} 个通过目标区域筛选。未收到：'+(', '.join(result['missing_samples']) or '无')
    lines=[f'# {config["batch_name"]}', '',summary,'',
        '## 可选克隆','','推荐来自本批目标区域的双向互补证据，不能代替全质粒测序。', '',
        '|构建|收到|通过|优先克隆|其他可选（按证据排序）|待复核|不选|',
        '|---|---:|---:|---|---|---|---|']
    for s in summaries:
        lines.append('|'+ '|'.join(str(s[k]) or '—' for k in ['construct','received','passed','recommended','backups','review','rejected'])+'|')
    lines+=['','## 判定与异常','','- 采用配置指定的实际实验参考，仅应用 batch.json 中明确声明的预期替换；未指定编辑时直接验证原参考。原始参考和测序输入均未修改。',
        f"- 每个 required_region 位点至少有一个 Q{parameters['coverage_quality']} 的预期碱基，至少 {parameters['min_bidirectional_fraction']*100:g}% 位点有正反两方向 Q{parameters['base_quality']} 一致；每个预期突变必须双向 Q{parameters['mutation_quality']}，且突变位点没有达到本工作流线索阈值的副峰。",
        '- PASS_COMPLEMENTARY 是合并证据后通过，局部低质量碱基由另一方向补足；PASS_BIDIRECTIONAL 才表示全部 required_region 碱基均有双向 Q20。两类都不表示质粒全长经过验证。',
        '- 单向高质量非预期碱基和高质量 indel 均不直接放行；两方向支持的非预期等位基因/indel 标为明确差异。任何目标区外观察到的 Q30 差异也单列复核，避免已知疑点样本被优先选中。',
        '- 副峰筛查读取 ABIF FWO_1、DATA9–12、PLOC2：在碱基峰位 ±1 样点取各通道最大值，主峰≥100、次峰/主峰≥0.33 时记线索。孤立单向副峰不等于混样，保留用于排序/人工核对；两方向同位点同一副等位基因出现则需复核。峰高比例不是混合样本比例。',
        '- 质量线、双向覆盖比例和副峰线索阈值均为配置中明确记录的本地筛选策略，不能称为供应商通用标准。',
        '- expected_references.fasta 是期望序列；observed_validated_regions.fasta 来自实际读段，未知/冲突为 N，删除位点为 N，插入另见 variants/evidence，不作为已测全质粒共识。','',
        '### 未通过或需复核样本','']
    for c in clones:
        if c['status'].startswith('PASS_'):continue
        detail='; '.join(dict.fromkeys(describe(i) for i in c['issues']+c['outside_region_findings']))
        if c['missing_positions1']:detail+=f'；Q30 覆盖不足 {len(c["missing_positions1"])} 位点'
        if not all(e['bidirectional_hq'] for e in c['edit_results']):detail+='；预期突变未满足双向 Q30'
        lines.append(f'- **{c["sample"]}**：{LABELS[c["status"]]}。{detail}')
    lines+=['','## 方法与复用','','- 主输入为原始 AB1（序列、Phred、四通道峰图）；同名 .seq 只用于一致性检查，不当作独立测序证据。',
        '- Richard Mott 风格最大累计得分修剪（cutoff=0.05），保留原始峰图索引；按配置引物的唯一全质粒结合位点确定方向，做引物后窗口的半全局比对，并统一到参考正链坐标。',
        '- 比对：match +2、mismatch −4、gap-open −8、gap-extend −1；两端残余可自由裁切，内部 gap 单列检查。默认至少 150 个比对碱基、identity ≥0.85。实际 indel 的等价重复区位置可能不唯一，当前输出为比对位置；单向或位置不一致者保留复核。',
        '- 双向不一致不靠多数投票掩盖；低质量 N、缺读段、缺样本、错误序列/输入不能自动补成预期参考。',
        '- '+config.get('expected_samples_note','缺失项按配置中的 expected_samples 报告；未配置时不推断应收样本数。'),' ',
        '```sh',f'uv run scripts/analyze_sanger_batch.py --config "{display_path(config_path)}" --output "{display_path(output)}"',
        '```','',
        '首次输出须为空目录；同一配置管理的结果可显式加 --overwrite 重跑。未来换批次主要修改 batch.json 的原始目录、命名正则、引物、逐构建参考、预期突变和验证区域；不必修改脚本。',
        '', '参考：',
        '- [Biopython ABI parser / Mott trimming](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py)',
        '- [Thermo Fisher：Sanger 数据 QV20 指标](https://www.thermofisher.com/blog/behindthebench/your-sanger-sequencing-data-analysis-question-seq-it-out/)',
        '- [Thermo Fisher：混峰与背景排查](https://www.thermofisher.com/TFS-Assets/LSG/manuals/MAN0014435_Trbleshoot_Sanger_seq_data_UB.pdf)',
        '', '完整逐碱基证据在 evidence.json.gz；输入、参考及代码哈希在 manifest.json。report.html 可直接离线查看推荐、原始峰图和逐碱基质量。','']
    markdown='\n'.join(lines)
    markdown=markdown.replace('副峰筛查读取 ABIF FWO_1、DATA9–12、PLOC2：在碱基峰位 ±1 样点取各通道最大值，主峰≥100、次峰/主峰≥0.33 时记线索。',
        f"副峰筛查读取 ABIF FWO_1、DATA9–12、PLOC2：在碱基峰位 ±{parameters['trace_radius']} 样点取各通道最大值，主峰≥{parameters['min_peak_signal']}、次峰/主峰≥{parameters['secondary_peak_ratio']} 时记线索。")
    markdown+=f"\nIndel 的质量取 gap 两侧各 {parameters['indel_context_bases']} bp（插入时包含插入碱基）的最低质量；达到覆盖质量线才作为高质量候选，外侧 indel 同样检查。远端低质量 indel 不作为已证实变异，但保留于逐读段证据中。\n"
    markdown=markdown.replace('Q30 覆盖不足',f"Q{parameters['coverage_quality']} 覆盖不足")
    markdown=markdown.replace('未满足双向 Q30',f"未满足双向 Q{parameters['mutation_quality']}")
    markdown=markdown.replace('均有双向 Q20',f"均有双向 Q{parameters['base_quality']}")
    markdown=markdown.replace('任何目标区外观察到的 Q30 差异',f"任何目标区外观察到的 Q{parameters['coverage_quality']} 差异")
    markdown=markdown.replace('至少 150 个比对碱基、identity ≥0.85',f"至少 {parameters['min_aligned_bases']} 个比对碱基、identity ≥{parameters['min_identity']}")
    markdown=markdown.replace('cutoff=0.05',f"cutoff={parameters['trim_error_probability']}")
    (output/'report.md').write_text(markdown)
    render_report_html(result,summary,output/'report.html')
    for path,digest in result['input_hashes'].items():
        if sha256(path)!=digest:raise RuntimeError(f'Input changed during run: {path}')
    for reference in result['references'].values():
        if sha256(reference['path'])!=reference['source_sha256']:raise RuntimeError('Reference changed during run')
    if sha256(config_path)!=config_digest:raise RuntimeError('Configuration changed during run')
    print(summary)
    for s in summaries:print(f"{s['construct']}: {s['recommended'] or '无自动推荐'}; {s['passed']}/{s['received']} passed")
    print(f'Report: {output / "report.html"}')


if __name__=='__main__':
    main()
