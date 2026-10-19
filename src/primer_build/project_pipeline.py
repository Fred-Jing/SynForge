from __future__ import annotations

import re
import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
import signal
import time
from dataclasses import dataclass
from multiprocessing import TimeoutError as MultiprocessingTimeoutError
from multiprocessing import get_context
from pathlib import Path

from .backbone import BackboneModel, reconstruct_backbone
from .codon_opt import optimize_codons
from .expression import rank_expression_candidates
from .io import (
    get_poi_feature,
    load_rules,
    read_snapgene,
    read_structure_design,
    validate_structure_map_pair,
)
from .models import AssemblyPlan, Oligo, StructureDesign
from .planner import plan_assembly
from .snapgene import write_clean_snapgene
from .validate import validate_plan
from .workbooks import build_project_workbooks


@dataclass(frozen=True)
class NumberingScheme:
    map_prefix: str
    map_separator: str
    map_start: int
    primer_prefix: str
    primer_start: int

    def construct_id(self, offset: int) -> str:
        return f"{self.map_prefix}{self.map_separator}{self.map_start + offset}"

    def primer_design_id(self, offset: int) -> str:
        return f"{self.primer_prefix}{self.primer_start + offset}"


def parse_numbering(map_example: str, primer_example: str) -> NumberingScheme:
    map_match = re.fullmatch(r"(.+?)(-?)(\d+)", map_example.strip())
    if not map_match:
        raise ValueError(
            "--map-start must look like KYN-1 or ARS1: a prefix and integer"
        )
    primer_match = re.fullmatch(
        r"(.+?)(\d+)-1-1-FW",
        primer_example.strip(),
        flags=re.IGNORECASE,
    )
    if not primer_match:
        raise ValueError(
            "--primer-start must look like K1-1-1-FW; the numeric design "
            "identifier will be incremented for subsequent structures"
        )
    return NumberingScheme(
        map_prefix=map_match.group(1),
        map_separator=map_match.group(2),
        map_start=int(map_match.group(3)),
        primer_prefix=primer_match.group(1),
        primer_start=int(primer_match.group(2)),
    )


def _natural_key(path: Path) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", path.name)
    )


def discover_structures(project_dir: Path) -> list[Path]:
    structures_dir = project_dir / "structures"
    if not structures_dir.is_dir():
        raise ValueError(f"Missing structures directory: {structures_dir}")
    paths = sorted(
        (
            path
            for path in structures_dir.rglob("*")
            if path.is_file() and path.suffix.lower() in {".pdb", ".cif", ".mmcif"}
        ),
        key=_natural_key,
    )
    if not paths:
        raise ValueError(f"No PDB/CIF structures found under {structures_dir}")
    return paths


def _all_fragment_oligos(plan: AssemblyPlan) -> list[Oligo]:
    return [
        oligo
        for fragment in plan.synthetic_fragments
        for oligo in (
            fragment.core_forward,
            fragment.core_reverse,
            fragment.outer_forward,
            fragment.outer_reverse,
        )
    ]


def _plan_selection_key(plan: AssemblyPlan) -> tuple[object, ...]:
    fragment_oligos = _all_fragment_oligos(plan)
    flag_count = sum(bool(oligo.risk_flags) for oligo in fragment_oligos)
    long_count = sum(oligo.tier == "long" for oligo in fragment_oligos)
    max_effective_length = max(
        len(oligo.effective_sequence)
        for fragment in plan.synthetic_fragments
        for oligo in (fragment.outer_forward, fragment.outer_reverse)
    )
    return (
        len(plan.synthetic_fragments),
        flag_count,
        long_count,
        max_effective_length,
        plan.candidate_rank,
        plan.score,
    )


def _design_one(
    structure_path: Path,
    construct_id: str,
    primer_design_id: str,
    chain: str,
    backbone,
    rules: dict,
    max_candidates: int,
) -> tuple[StructureDesign, AssemblyPlan]:
    design_started = time.perf_counter()
    design = read_structure_design(
        structure_path,
        construct_id,
        protein_chain=chain,
        prepend_initiator_methionine=rules["input"][
            "prepend_initiator_methionine"
        ],
        ignored_comp_ids=rules["input"]["ignored_structure_comp_ids"],
        forced_polymer_comp_ids=rules["input"]["forced_polymer_comp_ids"],
        component_translation_overrides=rules["input"][
            "component_translation_overrides"
        ],
    )
    codon_candidates = optimize_codons(
        design.expression_sequence,
        rules["codon_optimization"],
        ncaa_codon=rules["input"]["ncaa_codon"],
        sd_like_motifs=rules["translation_initiation"]["sd_like_motifs"],
        sd_like_search_penalty=rules["codon_optimization"][
            "sd_like_search_penalty"
        ],
    )
    candidates = rank_expression_candidates(
        codon_candidates,
        backbone.sequence_without_poi[: backbone.poi_insertion_index],
        rules["translation_initiation"],
    )
    candidate_total = min(max_candidates, len(candidates))
    optimization_elapsed = time.perf_counter() - design_started
    print(
        f"[{construct_id}] sequence: {len(design.expression_sequence)} aa | "
        f"{len(design.expression_sequence) * 3} bp CDS | "
        f"nonstandard={','.join(f'{residue.component_id}@{residue.expression_position}' for residue in design.nonstandard_residues) or 'none'} | "
        f"evaluating {candidate_total}/{len(candidates)} codon candidates "
        f"(generated in {optimization_elapsed:.1f} s)",
        flush=True,
    )
    feasible: list[AssemblyPlan] = []

    def evaluate(candidates_to_try, strategy: str) -> int:
        total = min(max_candidates, len(candidates_to_try))
        for rank, candidate in enumerate(candidates_to_try[:total], start=1):
            candidate_started = time.perf_counter()
            print(
                f"[{construct_id}] {strategy} candidate {rank}/{total}: "
                f"5'-MFE={candidate.metrics.five_prime_mfe_kcal_mol:.1f} | "
                f"init-open={candidate.metrics.initiation_unpaired_probability:.2f} | "
                f"local-min={candidate.metrics.local_mfe_min_kcal_mol:.1f} | "
                f"SD-like={len(candidate.metrics.sd_like_sites)} | "
                "searching Gibson boundaries, overlaps, and primers...",
                flush=True,
            )
            plasmid = backbone.build_plasmid(candidate.cds)
            plan = plan_assembly(
                design_id=construct_id,
                protein_sequence=design.expression_sequence,
                cds=candidate.cds,
                final_plasmid=plasmid,
                poi_start=backbone.poi_insertion_index,
                intrinsic_score=candidate.combined_score,
                rules=rules,
                primer_design_id=primer_design_id,
            )
            if plan is None:
                print(
                    f"[{construct_id}] {strategy} candidate {rank}/{total}: "
                    f"no feasible plan ({time.perf_counter() - candidate_started:.1f} s)",
                    flush=True,
                )
                continue
            assembly_choice = "preferred_fragment_count"
            if (
                rules["oligos"].get("prefer_short_long_oligos", False)
                and len(plan.synthetic_fragments) < int(rules["assembly"]["max_synthetic_fragments"])
            ):
                current_oligos = _all_fragment_oligos(plan)
                current_long = sum(len(oligo.sequence) > 60 for oligo in current_oligos)
                if current_long > 2:
                    four_fragment_rules = {
                        **rules,
                        "assembly": {
                            **rules["assembly"],
                            "min_synthetic_fragments": 4,
                            "max_synthetic_fragments": 4,
                            "preferred_synthetic_fragments": 4,
                        },
                    }
                    alternative = plan_assembly(
                        design_id=construct_id,
                        protein_sequence=design.expression_sequence,
                        cds=candidate.cds,
                        final_plasmid=plasmid,
                        poi_start=backbone.poi_insertion_index,
                        intrinsic_score=candidate.combined_score,
                        rules=four_fragment_rules,
                        primer_design_id=primer_design_id,
                    )
                    if alternative is not None:
                        alternative_oligos = _all_fragment_oligos(alternative)
                        alternative_long = sum(len(oligo.sequence) > 60 for oligo in alternative_oligos)
                        alternative_flags = sum(bool(oligo.risk_flags) for oligo in alternative_oligos)
                        current_flags = sum(bool(oligo.risk_flags) for oligo in current_oligos)
                        if alternative_flags <= current_flags and alternative_long < current_long:
                            plan = alternative
                            assembly_choice = f"four_fragments_reduce_long_oligos_{current_long}_to_{alternative_long}"
            plan.candidate_rank = rank
            plan.expression_metrics = candidate.metrics.to_dict()
            plan.expression_metrics["codon_strategy"] = strategy
            plan.expression_metrics["assembly_choice"] = assembly_choice
            plan.validation = validate_plan(
                plan,
                design,
                backbone.sequence_without_poi,
                backbone.poi_insertion_index,
                rules["expression"]["fixed_c_terminal_dna"],
                rules,
            )
            feasible.append(plan)
            screening_flags = len(plan.validation["screening_flags"])
            long_oligos = len(plan.validation["long_tier_oligos"])
            print(
                f"[{construct_id}] {strategy} candidate {rank}/{total}: feasible | "
                f"{len(plan.synthetic_fragments)} fragments | "
                f"{4 * len(plan.synthetic_fragments)} design primers | "
                f"long-tier={long_oligos} | flags={screening_flags} | "
                f"{time.perf_counter() - candidate_started:.1f} s",
                flush=True,
            )
        return total

    evaluate(candidates, "default")
    fallback_overrides = rules["codon_optimization"].get("assembly_fallback")
    if not feasible and fallback_overrides:
        fallback_rules = {**rules["codon_optimization"], **fallback_overrides}
        fallback_candidates = rank_expression_candidates(
            optimize_codons(
                design.expression_sequence,
                fallback_rules,
                ncaa_codon=rules["input"]["ncaa_codon"],
                sd_like_motifs=rules["translation_initiation"]["sd_like_motifs"],
                sd_like_search_penalty=fallback_rules["sd_like_search_penalty"],
            ),
            backbone.sequence_without_poi[: backbone.poi_insertion_index],
            rules["translation_initiation"],
        )
        print(f"[{construct_id}] trying GC-balanced assembly fallback", flush=True)
        evaluate(fallback_candidates, "gc_balanced")
    if not feasible:
        raise RuntimeError(
            f"{construct_id}: no feasible plan among {candidate_total} "
            "synonymous CDS candidates"
        )
    plan = min(feasible, key=_plan_selection_key)
    warnings: list[str] = []
    if plan.expression_metrics["codon_strategy"] != "default":
        warnings.append("GC-balanced synonymous codons required for feasible assembly")
    if len(plan.synthetic_fragments) > int(
        rules["assembly"]["preferred_synthetic_fragments"]
    ):
        if plan.expression_metrics["assembly_choice"].startswith("four_fragments_reduce"):
            warnings.append("four fragments selected to reduce long oligonucleotides")
        else:
            warnings.append(
                f"three-fragment plan failed hard constraints; "
                f"used {len(plan.synthetic_fragments)} fragments"
            )
    long_names = plan.validation["long_tier_oligos"]
    if long_names:
        warnings.append("long-tier: " + ", ".join(long_names))
    flagged = plan.validation["screening_flags"]
    if flagged:
        warnings.append("manual structure review: " + ", ".join(flagged))
    expression_flags = plan.expression_metrics["risk_flags"]
    if expression_flags:
        warnings.append("expression review: " + ", ".join(expression_flags))
    plan.warnings = warnings
    expression_metrics = plan.expression_metrics
    print(
        f"[{construct_id}] selected candidate {plan.candidate_rank}/"
        f"{candidate_total} ({expression_metrics['codon_strategy']}) | "
        f"{len(plan.synthetic_fragments)} fragments | "
        f"5'-MFE={expression_metrics['five_prime_mfe_kcal_mol']:.1f} | "
        f"init-open={expression_metrics['initiation_unpaired_probability']:.2f} | "
        f"SD-like={len(expression_metrics['sd_like_sites'])} | "
        f"design elapsed {time.perf_counter() - design_started:.1f} s",
        flush=True,
    )
    return design, plan


def _design_worker(
    job: tuple[
        int,
        Path,
        str,
        str,
        str,
        BackboneModel,
        dict,
        int,
    ],
) -> tuple[int, StructureDesign, AssemblyPlan, float]:
    (
        index,
        structure_path,
        construct_id,
        primer_design_id,
        chain,
        backbone,
        rules,
        max_candidates,
    ) = job
    started = time.perf_counter()
    design, plan = _design_one(
        structure_path,
        construct_id,
        primer_design_id,
        chain,
        backbone,
        rules,
        max_candidates,
    )
    return index, design, plan, time.perf_counter() - started


def _ignore_worker_interrupts() -> None:
    """Let the parent process own Ctrl+C and terminate the pool cleanly."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def deduplicate_primer_rows(
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    unique_rows: list[dict[str, str]] = []
    seen_sequences: set[str] = set()
    for row in rows:
        sequence = "".join(str(row["sequence"]).upper().split())
        if sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)
        unique_rows.append(
            {"name": str(row["name"]), "sequence": sequence}
        )
    return unique_rows


def _usage_prefix(numbers: set[int], all_numbers: set[int]) -> str:
    ordered = sorted(numbers)
    if numbers == all_numbers and len(ordered) >= 5 and ordered == list(range(ordered[0], ordered[-1] + 1)):
        return f"{ordered[0]}to{ordered[-1]}"
    return "_".join(map(str, ordered))


def _primer_rows(
    plans: list[AssemblyPlan], *, include_usage_in_name: bool = False
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    usage: dict[str, set[int]] = {}
    suffixes: dict[str, str] = {}
    numbers: set[int] = set()
    if plans:
        rows.extend(
            {"name": oligo.name, "sequence": oligo.sequence.upper()}
            for oligo in plans[0].vector_oligos
        )
    for plan in plans:
        if include_usage_in_name:
            match = re.search(r"(\d+)$", plan.design_id)
            if match is None:
                raise ValueError(f"Construct ID has no numeric suffix: {plan.design_id}")
            number = int(match.group(1))
            if number in numbers:
                raise ValueError(f"Duplicate construct number for primer naming: {number}")
            numbers.add(number)
            for oligo in plan.vector_oligos:
                usage.setdefault(oligo.sequence.upper(), set()).add(number)
        for fragment in plan.synthetic_fragments:
            for oligo in (
                fragment.outer_forward,
                fragment.outer_reverse,
                fragment.core_forward,
                fragment.core_reverse,
            ):
                rows.append(
                    {"name": oligo.name, "sequence": oligo.sequence.upper()}
                )
                if include_usage_in_name:
                    if not oligo.name.startswith(f"{plan.design_id}-"):
                        raise ValueError(f"Primer name does not match construct ID: {oligo.name}")
                    sequence = oligo.sequence.upper()
                    usage.setdefault(sequence, set()).add(number)
                    suffix = oligo.name[len(plan.design_id) + 1:]
                    if sequence in suffixes and suffixes[sequence] != suffix:
                        raise ValueError(f"Reused primer has different roles: {oligo.name}")
                    suffixes[sequence] = suffix
    unique = deduplicate_primer_rows(rows)
    if include_usage_in_name:
        for row in unique:
            sequence = row["sequence"]
            suffix = row["name"] if row["name"].startswith("vec-") else suffixes[sequence]
            row["name"] = f"{_usage_prefix(usage[sequence], numbers)}-{suffix}"
        names = [row["name"] for row in unique]
        if len(set(names)) != len(names):
            raise ValueError("Shared primer naming produced duplicate names")
    return unique


def _summary_row(
    structure: Path,
    chain: str,
    plan: AssemblyPlan,
    plasmid_prefix: str,
) -> dict[str, object]:
    validation = plan.validation
    expression = plan.expression_metrics
    sd_sites = expression["sd_like_sites"]
    return {
        "map": f"{plasmid_prefix}-{plan.design_id}.dna",
        "structure": structure.name,
        "chain": chain,
        "proteinLength": len(plan.protein_sequence),
        "cdsLength": len(plan.cds),
        "nonstandardResidues": ";".join(
            f"{residue['component_id']}@{residue['expression_position']}"
            for residue in validation["nonstandard_residues"]
        ),
        "fragments": len(plan.synthetic_fragments),
        "designPrimers": 4 * len(plan.synthetic_fragments),
        "longTierPrimers": len(validation["long_tier_oligos"]),
        "oeLengths": ";".join(
            str(value) for value in validation["oe_overlap_lengths_bp"]
        ),
        "oeTm": ";".join(str(value) for value in validation["oe_overlap_tm_c"]),
        "gibsonLengths": ";".join(
            str(value) for value in validation["gibson_overlap_lengths_bp"]
        ),
        "gibsonTm": ";".join(
            str(value) for value in validation["gibson_overlap_tm_c"]
        ),
        "candidate": plan.candidate_rank,
        "fivePrimeMfe": expression["five_prime_mfe_kcal_mol"],
        "initiationOpen": expression["initiation_unpaired_probability"],
        "localMfeMin": expression["local_mfe_min_kcal_mol"],
        "localMfeStart": expression["local_mfe_min_start_nt"],
        "sdLikeSites": ";".join(
            f"{site['motif']}@{site['start_nt']}"
            for site in sd_sites
        ),
        "expressionPenalty": expression["expression_penalty"],
        "expressionFlags": ";".join(expression["risk_flags"]),
        "status": validation["status"],
        "notes": " | ".join(plan.warnings),
    }


def _write_and_validate_design(
    backbone: BackboneModel,
    maps_dir: Path,
    structure_path: Path,
    chain: str,
    design: StructureDesign,
    plan: AssemblyPlan,
    plasmid_prefix: str,
) -> dict[str, object]:
    template = backbone.template_by_poi_length.get(
        len(plan.cds),
        next(iter(backbone.template_by_poi_length.values())),
    )
    map_path = maps_dir / f"{plasmid_prefix}-{plan.design_id}.dna"
    write_clean_snapgene(template, map_path, plan, plasmid_prefix)
    record = read_snapgene(map_path)
    validate_structure_map_pair(design, record)
    expected_primers = 4 + 4 * len(plan.synthetic_fragments)
    observed_primers = sum(
        feature.type == "primer_bind" for feature in record.features
    )
    if observed_primers != expected_primers:
        raise ValueError(
            f"{plan.design_id}: SnapGene round trip found "
            f"{observed_primers} primers; expected {expected_primers}"
        )
    poi = get_poi_feature(record, expected_label=plan.design_id)
    if int(poi.location.start) != backbone.poi_insertion_index:
        raise ValueError(f"{plan.design_id}: POI feature start moved")
    return _summary_row(structure_path, chain, plan, plasmid_prefix)


def run_project(
    primer_build_root: Path,
    project_dir: Path,
    chain: str,
    map_start: str,
    primer_start: str,
    max_candidates: int = 2,
    workers: int = 1,
    overwrite: bool = False,
    rules_path: Path | None = None,
) -> list[AssemblyPlan]:
    project_started = time.perf_counter()
    if max_candidates < 1:
        raise ValueError("--candidates must be at least 1")
    if workers < 1:
        raise ValueError("--workers must be at least 1")
    project_dir = project_dir.resolve()
    golden = (primer_build_root / "ADN round1").resolve()
    if project_dir == golden or golden in project_dir.parents:
        raise ValueError("ADN round1 is a read-only golden reference; copy structures to a new project")
    numbering = parse_numbering(map_start, primer_start)
    rules_path = (rules_path or primer_build_root / "config" / "design_rules.json").resolve()
    rules = load_rules(rules_path)
    structures = discover_structures(project_dir)
    template_name = rules["expression"].get("backbone_template", "pet22b_expression_backbone.dna")
    plasmid_prefix = rules["expression"].get("plasmid_name_prefix", "PET22B")
    if not re.fullmatch(r"[A-Za-z0-9-]+", plasmid_prefix):
        raise ValueError("plasmid_name_prefix must contain only letters, digits, and hyphens")
    if Path(template_name).name != template_name:
        raise ValueError("backbone_template must be a filename within templates/")
    template_path = primer_build_root / "templates" / template_name
    if not template_path.is_file():
        raise FileNotFoundError(
            "Bundled pET22B expression-backbone template is missing: "
            f"{template_path}"
        )
    backbone = reconstruct_backbone(
        [template_path],
        rules["expression"]["fixed_c_terminal_dna"],
    )
    maps_dir = project_dir / "maps"
    primer_workbook = project_dir / "all_primers.xlsx"
    summary_workbook = maps_dir / "statistics_summary.xlsx"
    manifest_path = project_dir / "run_manifest.json"
    targets = [
        maps_dir / f"{plasmid_prefix}-{numbering.construct_id(index)}.dna"
        for index in range(len(structures))
    ] + [primer_workbook, summary_workbook, manifest_path]
    if any(path.resolve() == golden or golden in path.resolve().parents for path in targets):
        raise ValueError("Output paths must not resolve into the read-only golden reference")
    collisions = [path for path in targets if path.exists()]
    if collisions and not overwrite:
        joined = "\n".join(str(path) for path in collisions)
        raise FileExistsError(
            "Output files already exist; rerun with --overwrite to replace "
            f"only these generated targets:\n{joined}"
        )
    effective_workers = min(workers, len(structures))
    mode = "parallel" if effective_workers > 1 else "serial"
    print(
        f"Project: {project_dir.name} | structures={len(structures)} | "
        f"chain={chain} | candidates/design={max_candidates} | "
        f"workers={effective_workers} ({mode})",
        flush=True,
    )
    if workers > effective_workers:
        print(
            f"Requested {workers} workers, limited to {effective_workers} "
            "because there are fewer structures.",
            flush=True,
        )
    print(f"Output directory: {maps_dir}", flush=True)

    jobs = [
        (
            index,
            structure_path,
            numbering.construct_id(index),
            numbering.primer_design_id(index),
            chain,
            backbone,
            rules,
            max_candidates,
        )
        for index, structure_path in enumerate(structures)
    ]
    plan_results: list[AssemblyPlan | None] = [None] * len(structures)
    summary_results: list[dict[str, object] | None] = [None] * len(structures)
    completed = 0

    def accept_result(
        index: int,
        design: StructureDesign,
        plan: AssemblyPlan,
        design_elapsed: float,
    ) -> None:
        nonlocal completed
        structure_path = structures[index]
        summary_results[index] = _write_and_validate_design(
            backbone,
            maps_dir,
            structure_path,
            chain,
            design,
            plan,
            plasmid_prefix,
        )
        plan_results[index] = plan
        completed += 1
        map_name = f"{plasmid_prefix}-{plan.design_id}.dna"
        print(
            f"[completed {completed}/{len(structures)}] [{plan.design_id}] "
            f"PASS | {len(plan.synthetic_fragments)} fragments | "
            f"{4 * len(plan.synthetic_fragments)} design primers | "
            f"map={map_name} | {design_elapsed:.1f} s",
            flush=True,
        )

    if effective_workers == 1:
        for job in jobs:
            index, structure_path, construct_id = job[:3]
            print(
                f"\n[{index + 1}/{len(structures)}] [{construct_id}] "
                f"{structure_path.name} | chain {chain} | starting",
                flush=True,
            )
            result = _design_worker(job)
            accept_result(*result)
    else:
        print(
            f"Starting {effective_workers} worker processes; progress lines "
            "may be interleaved, and every line includes its design ID.",
            flush=True,
        )
        for job in jobs:
            index, structure_path, construct_id = job[:3]
            print(
                f"[{index + 1}/{len(structures)}] [{construct_id}] "
                f"{structure_path.name} | queued",
                flush=True,
            )
        pool = get_context("spawn").Pool(
            processes=effective_workers,
            initializer=_ignore_worker_interrupts,
        )
        try:
            result_iterator = pool.imap_unordered(
                _design_worker,
                jobs,
                chunksize=1,
            )
            while completed < len(jobs):
                try:
                    result = result_iterator.next(timeout=0.5)
                except MultiprocessingTimeoutError:
                    continue
                accept_result(*result)
        except BaseException:
            pool.terminate()
            raise
        else:
            pool.close()
        finally:
            pool.join()

    if any(plan is None for plan in plan_results):
        raise RuntimeError("Internal error: one or more design results are missing")
    if any(summary is None for summary in summary_results):
        raise RuntimeError("Internal error: one or more summary rows are missing")
    plans = [plan for plan in plan_results if plan is not None]
    summary_rows = [
        summary for summary in summary_results if summary is not None
    ]

    primer_rows = _primer_rows(
        plans,
        include_usage_in_name=bool(rules["oligos"].get("order_names_include_construct_usage", False)),
    )
    primer_alias_count = (
        len(plans[0].vector_oligos)
        + sum(4 * len(plan.synthetic_fragments) for plan in plans)
        - len(primer_rows)
    )
    print(
        "\nBuilding primer and statistics workbooks... "
        f"unique oligos={len(primer_rows)} | "
        f"reused sequence aliases removed={primer_alias_count}",
        flush=True,
    )
    build_project_workbooks(
        primer_rows=primer_rows,
        summary_rows=summary_rows,
        primer_output=primer_workbook,
        summary_output=summary_workbook,
    )
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "platform": platform.platform(),
        "architecture": platform.machine(),
        "dependencies": {name: version(name) for name in ("biopython", "ViennaRNA", "openpyxl", "numpy")},
        "parameters": {"chain": chain, "map_start": map_start, "primer_start": primer_start,
                       "candidates": max_candidates, "workers": effective_workers},
        "elapsed_seconds": round(time.perf_counter() - project_started, 3),
        "selected_candidates": {plan.design_id: plan.candidate_rank for plan in plans},
        "selected_codon_strategies": {
            plan.design_id: plan.expression_metrics["codon_strategy"] for plan in plans
        },
        "input_sha256": {str(path.relative_to(project_dir)): hashlib.sha256(path.read_bytes()).hexdigest() for path in structures},
        "rules_file": str(rules_path),
        "rules_sha256": hashlib.sha256(rules_path.read_bytes()).hexdigest(),
        "template_sha256": hashlib.sha256(template_path.read_bytes()).hexdigest(),
        "source_sha256": {str(path.relative_to(primer_build_root)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in sorted((primer_build_root / "src/primer_build").glob("*.py"))},
        "output_sha256": {str(path.relative_to(project_dir)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in targets if path != manifest_path},
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        f"COMPLETE | designs={len(plans)} | maps={maps_dir} | "
        f"primers={primer_workbook.name} | "
        f"total={time.perf_counter() - project_started:.1f} s",
        flush=True,
    )
    return plans
