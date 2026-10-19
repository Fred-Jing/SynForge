from __future__ import annotations

from Bio.Seq import Seq

from .models import AssemblyPlan, Oligo, StructureDesign
from .thermo import reverse_complement


def _translate_with_ncaa(cds: str, ncaa_codon: str) -> str:
    protein: list[str] = []
    for offset in range(0, len(cds), 3):
        codon = cds[offset : offset + 3]
        if codon == ncaa_codon:
            protein.append("X")
        else:
            amino_acid = str(Seq(codon).translate())
            if amino_acid == "*":
                raise ValueError(f"Unexpected stop codon {codon} at codon {offset // 3 + 1}")
            protein.append(amino_acid)
    return "".join(protein)


def _validate_oligo_mapping(
    oligo: Oligo,
    reference: str,
    coordinate_offset: int = 0,
) -> None:
    observed = reference[
        coordinate_offset + oligo.start : coordinate_offset + oligo.end
    ]
    expected = observed if oligo.orientation == "FW" else reverse_complement(observed)
    if oligo.sequence != expected:
        raise ValueError(f"{oligo.name}: oligo does not match its declared interval")


def validate_plan(
    plan: AssemblyPlan,
    design: StructureDesign,
    sequence_without_poi: str,
    poi_insertion_index: int,
    fixed_c_terminal_dna: str,
    rules: dict,
) -> dict[str, object]:
    ncaa_codon = rules["input"]["ncaa_codon"]
    if len(plan.cds) % 3:
        raise ValueError(f"{plan.design_id}: CDS length is not divisible by three")
    translated = _translate_with_ncaa(plan.cds, ncaa_codon)
    if translated != design.expression_sequence:
        raise ValueError(
            f"{plan.design_id}: CDS translation mismatch\n"
            f"expected {design.expression_sequence}\nobserved {translated}"
        )
    in_frame_ncaa_count = sum(
        plan.cds[offset : offset + 3] == ncaa_codon
        for offset in range(0, len(plan.cds), 3)
    )
    if in_frame_ncaa_count != len(design.ncaa_expression_positions):
        raise ValueError(f"{plan.design_id}: unexpected TAG count")
    for position in design.ncaa_expression_positions:
        if plan.cds[(position - 1) * 3 : position * 3] != ncaa_codon:
            raise ValueError(f"{plan.design_id}: missing TAG at codon {position}")

    expected_tail_translation = rules["expression"]["fixed_c_terminal_translation"]
    observed_tail_translation = str(Seq(fixed_c_terminal_dna).translate())
    if observed_tail_translation != expected_tail_translation:
        raise ValueError(
            f"{plan.design_id}: C-terminal DNA translates to {observed_tail_translation}, "
            f"expected {expected_tail_translation}"
        )

    expected_plasmid = (
        sequence_without_poi[:poi_insertion_index]
        + plan.cds
        + sequence_without_poi[poi_insertion_index:]
    )
    if plan.final_plasmid_sequence != expected_plasmid:
        raise ValueError(f"{plan.design_id}: final plasmid does not match backbone + CDS")
    if not plan.final_plasmid_sequence[
        poi_insertion_index + len(plan.cds) :
    ].startswith(fixed_c_terminal_dna):
        raise ValueError(f"{plan.design_id}: configured C-terminal cassette was altered")
    if plan.insert_sequence != plan.final_plasmid_sequence[
        plan.insert_start : plan.insert_end
    ]:
        raise ValueError(f"{plan.design_id}: insert coordinate mismatch")

    fragment_sequences: list[str] = []
    for fragment in plan.synthetic_fragments:
        sequence = plan.insert_sequence[fragment.start : fragment.end]
        fragment_sequences.append(sequence)
        for oligo in (
            fragment.core_forward,
            fragment.core_reverse,
            fragment.outer_forward,
            fragment.outer_reverse,
        ):
            _validate_oligo_mapping(oligo, plan.insert_sequence)
            if len(oligo.sequence) > int(rules["oligos"]["absolute_max_nt"]):
                raise ValueError(f"{oligo.name}: exceeds absolute oligo length limit")
        for outer in (fragment.outer_forward, fragment.outer_reverse):
            if len(outer.sequence) > int(rules["oligos"]["standard_tier_max_nt"]):
                raise ValueError(f"{outer.name}: outer primer exceeds 59 nt")
            if not (
                float(rules["thermodynamics"]["annealing_tm_min_c"])
                <= float(outer.effective_tm_c)
                <= float(rules["thermodynamics"]["annealing_tm_max_c"])
            ):
                raise ValueError(f"{outer.name}: effective Tm outside 62–66 C")

    rebuilt = fragment_sequences[0]
    for previous, current, junction in zip(
        plan.synthetic_fragments,
        plan.synthetic_fragments[1:],
        plan.fragment_junctions,
    ):
        if (junction.start, junction.end) != (current.start, previous.end):
            raise ValueError(f"{junction.name}: inconsistent fragment overlap coordinates")
        if plan.insert_sequence[junction.start : junction.end] != junction.sequence:
            raise ValueError(f"{junction.name}: overlap sequence mismatch")
        rebuilt += plan.insert_sequence[previous.end : current.end]
    if rebuilt != plan.insert_sequence:
        raise ValueError(f"{plan.design_id}: OE fragment reconstruction failed")

    for junction in plan.fragment_junctions:
        if not 20 <= junction.length <= 35:
            raise ValueError(f"{junction.name}: OE overlap length outside 20–35 bp")
        if not 62.0 <= junction.tm_c <= 66.0:
            raise ValueError(f"{junction.name}: OE overlap Tm outside 62–66 C")
    for junction in plan.gibson_junctions:
        if not 20 <= junction.length <= 35:
            raise ValueError(f"{junction.name}: Gibson overlap length outside 20–35 bp")
        if junction.tm_c < 60.0:
            raise ValueError(f"{junction.name}: Gibson overlap Tm below 60 C")
        motif = plan.final_plasmid_sequence[junction.start : junction.end]
        if motif != junction.sequence:
            raise ValueError(f"{junction.name}: Gibson coordinate mismatch")

    for oligo in plan.vector_oligos:
        _validate_oligo_mapping(oligo, plan.final_plasmid_sequence)
        if len(oligo.sequence) > int(rules["oligos"]["standard_tier_max_nt"]):
            raise ValueError(f"{oligo.name}: vector primer exceeds 59 nt")

    left, right, split = plan.gibson_junctions
    vector_amplicon_1_bp = left.end - split.start
    vector_amplicon_2_bp = (
        len(plan.final_plasmid_sequence) - right.start + split.end
    )
    amplicon_min = int(rules["assembly"]["vector_amplicon_min_bp"])
    amplicon_max = int(rules["assembly"]["vector_amplicon_max_bp"])
    if not (
        amplicon_min <= vector_amplicon_1_bp <= amplicon_max
        and amplicon_min <= vector_amplicon_2_bp <= amplicon_max
    ):
        raise ValueError(f"{plan.design_id}: vector amplicon length outside limits")
    assembled_length = (
        len(plan.insert_sequence)
        + vector_amplicon_1_bp
        + vector_amplicon_2_bp
        - sum(junction.length for junction in plan.gibson_junctions)
    )
    if assembled_length != len(plan.final_plasmid_sequence):
        raise ValueError(f"{plan.design_id}: Gibson length accounting failed")

    long_oligos = [
        oligo.name
        for fragment in plan.synthetic_fragments
        for oligo in (fragment.core_forward, fragment.core_reverse)
        if oligo.tier == "long"
    ]
    flagged_oligos = sorted(
        {
            oligo.name
            for fragment in plan.synthetic_fragments
            for oligo in (
                fragment.core_forward,
                fragment.core_reverse,
                fragment.outer_forward,
                fragment.outer_reverse,
            )
            if oligo.risk_flags
        }
        | {oligo.name for oligo in plan.vector_oligos if oligo.risk_flags}
    )
    return {
        "status": "PASS",
        "plasmid_length_bp": len(plan.final_plasmid_sequence),
        "poi_length_bp": len(plan.cds),
        "protein_length_aa": len(design.expression_sequence),
        "ncaa_codon": ncaa_codon,
        "ncaa_expression_positions": list(design.ncaa_expression_positions),
        "nonstandard_residues": [
            {
                "component_id": residue.component_id,
                "structure_position": residue.structure_position,
                "expression_position": residue.expression_position,
                "parent_component_id": residue.parent_component_id,
            }
            for residue in design.nonstandard_residues
        ],
        "synthetic_fragment_count": len(plan.synthetic_fragments),
        "oe_overlap_lengths_bp": [j.length for j in plan.fragment_junctions],
        "oe_overlap_tm_c": [round(j.tm_c, 2) for j in plan.fragment_junctions],
        "gibson_overlap_lengths_bp": [j.length for j in plan.gibson_junctions],
        "gibson_overlap_tm_c": [round(j.tm_c, 2) for j in plan.gibson_junctions],
        "vector_amplicon_lengths_bp": [
            vector_amplicon_1_bp,
            vector_amplicon_2_bp,
        ],
        "long_tier_oligos": long_oligos,
        "screening_flags": flagged_oligos,
    }
