from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class NonstandardResidue:
    component_id: str
    structure_position: int
    expression_position: int
    parent_component_id: str | None = None


@dataclass(frozen=True)
class StructureDesign:
    design_id: str
    cif_path: Path
    structure_sequence: str
    expression_sequence: str
    ncaa_structure_positions: tuple[int, ...]
    ncaa_expression_positions: tuple[int, ...]
    nonstandard_residues: tuple[NonstandardResidue, ...] = ()
    protein_chain: str = "A"


@dataclass(frozen=True)
class Oligo:
    name: str
    sequence: str
    purpose: str
    orientation: str
    start: int
    end: int
    effective_sequence: str = ""
    effective_tm_c: float | None = None
    full_tm_c: float | None = None
    tier: str = "standard"
    risk_flags: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class SyntheticFragment:
    name: str
    start: int
    end: int
    core_start: int
    core_end: int
    core_forward: Oligo
    core_reverse: Oligo
    outer_forward: Oligo
    outer_reverse: Oligo
    score: float
    risk_flags: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Junction:
    name: str
    start: int
    end: int
    sequence: str
    tm_c: float
    purpose: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class AssemblyPlan:
    design_id: str
    protein_sequence: str
    cds: str
    final_plasmid_sequence: str
    insert_start: int
    insert_end: int
    insert_sequence: str
    synthetic_fragments: list[SyntheticFragment]
    fragment_junctions: list[Junction]
    gibson_junctions: list[Junction]
    vector_oligos: list[Oligo]
    score: float
    candidate_rank: int = 0
    expression_metrics: dict[str, Any] = field(default_factory=dict)
    synonymous_changes: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    validation: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
