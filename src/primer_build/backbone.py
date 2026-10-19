from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq

from .io import get_poi_feature, read_snapgene


@dataclass(frozen=True)
class BackboneModel:
    sequence_without_poi: str
    poi_insertion_index: int
    source_maps: tuple[Path, ...]
    template_by_poi_length: dict[int, Path]

    def build_plasmid(self, cds: str) -> str:
        return (
            self.sequence_without_poi[: self.poi_insertion_index]
            + cds.upper()
            + self.sequence_without_poi[self.poi_insertion_index :]
        )


def reconstruct_backbone(map_paths: list[Path], fixed_c_terminal_dna: str) -> BackboneModel:
    if not map_paths:
        raise ValueError("No SnapGene maps supplied")
    stripped_sequences: list[str] = []
    insertion_indices: list[int] = []
    templates: dict[int, Path] = {}

    for path in map_paths:
        record = read_snapgene(path)
        poi = get_poi_feature(record)
        start = int(poi.location.start)
        end = int(poi.location.end)
        sequence = str(record.seq).upper()
        if not sequence[end:].startswith(fixed_c_terminal_dna.upper()):
            observed = sequence[end : end + len(fixed_c_terminal_dna)]
            raise ValueError(
                f"{path.name}: unexpected C-terminal cassette {observed}; "
                f"expected {fixed_c_terminal_dna}"
            )
        stripped_sequences.append(sequence[:start] + sequence[end:])
        insertion_indices.append(start)
        templates[end - start] = path

    if len(set(insertion_indices)) != 1:
        raise ValueError(f"POI insertion indices disagree: {insertion_indices}")
    first = stripped_sequences[0]
    disagreements = [
        map_paths[index].name
        for index, sequence in enumerate(stripped_sequences)
        if sequence != first
    ]
    if disagreements:
        raise ValueError(f"Backbone sequences disagree outside POI: {disagreements}")
    return BackboneModel(
        sequence_without_poi=first,
        poi_insertion_index=insertion_indices[0],
        source_maps=tuple(map_paths),
        template_by_poi_length=templates,
    )


def translate_expression_cassette(cds: str, fixed_c_terminal_dna: str) -> str:
    sequence = cds.upper() + fixed_c_terminal_dna.upper()
    if len(sequence) % 3:
        raise ValueError("Expression cassette length is not divisible by 3")
    return str(Seq(sequence).translate())
