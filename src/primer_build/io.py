from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable

from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import PDBParser, is_aa

from .models import NonstandardResidue, StructureDesign


STANDARD_AMINO_ACIDS: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def load_rules(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def _cif_parent_components(data: dict[str, Any]) -> dict[str, str]:
    component_ids = _as_list(data.get("_chem_comp.id"))
    parents = _as_list(data.get("_chem_comp.mon_nstd_parent_comp_id"))
    if len(component_ids) != len(parents):
        return {}
    return {
        component.upper(): parent.upper()
        for component, parent in zip(component_ids, parents)
        if parent not in {"?", ".", ""}
    }


def _target_cif_polymer_entity(
    data: dict[str, Any],
    protein_chain: str,
) -> str | None:
    atom_chains = _as_list(
        data.get("_atom_site.auth_asym_id", data.get("_atom_site.label_asym_id"))
    )
    atom_entities = _as_list(data.get("_atom_site.label_entity_id"))
    polymer_entities = set(_as_list(data.get("_entity_poly_seq.entity_id")))
    candidates = {
        entity
        for chain, entity in zip(atom_chains, atom_entities)
        if chain == protein_chain and entity in polymer_entities
    }

    if not candidates:
        entity_ids = _as_list(data.get("_entity_poly.entity_id"))
        strand_ids = _as_list(data.get("_entity_poly.pdbx_strand_id"))
        for entity, strands in zip(entity_ids, strand_ids):
            aliases = {value.strip() for value in strands.split(",")}
            if protein_chain in aliases and entity in polymer_entities:
                candidates.add(entity)

    if len(candidates) > 1:
        raise ValueError(
            f"chain {protein_chain!r} maps to multiple polymer entities: "
            + ", ".join(sorted(candidates))
        )
    return next(iter(candidates), None)


def _read_cif_polymer_components(
    path: Path,
    protein_chain: str,
    ignored_components: set[str],
    forced_components: set[str],
) -> tuple[list[str], dict[str, str]]:
    data = MMCIF2Dict(str(path))
    parent_components = _cif_parent_components(data)
    entity = _target_cif_polymer_entity(data, protein_chain)
    if entity is not None:
        entities = _as_list(data.get("_entity_poly_seq.entity_id"))
        components = _as_list(data.get("_entity_poly_seq.mon_id"))
        sequence = [
            component.upper()
            for observed_entity, component in zip(entities, components)
            if observed_entity == entity
            and component.upper() not in ignored_components
        ]
        if sequence:
            return sequence, parent_components

    # Fallback for minimal mmCIF files without _entity_poly_seq. Restrict
    # HETATM components to those declared peptide-linking or explicitly forced.
    chains = _as_list(
        data.get("_atom_site.auth_asym_id", data.get("_atom_site.label_asym_id"))
    )
    components = _as_list(data.get("_atom_site.label_comp_id"))
    residue_numbers = _as_list(
        data.get("_atom_site.auth_seq_id", data.get("_atom_site.label_seq_id"))
    )
    insertion_codes = _as_list(data.get("_atom_site.pdbx_PDB_ins_code"))
    if not insertion_codes:
        insertion_codes = ["?"] * len(residue_numbers)
    groups = _as_list(data.get("_atom_site.group_PDB"))
    if not groups:
        groups = ["ATOM"] * len(residue_numbers)
    chem_ids = _as_list(data.get("_chem_comp.id"))
    chem_types = _as_list(data.get("_chem_comp.type"))
    component_types = {
        component.upper(): component_type.lower()
        for component, component_type in zip(chem_ids, chem_types)
    }

    residues: OrderedDict[str, str] = OrderedDict()
    for chain, component, residue_number, insertion_code, group in zip(
        chains,
        components,
        residue_numbers,
        insertion_codes,
        groups,
    ):
        component = component.upper()
        if chain != protein_chain or component in ignored_components:
            continue
        polymer_like = (
            component in STANDARD_AMINO_ACIDS
            or component in forced_components
            or "peptide linking" in component_types.get(component, "")
            or group.upper() == "ATOM"
        )
        if polymer_like:
            residues.setdefault(f"{residue_number}:{insertion_code}", component)
    return list(residues.values()), parent_components


def _read_pdb_seqres_components(
    path: Path,
    protein_chain: str,
    ignored_components: set[str],
) -> list[str]:
    components: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("SEQRES"):
                continue
            chain = line[11].strip()
            if chain == protein_chain:
                components.extend(
                    component.upper()
                    for component in line[19:70].split()
                    if component.upper() not in ignored_components
                )
    return components


def _peptide_bonded(left, right) -> bool:
    if "C" not in left or "N" not in right:
        return False
    return float(left["C"] - right["N"]) <= 2.2


def _read_pdb_polymer_components(
    path: Path,
    protein_chain: str,
    ignored_components: set[str],
    forced_components: set[str],
) -> list[str]:
    seqres = _read_pdb_seqres_components(
        path,
        protein_chain,
        ignored_components,
    )
    if seqres:
        return seqres

    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    model = next(structure.get_models())
    if protein_chain not in model:
        raise ValueError(
            f"{path.name}: protein chain {protein_chain!r} not found"
        )
    residues = sorted(
        (
            residue
            for residue in model[protein_chain].get_residues()
            if residue.get_resname().strip().upper() not in ignored_components
        ),
        key=lambda residue: (residue.id[1], residue.id[2], residue.id[0]),
    )
    components: list[str] = []
    for index, residue in enumerate(residues):
        component = residue.get_resname().strip().upper()
        hetero_flag = str(residue.id[0]).strip()
        previous_bond = index > 0 and _peptide_bonded(
            residues[index - 1],
            residue,
        )
        next_bond = index + 1 < len(residues) and _peptide_bonded(
            residue,
            residues[index + 1],
        )
        polymer_like = (
            not hetero_flag
            or is_aa(residue, standard=False)
            or component in forced_components
            or previous_bond
            or next_bond
        )
        if polymer_like:
            components.append(component)
    return components


def read_structure_design(
    cif_path: Path,
    design_id: str,
    protein_chain: str = "A",
    prepend_initiator_methionine: bool = True,
    ignored_comp_ids: Iterable[str] = ("HOH", "WAT"),
    forced_polymer_comp_ids: Iterable[str] = (),
    component_translation_overrides: dict[str, str] | None = None,
) -> StructureDesign:
    suffix = cif_path.suffix.lower()
    ignored_components = {
        component.upper() for component in ignored_comp_ids
    }
    forced_components = {
        component.upper() for component in forced_polymer_comp_ids
    }
    overrides = {
        component.upper(): one_letter.upper()
        for component, one_letter in (
            component_translation_overrides or {}
        ).items()
    }
    invalid_overrides = {
        component: one_letter
        for component, one_letter in overrides.items()
        if len(one_letter) != 1
        or one_letter not in set(STANDARD_AMINO_ACIDS.values()) | {"X"}
    }
    if invalid_overrides:
        raise ValueError(
            "component_translation_overrides values must be one of the "
            "20 standard amino-acid letters or X: "
            + ", ".join(
                f"{component}={value}"
                for component, value in sorted(invalid_overrides.items())
            )
        )

    parent_components: dict[str, str] = {}
    if suffix in {".cif", ".mmcif"}:
        components, parent_components = _read_cif_polymer_components(
            cif_path,
            protein_chain,
            ignored_components,
            forced_components,
        )
    elif suffix == ".pdb":
        components = _read_pdb_polymer_components(
            cif_path,
            protein_chain,
            ignored_components,
            forced_components,
        )
    else:
        raise ValueError(
            f"{cif_path.name}: unsupported structure format {cif_path.suffix!r}; "
            "expected .pdb, .cif, or .mmcif"
        )
    if not components:
        raise ValueError(f"{cif_path.name}: protein chain {protein_chain!r} not found")

    sequence_parts: list[str] = []
    ncaa_positions: list[int] = []
    nonstandard_components: list[tuple[str, int, str | None]] = []
    for position, component in enumerate(components, start=1):
        component = component.upper()
        if component in overrides:
            one_letter = overrides[component]
        else:
            one_letter = STANDARD_AMINO_ACIDS.get(component, "X")
        sequence_parts.append(one_letter)
        if one_letter == "X":
            ncaa_positions.append(position)
            nonstandard_components.append(
                (
                    component,
                    position,
                    parent_components.get(component),
                )
            )

    structure_sequence = "".join(sequence_parts)
    add_met = prepend_initiator_methionine and not structure_sequence.startswith("M")
    expression_sequence = ("M" if add_met else "") + structure_sequence
    offset = 1 if add_met else 0
    expression_positions = tuple(position + offset for position in ncaa_positions)
    nonstandard_residues = tuple(
        NonstandardResidue(
            component_id=component,
            structure_position=position,
            expression_position=position + offset,
            parent_component_id=parent,
        )
        for component, position, parent in nonstandard_components
    )
    return StructureDesign(
        design_id=design_id,
        cif_path=cif_path,
        structure_sequence=structure_sequence,
        expression_sequence=expression_sequence,
        ncaa_structure_positions=tuple(ncaa_positions),
        ncaa_expression_positions=expression_positions,
        nonstandard_residues=nonstandard_residues,
        protein_chain=protein_chain,
    )


def read_snapgene(path: Path):
    return SeqIO.read(path, "snapgene")


def get_poi_feature(record, expected_label: str | None = None):
    cds_features = [
        feature
        for feature in record.features
        if feature.type == "CDS"
    ]
    if expected_label is not None:
        matches = [
            feature
            for feature in cds_features
            if expected_label
            in {str(value) for value in feature.qualifiers.get("label", [])}
        ]
    else:
        matches = [
            feature
            for feature in cds_features
            if any(
                (
                    "intentional amber codon" in str(value).lower()
                    or "internal tag codon" in str(value).lower()
                )
                for value in feature.qualifiers.get("note", [])
            )
        ]
    if not matches:
        # Backward-compatible fallback for a renamed construct without an
        # explicit note: the POI is the substantial CDS beginning at the
        # validated backbone insertion coordinate 5200.
        matches = [
            feature
            for feature in cds_features
            if int(feature.location.start) == 5200
            and int(feature.location.end) - int(feature.location.start) >= 60
        ]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one POI CDS feature, found {len(matches)}")
    return matches[0]


def map_expression_sequence(record, expected_label: str | None = None) -> str:
    feature = get_poi_feature(record, expected_label=expected_label)
    return str(feature.extract(record.seq).translate()).replace("*", "X")


def validate_structure_map_pair(design: StructureDesign, record) -> None:
    observed = map_expression_sequence(record, expected_label=design.design_id)
    if observed != design.expression_sequence:
        raise ValueError(
            f"{design.design_id}: CIF/map protein mismatch\n"
            f"CIF: {design.expression_sequence}\nMAP: {observed}"
        )


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, default=str)
