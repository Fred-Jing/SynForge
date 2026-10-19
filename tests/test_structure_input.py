from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from primer_build.io import read_structure_design  # noqa: E402
from primer_build.project_pipeline import _usage_prefix, deduplicate_primer_rows, parse_numbering  # noqa: E402


class StructureInputTests(unittest.TestCase):
    def test_map_numbering_accepts_both_hyphenated_and_plain_ids(self) -> None:
        self.assertEqual(parse_numbering('ADN-1', 'A1-1-1-FW').construct_id(2), 'ADN-3')
        self.assertEqual(parse_numbering('ARS1', 'ARS1-1-1-FW').construct_id(11), 'ARS12')

    def test_pdb_hetero_residue_uses_residue_number_not_record_order(self) -> None:
        def atom(serial: int, name: str, component: str, position: int, x: float, hetero: bool = False) -> str:
            group = "HETATM" if hetero else "ATOM  "
            return (
                f"{group}{serial:5d} {name:<4s} {component:>3s} A{position:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {name[0]:>2s}  \n"
            )

        pdb_text = "".join(
            [
                atom(1, "N", "ALA", 1, 0.0), atom(2, "C", "ALA", 1, 1.3),
                atom(3, "N", "GLY", 3, 4.0), atom(4, "C", "GLY", 3, 5.3),
                atom(5, "N", "7N8", 2, 2.0, True), atom(6, "C", "7N8", 2, 3.3, True),
                "END\n",
            ]
        )
        with tempfile.TemporaryDirectory() as temp_name:
            path = Path(temp_name) / "out_of_order.pdb"
            path.write_text(pdb_text, encoding="utf-8")
            design = read_structure_design(path, "example", protein_chain="A")

        self.assertEqual(design.structure_sequence, "AXG")
        self.assertEqual(design.ncaa_expression_positions, (3,))

    def test_cif_polymer_nonstandard_residues_are_generic(self) -> None:
        cif_text = """\
data_generic
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 ALA
1 2 7N8
1 3 SEP
1 4 GLY
loop_
_atom_site.auth_asym_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_comp_id
_atom_site.auth_seq_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
A A 1 ALA 1 1 ?
A A 1 7N8 2 2 ?
A A 1 SEP 3 3 ?
A A 1 GLY 4 4 ?
loop_
_chem_comp.id
_chem_comp.mon_nstd_parent_comp_id
ALA ?
7N8 PHE
SEP SER
GLY ?
"""
        with tempfile.TemporaryDirectory() as temp_name:
            path = Path(temp_name) / "generic.cif"
            path.write_text(cif_text, encoding="utf-8")
            design = read_structure_design(path, "generic", protein_chain="A")

        self.assertEqual(design.structure_sequence, "AXXG")
        self.assertEqual(design.expression_sequence, "MAXXG")
        self.assertEqual(design.ncaa_structure_positions, (2, 3))
        self.assertEqual(design.ncaa_expression_positions, (3, 4))
        self.assertEqual(
            [
                (
                    residue.component_id,
                    residue.structure_position,
                    residue.expression_position,
                    residue.parent_component_id,
                )
                for residue in design.nonstandard_residues
            ],
            [
                ("7N8", 2, 3, "PHE"),
                ("SEP", 3, 4, "SER"),
            ],
        )

    def test_pdb_seqres_accepts_arbitrary_modified_residue_names(self) -> None:
        residues = "ALA SEP NLE GLY"
        pdb_text = f"SEQRES {1:>3} A{4:>5}  {residues}\nEND\n"
        with tempfile.TemporaryDirectory() as temp_name:
            path = Path(temp_name) / "generic.pdb"
            path.write_text(pdb_text, encoding="utf-8")
            design = read_structure_design(path, "generic", protein_chain="A")

        self.assertEqual(design.structure_sequence, "AXXG")
        self.assertEqual(
            [residue.component_id for residue in design.nonstandard_residues],
            ["SEP", "NLE"],
        )

    def test_component_override_can_map_a_modified_code_to_standard(self) -> None:
        residues = "MET MSE GLY"
        pdb_text = f"SEQRES {1:>3} A{3:>5}  {residues}\nEND\n"
        with tempfile.TemporaryDirectory() as temp_name:
            path = Path(temp_name) / "override.pdb"
            path.write_text(pdb_text, encoding="utf-8")
            design = read_structure_design(
                path,
                "override",
                protein_chain="A",
                component_translation_overrides={"MSE": "M"},
            )

        self.assertEqual(design.structure_sequence, "MMG")
        self.assertEqual(design.nonstandard_residues, ())

    def test_primer_order_rows_are_deduplicated_by_sequence(self) -> None:
        rows = [
            {"name": "A1-1-1-FW", "sequence": "acgt acgt"},
            {"name": "A2-1-1-FW", "sequence": "ACGTACGT"},
            {"name": "A1-1-RV", "sequence": "TTGGCCAA"},
        ]
        observed = deduplicate_primer_rows(rows)
        self.assertEqual(
            observed,
            [
                {"name": "A1-1-1-FW", "sequence": "ACGTACGT"},
                {"name": "A1-1-RV", "sequence": "TTGGCCAA"},
            ],
        )

    def test_shared_primer_name_shows_using_constructs_compactly(self) -> None:
        all_numbers = set(range(1, 13))
        self.assertEqual(_usage_prefix({1, 2, 4}, all_numbers), "1_2_4")
        self.assertEqual(_usage_prefix(all_numbers, all_numbers), "1to12")
        self.assertEqual(_usage_prefix({5, 6, 7, 8}, all_numbers), "5_6_7_8")


if __name__ == "__main__":
    unittest.main()
