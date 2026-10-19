from __future__ import annotations

import copy
import struct
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

from .models import AssemblyPlan, Oligo


COOKIE = 0x09
DNA = 0x00
PRIMERS = 0x05
NOTES = 0x06
FEATURES = 0x0A


def read_packets(path: Path) -> list[tuple[int, bytes]]:
    data = path.read_bytes()
    packets: list[tuple[int, bytes]] = []
    offset = 0
    while offset < len(data):
        if offset + 5 > len(data):
            raise ValueError(f"{path.name}: truncated SnapGene packet header")
        packet_type = data[offset]
        length = struct.unpack(">I", data[offset + 1 : offset + 5])[0]
        start = offset + 5
        end = start + length
        if end > len(data):
            raise ValueError(f"{path.name}: truncated SnapGene packet payload")
        packets.append((packet_type, data[start:end]))
        offset = end
    if not packets or packets[0][0] != COOKIE:
        raise ValueError(f"{path.name}: missing SnapGene cookie")
    return packets


def _pack(packet_type: int, payload: bytes) -> bytes:
    return bytes([packet_type]) + struct.pack(">I", len(payload)) + payload


def _first_packet(packets: list[tuple[int, bytes]], packet_type: int) -> bytes:
    for observed_type, payload in packets:
        if observed_type == packet_type:
            return payload
    raise ValueError(f"Template is missing packet type 0x{packet_type:02x}")


def _replace_text_qualifier(feature: ET.Element, name: str, text: str) -> None:
    for qualifier in feature.findall("Q"):
        if qualifier.get("name") == name:
            for child in list(qualifier):
                qualifier.remove(child)
            ET.SubElement(qualifier, "V", {"text": text})
            return
    qualifier = ET.SubElement(feature, "Q", {"name": name})
    ET.SubElement(qualifier, "V", {"text": text})


def _features_payload(
    template_payload: bytes,
    design_id: str,
    protein_sequence: str,
    cds_length: int,
) -> bytes:
    root = ET.fromstring(template_payload.decode("utf-8"))
    cds_features = [
        feature
        for feature in root.findall("Feature")
        if feature.get("type") == "CDS"
    ]
    linker_features = [
        feature
        for feature in cds_features
        if (feature.get("name") or "").strip().lower() == "gs linker"
    ]
    if len(linker_features) != 1:
        raise ValueError(
            "Expected exactly one GS linker CDS in the backbone template, "
            f"found {len(linker_features)}"
        )
    linker_segments = linker_features[0].findall("Segment")
    if len(linker_segments) != 1:
        raise ValueError("Template GS linker CDS must contain exactly one segment")
    linker_start_1 = int(
        linker_segments[0].get("range", "").split("-")[0]
    )
    poi_features = []
    for feature in cds_features:
        segments = feature.findall("Segment")
        if len(segments) != 1:
            continue
        _, feature_end_1 = (
            int(value)
            for value in segments[0].get("range", "").split("-")
        )
        if feature_end_1 + 1 == linker_start_1:
            poi_features.append(feature)
    if len(poi_features) != 1:
        raise ValueError(
            "Expected exactly one POI CDS immediately upstream of the GS "
            f"linker, found {len(poi_features)}"
        )
    poi = poi_features[0]
    poi_segments = poi.findall("Segment")
    if len(poi_segments) != 1:
        raise ValueError("Template POI CDS must contain exactly one segment")
    old_start_1, old_end_1 = (
        int(value) for value in poi_segments[0].get("range", "").split("-")
    )
    old_length = old_end_1 - old_start_1 + 1
    delta = cds_length - old_length

    for feature in root.findall("Feature"):
        for attribute in ("maxRunOn", "maxFusedRunOn"):
            if feature.get(attribute) is not None:
                value = int(feature.get(attribute, "0"))
                if value > old_end_1:
                    feature.set(attribute, str(value + delta))
        for segment in feature.findall("Segment"):
            range_spec = segment.get("range")
            if not range_spec:
                continue
            segment_start, segment_end = (
                int(value) for value in range_spec.split("-")
            )
            if feature is poi:
                segment.set(
                    "range",
                    f"{old_start_1}-{old_start_1 + cds_length - 1}",
                )
            elif segment_start > old_end_1:
                segment.set(
                    "range",
                    f"{segment_start + delta}-{segment_end + delta}",
                )
            elif segment_end >= old_start_1:
                raise ValueError(
                    "Template contains a non-POI feature overlapping the POI CDS; "
                    "automatic coordinate shifting would be ambiguous"
                )

    poi.set("name", design_id)
    _replace_text_qualifier(poi, "translation", protein_sequence.replace("X", "*"))
    _replace_text_qualifier(
        poi,
        "note",
        "Internal TAG codon(s) intentionally represent nonstandard "
        "structure residue(s) for ncAA/PTM incorporation.",
    )
    return (
        b'<?xml version="1.0"?>'
        + ET.tostring(root, encoding="utf-8", short_empty_elements=True)
    )


def _notes_payload(template_payload: bytes, design_id: str, plasmid_prefix: str) -> bytes:
    root = ET.fromstring(template_payload.decode("utf-8"))

    def set_child(tag: str, value: str) -> None:
        child = root.find(tag)
        if child is None:
            child = ET.SubElement(root, tag)
        child.text = value

    today = datetime.now(timezone.utc).strftime("%Y.%m.%d")
    set_child("UUID", str(uuid.uuid4()))
    set_child("Type", "Synthetic")
    set_child("ConfirmedExperimentally", "0")
    set_child("LastModified", today)
    set_child(
        "Comments",
        f"{plasmid_prefix}-{design_id}; automated oligo assembly, overlap extension, and Gibson plan.",
    )
    return ET.tostring(root, encoding="utf-8", short_empty_elements=True)


def _all_plan_oligos(plan: AssemblyPlan) -> list[tuple[Oligo, int]]:
    oligos: list[tuple[Oligo, int]] = [(oligo, 0) for oligo in plan.vector_oligos]
    for fragment in plan.synthetic_fragments:
        for oligo in (
            fragment.core_forward,
            fragment.core_reverse,
            fragment.outer_forward,
            fragment.outer_reverse,
        ):
            oligos.append((oligo, plan.insert_start))
    return oligos


def _primers_payload(plan: AssemblyPlan) -> bytes:
    dated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    oligos = _all_plan_oligos(plan)
    root = ET.Element("Primers", {"nextValidID": str(len(oligos))})
    ET.SubElement(
        root,
        "HybridizationParams",
        {
            "minContinuousMatchLen": "10",
            "allowMismatch": "1",
            "minMeltingTemperature": "40",
            "showAdditionalFivePrimeMatches": "1",
            "minimumFivePrimeAnnealing": "15",
        },
    )
    for recent_id, (oligo, offset) in enumerate(oligos):
        start = oligo.start + offset
        end = oligo.end + offset
        location = f"{start}-{end - 1}"
        full_tm = oligo.full_tm_c if oligo.full_tm_c is not None else oligo.effective_tm_c
        description = (
            f"purpose={oligo.purpose}; effective_Tm={oligo.effective_tm_c:.1f} C; "
            f"tier={oligo.tier}"
            if oligo.effective_tm_c is not None
            else f"purpose={oligo.purpose}; tier={oligo.tier}"
        )
        if oligo.risk_flags:
            description += "; flags=" + ",".join(oligo.risk_flags)
        primer = ET.SubElement(
            root,
            "Primer",
            {
                "recentID": str(recent_id),
                "name": oligo.name,
                "sequence": oligo.sequence,
                "description": f"<html><body>{description}</body></html>",
                "dateAdded": dated,
            },
        )
        strand = "0" if oligo.orientation == "FW" else "1"
        site_attributes = {
            "location": location,
            "boundStrand": strand,
            "annealedBases": oligo.sequence,
            "meltingTemperature": str(int(round(full_tm or 0.0))),
        }
        site = ET.SubElement(primer, "BindingSite", site_attributes)
        ET.SubElement(
            site,
            "Component",
            {"hybridizedRange": location, "bases": oligo.sequence},
        )
        simplified = copy.deepcopy(site)
        simplified.set("simplified", "1")
        primer.append(simplified)
    return (
        b'<?xml version="1.0"?>'
        + ET.tostring(root, encoding="utf-8", short_empty_elements=True)
    )


def write_clean_snapgene(
    template_path: Path,
    output_path: Path,
    plan: AssemblyPlan,
    plasmid_prefix: str = "PET22B",
) -> None:
    """Write a compact SnapGene file without the template's edit history."""
    packets = read_packets(template_path)
    cookie = _first_packet(packets, COOKIE)
    template_dna = _first_packet(packets, DNA)
    flags = template_dna[:1]
    features = _features_payload(
        _first_packet(packets, FEATURES),
        plan.design_id,
        plan.protein_sequence,
        len(plan.cds),
    )
    notes = _notes_payload(_first_packet(packets, NOTES), plan.design_id, plasmid_prefix)
    primers = _primers_payload(plan)
    payloads = [
        (COOKIE, cookie),
        (DNA, flags + plan.final_plasmid_sequence.lower().encode("ascii")),
        (FEATURES, features),
        (PRIMERS, primers),
        (NOTES, notes),
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(b"".join(_pack(packet_type, data) for packet_type, data in payloads))
