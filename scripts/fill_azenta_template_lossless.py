"""Fill Azenta's original XLSX without rewriting its vendor-specific OOXML.

The vendor upload rejects a workbook round-tripped through a generic writer.
Only sheet1 cell contents and the shared-string table change here; hidden
sheets, x14 validations, drawings, metadata, relationships and calc chain are
copied directly from the original template.
"""
from __future__ import annotations

import argparse
from datetime import date
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape
from zipfile import ZipFile

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
MAIN = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def cell_pattern(ref):
    start = rf'<c r="{re.escape(ref)}"'
    return re.compile(rf'{start}[^>]*?/>|{start}[^>]*?>.*?</c>', re.DOTALL)


def replace_cell(row_xml, ref, value, kind, string_index):
    match = cell_pattern(ref).search(row_xml)
    if match is None:
        raise ValueError(f'Missing existing styled cell {ref}')
    old = match.group(0)
    style = re.search(r'\bs="(\d+)"', old)
    if style is None:
        raise ValueError(f'Missing original cell style {ref}')
    prefix = f'<c r="{ref}" s="{style.group(1)}"'
    if kind == 'string':
        replacement = f'{prefix} t="s"><v>{string_index(value)}</v></c>'
    elif kind == 'number':
        replacement = f'{prefix}><v>{value}</v></c>'
    elif kind == 'formula_cache':
        if '<f' not in old or '<v>0</v>' not in old:
            raise ValueError(f'Unexpected formula cell {ref}')
        replacement = old.replace('<v>0</v>', f'<v>{value}</v>', 1)
    else:
        raise ValueError(kind)
    return row_xml[:match.start()] + replacement + row_xml[match.end():]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--orders', type=Path, required=True)
    parser.add_argument('--contact-source', type=Path, required=True, help='Workbook with contact fields B5:B12 on 订购单; never published')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--order-date', default=date.today().isoformat())
    parser.add_argument('--expected-count', type=int)
    parser.add_argument('--max-length', type=int, default=89)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Output already exists; choose a new file')
    orders = json.loads(args.orders.read_text())
    if not orders:
        raise ValueError('Order list must not be empty')
    if args.expected_count is not None and len(orders) != args.expected_count:
        raise ValueError(f'Expected {args.expected_count} primers, found {len(orders)}')
    if args.output.resolve() in {args.template.resolve(), args.contact_source.resolve()}:
        raise ValueError('Output must not overwrite the original template or contact source')
    date.fromisoformat(args.order_date)
    prior = load_workbook(args.contact_source, read_only=True, data_only=True)['订购单']
    contact = {
        3: args.order_date,
        4: prior['B5'].value,
        5: prior['B6'].value,
        6: prior['B7'].value,
        7: prior['B8'].value,
        9: prior['B10'].value,
        10: prior['B11'].value,
        11: prior['B12'].value,
    }
    assert all(value is not None for value in contact.values())
    with ZipFile(args.template) as src:
        sheet = src.read('xl/worksheets/sheet1.xml').decode('utf-8')
        shared = src.read('xl/sharedStrings.xml').decode('utf-8')
        root = ET.fromstring(shared)
        old_count = int(root.attrib['count'])
        old_unique = int(root.attrib['uniqueCount'])
        existing = {}
        for index, si in enumerate(root.findall(f'{MAIN}si')):
            text = ''.join(t.text or '' for t in si.iter(f'{MAIN}t'))
            existing.setdefault(text, index)
        assert len(root.findall(f'{MAIN}si')) == old_unique
        additions = []
        string_uses = 0

        def string_index(value):
            nonlocal string_uses
            string_uses += 1
            value = str(value)
            if value not in existing:
                existing[value] = old_unique + len(additions)
                additions.append(f'<si><t>{escape(value)}</t></si>')
            return existing[value]

        changes = {row: [('B', value, 'string')] for row, value in contact.items()}
        for offset, order in enumerate(orders):
            row = 17 + offset
            sequence = order['sequence_5to3']
            if (len(sequence) != order['length'] or not 1 <= len(sequence) <= args.max_length
                    or re.fullmatch('[ACGT]+', sequence) is None):
                raise ValueError(f'Invalid sequence or length for order row {row}')
            changes[row] = [
                ('B', order['order_name'], 'string'),
                ('C', sequence, 'string'),
                ('D', len(sequence), 'formula_cache'),
                ('E', 2, 'number'),
                ('F', 1, 'number'),
                ('G', 'tPAGE', 'string'),
                ('J', 'dd H₂O', 'string'),
                ('K', 100, 'number'),
            ]

        for row, edits in changes.items():
            pattern = re.compile(rf'<row r="{row}"[^>]*>.*?</row>', re.DOTALL)
            found = pattern.search(sheet)
            if found is None:
                raise ValueError(f'Missing original template row {row}')
            updated = found.group(0)
            for col, value, kind in edits:
                updated = replace_cell(updated, f'{col}{row}', value, kind, string_index)
            sheet = sheet[:found.start()] + updated + sheet[found.end():]

        shared = shared.replace('</sst>', ''.join(additions) + '</sst>', 1)
        shared = re.sub(r'(<sst[^>]*\bcount=")(\d+)(")',
                        lambda m: f'{m.group(1)}{old_count + string_uses}{m.group(3)}', shared, count=1)
        shared = re.sub(r'(<sst[^>]*\buniqueCount=")(\d+)(")',
                        lambda m: f'{m.group(1)}{old_unique + len(additions)}{m.group(3)}', shared, count=1)
        ET.fromstring(sheet)
        ET.fromstring(shared)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with ZipFile(args.output, 'w') as dst:
            for item in src.infolist():
                data = src.read(item.filename)
                if item.filename == 'xl/worksheets/sheet1.xml':
                    data = sheet.encode('utf-8')
                elif item.filename == 'xl/sharedStrings.xml':
                    data = shared.encode('utf-8')
                dst.writestr(item, data)
    print(json.dumps({'output': str(args.output), 'order_rows': len(orders),
                      'unique_oligos': len({x['sequence_5to3'] for x in orders}),
                      'changed_zip_entries': ['xl/worksheets/sheet1.xml', 'xl/sharedStrings.xml']}))


if __name__ == '__main__':
    main()
