# Azenta workbook adapter

This adapter supports the specific web-order XLSX layout previously tested by the project. It preserves all ZIP/OOXML entries except the first worksheet and shared-string table. Generic workbook rewrites can remove vendor metadata and cause import failure.

Download the current template from the vendor yourself and inspect its layout. The repository does not distribute the vendor's workbook, filled real orders or customer details. A newer template is not automatically compatible with this adapter.

The order JSON is a list of `{ "order_name": "example-F", "sequence_5to3": "ACGT…", "length": 20 }` records. Sequences must contain A/C/G/T only, and lengths must match. Use full nucleotide sequences in actual JSON, not the abbreviated example.

```sh
uv run scripts/fill_azenta_template_lossless.py \
  --template orders/current_vendor_template.xlsx \
  --contact-source orders/local_contact_workbook.xlsx \
  --orders orders/primers.json \
  --output orders/new_order.xlsx \
  --expected-count 16 --max-length 59
```

All paths are mandatory; existing outputs are refused. The contact source must have the expected `订购单` sheet and contact fields at B5:B12. The original template and contact source must not be overwritten.

Current fixed order settings are **2 OD, one tube, tPAGE, dd H₂O, 100 µM**. The adapter does not choose purification or synthesis scale based on scientific optimality and is not an arbitrary supplier/modified-oligo exporter. Inspect these settings before using it for a different task.

Read back names, sequences, lengths, total rows, purification, modifications and delivery settings before submission. For a liquid-filled worksheet, choose a delivery form containing liquid on the vendor web page. Successful spreadsheet import does not itself establish final submission, acceptance or delivery date. Supplier schedules and prices must be confirmed for the actual order.
