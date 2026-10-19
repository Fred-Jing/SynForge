# SynForge — public core

- Use Python 3.12 and the project environment: `uv sync --locked` and `uv run …`.
- Keep changes focused. Laboratory data, contact information and real order files are not repository fixtures.
- Use explicit references, expected edits and validation regions. Never fill missing sequencing evidence with the expected reference.
- Distinguish predicted synthesis/expression metrics, local Sanger validation, whole-plasmid validation and measured protein expression.
- Preserve original input/reference files. Write generated artifacts into new output locations.
- Test with `uv run python -m unittest discover -s tests -v`. New scientific behavior requires correctness and relevant risk checks; improved scores alone do not establish improved wet-lab performance.
- PCA layout search is bounded and heuristic. Do not describe search failure as proof of infeasibility or output as globally optimal.
- The supplier adapter supports one known workbook layout. Preserve untouched OOXML components and require explicit user inputs.
