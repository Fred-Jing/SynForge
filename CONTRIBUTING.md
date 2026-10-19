# Contributing

Use Python 3.12 and `uv sync --locked`, then run:

```sh
uv run python -m unittest discover -s tests -v
```

Keep scientific changes small and reviewable. Explain the problem, biological assumptions, behavior change and validation. Synthetic examples and regression checks should cover the failure being addressed; do not tune fixtures to make a failing method appear successful.

Keep sequence identities, reference versions, expected edits and coordinate conventions explicit. Predicted scores cannot establish synthesis success or expression yield. Preserve evidence differences between complete paired coverage and complementary Sanger coverage.

Do not commit raw laboratory batches, customer contact fields, filled supplier workbooks, tokens, environments or generated run outputs. The example generator intentionally uses synthetic data. Report bugs with a minimal reproducible example, software versions and redacted configuration.
