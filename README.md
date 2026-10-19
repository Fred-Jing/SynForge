# Primer Build

A local Python toolkit for **gene assembly design, expression-construct preparation, and sequencing-based clone selection**. It grew out of a wet-lab workflow for converting protein designs into plasmid maps and synthesis oligos.

The priorities are robust synthesis and good sequence-level expression potential. Computational scores are screening proxies; they do not establish synthesis success, expression yield, solubility or activity.

[中文说明](docs/README.zh-CN.md) · [Sanger workflow](docs/SANGER_WORKFLOW.md) · [Design and method boundaries](docs/METHODS.md)

## What is included

|Component|Status|
|---|---|
|Structure → coding DNA → plasmid maps and oligo workbooks|CLI; specialized to the supplied pET22B-derived backbone and overlap-PCR/Gibson workflow|
|Alternating PCA layout design|Reusable heuristic engine and standalone CLI for one exact DNA target|
|Codon and expression-related screening|Codon search, GC/repeat checks and ViennaRNA-based initiation/local-folding metrics|
|Sanger analysis|Config-driven AB1 quality, trace and bidirectional evidence analysis; exports **all passing clones**|
|Whole-plasmid analysis primitives|Circular normalization, sequence comparison and BAM QC library functions; no generic batch CLI yet|
|Azenta workbook adapter|Preserves the original workbook's OOXML; explicit input paths required; supports one known vendor layout|

Point-mutation primer design is not yet a generic automated CLI. The PCA CLI does not design vector assembly, AF/AR primers, or synonymous sequence changes. Neither feature should be inferred from the presence of lower-level library functions.

## Install

Python **3.12**, using [uv](https://docs.astral.sh/uv/):

```sh
git clone https://github.com/Fred-Jing/primer-build.git
cd primer-build
uv sync --locked
```

Run from this checkout using `uv run`; this release is not a PyPI package. macOS Apple Silicon is the original development platform. A macOS/Linux CI recipe is supplied in `ci/github-actions.yml`; remote CI is not enabled yet because the current GitHub credentials cannot write Actions workflows. Windows native support is not claimed: `pysam` and other scientific dependencies may require Linux/WSL or additional setup.

## Try the synthetic examples

These commands generate artificial DNA and ideal synthetic ABIF traces. No laboratory data or vendor account is needed. Use a new output directory on each run.

```sh
uv run scripts/make_demo.py --output demo

# Bidirectional Sanger analysis
uv run scripts/analyze_sanger_batch.py --config demo/batch.json

# Alternating PCA oligos for one fixed DNA sequence
uv run scripts/design_pca.py \
  --target demo/pca_target.fasta --output demo/pca_layout.json

# Structure-sequence input to maps and overlap/Gibson oligos
uv run scripts/design_project.py \
  --project demo/design --chain A \
  --map-start DEMO-1 --primer-start DEMO-1-1-1-FW \
  --candidates 2 --workers 1
```

The mmCIF demo contains a synthetic polymer sequence, not a predicted structure. The design command uses sequence information, not 3D geometry. The expected Sanger outcome is **three passing clones**, one rejected parental allele, and one missing sample; inspect `demo/results/usable_samples.md` and `report.html`. Synthetic results demonstrate software behavior, not experimental accuracy.

## Use your own data

- **Structure design:** put PDB/mmCIF files in `<project>/structures/`; specify the chain and numbering explicitly. Outputs are `maps/*.dna`, `maps/statistics_summary.xlsx`, `all_primers.xlsx` and `run_manifest.json`.
- **PCA:** provide exactly one unambiguous DNA FASTA record. The wrapper supports oligos up to 60 nt, consistent with its Primer3 structural-screening limit. Output coordinates are zero-based, half-open; oligo sequences are always written 5′→3′. It does not alter the input DNA.
- **Sanger:** copy the generated `demo/batch.json`, replace sample mapping, references, expected substitutions and required regions. See the [configuration and quality workflow](docs/SANGER_WORKFLOW.md). References and raw reads stay unchanged.
- **Orders:** obtain the current vendor workbook yourself; templates, contacts and filled real orders are not bundled. See [workbook usage](docs/ORDER_UPLOAD.md).

Existing results are protected by default. Only use explicit overwrite options after checking the corresponding CLI; the PCA wrapper and example generator require new outputs.

## Biological defaults and limitations

The supplied backbone is a particular **pET22B-derived intracellular expression construct without pelB**, not a universal pET22B template. Defaults include starting Met, TAG for nonstandard polymer residues, and `POI-GGSGGSG-C-SLE-His6-TGA`. Additional internal stops, tags, host context and natural cysteines need task-specific review. The template contains a synthetic placeholder POI which the design pipeline replaces; it is not an experimental result.

Synonymous optimization preserves encoded targets under the declared nonstandard-residue policy. RNA folding, codon usage and primer thermodynamics cannot guarantee expression or assembly. PCA uses bounded beam search: a returned layout is not guaranteed globally optimal, and search failure does not prove infeasibility. Sanger PASS validates the configured high-quality covered region, **not the entire plasmid**.

## Development

```sh
uv run python -m unittest discover -s tests -v
```

Tests cover circular coordinates, variant handling, Sanger quality/coverage and trace gates, structural input parsing, expression-screening invariants, PCA screening, workbook output and synthetic CLI roundtrips. Real experimental batches and their dataset-dependent regressions remain outside this public core.

See [validation](docs/VALIDATION.md), [contributing](CONTRIBUTING.md) and [methods](docs/METHODS.md). Source code is MIT licensed; dependencies retain their respective licenses. SnapGene, GENEWIZ and Azenta are third-party names; this project is not affiliated with them.
