# Public-core validation

Date: 2026-10-19. Platform: macOS, Apple Silicon, Python 3.12.14. Dependency installation: `uv sync --locked`, uv 0.13.1.

- **49 tests passed**, including three synthetic CLI workflows: PCA layout, bidirectional Sanger clone selection, and structure-sequence to plasmid map/workbooks.
- The synthetic Sanger run exports all three passing clones, rejects the parental allele, reports the absent fifth sample and preserves the input AB1 hashes.
- PCA tests verify complete target coverage, 5′→3′ orientation, length and junction bounds, and refusal to overwrite or silently accept unsupported oligo lengths.
- The structure design roundtrip verifies the encoded POI and workbook/manifest production using the public template.
- Public template POI replaced by a same-length synthetic placeholder. Reconstructed backbone sequence outside the POI is identical to the laboratory template.
- Vendor-adapter smoke check used a locally supplied original template and artificial contact/order records. Exactly two OOXML entries changed; the other 24 entries were byte-identical. Vendor workbooks are not distributed. Compatibility still depends on the vendor layout.
- Candidate source/document links were checked; personal absolute paths and common credential patterns were not found. Generated runs, environments and laboratory batches are excluded from Git.

The test workflow is configured for macOS and Linux; remote CI results are recorded by GitHub Actions, not inferred from the local run. Synthetic traces are idealized fixtures, not empirical validation of difficult chromatograms. Dataset-specific regressions remain in the laboratory workspace, so this public suite is not claimed to reproduce every historical experimental batch.
