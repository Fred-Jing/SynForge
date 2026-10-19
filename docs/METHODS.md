# Methods and scope

## Design layers

1. Read a specified structure chain using sequence declarations where available. Handle modified polymer residues explicitly; do not interpret water or arbitrary ligands as amino acids.
2. Generate synonymous coding candidates while preserving declared amino-acid identities and TAG sites. Screen codon usage, local GC, homopolymers and repeated motifs.
3. Rank expression-related proxies with ViennaRNA: initiation-region accessibility, local folding and internal Shine–Dalgarno-like motifs. These are model outputs, not calibrated expression predictions.
4. Evaluate assembly constraints. The general structure CLI currently uses the specialized backbone's overlap-PCR/Gibson strategy. PCA is a separate target-DNA layout engine.
5. Retain explicit references and sequence/parameter provenance. Clone validation uses observed sequence evidence rather than the design score.

`config/design_rules.json` records defaults and units. It is a specialized configuration, not a universal protocol. Its overlap/Gibson Tm model differs from the explicit Primer3 concentration model in `primer_build.pca`; values must not be mixed without recalculation. Nonstandard TAG decoding requires the appropriate expression system outside this software.

## PCA

`primer_build.pca.design_layout` searches candidate junctions with a bounded beam. Layouts alternate forward and reverse oligos and span the exact supplied target. Screening covers intended junction Tm, individual hairpins/self-dimers, 3′ uniqueness, homopolymers and nonadjacent oligo interactions. Intended adjacent pairing is not classified as an off-target interaction.

The standalone CLI exposes target, oligo length, junction bounds and ionic/DNA concentrations. Primer3 structural screening limits it to 60 nt oligos. It rejects unsupported input instead of silently skipping structural checks. Search heuristics and thresholds can reject a sequence for which a different method or search would work; failure is not a feasibility proof.

This command produces assembly oligos only. The experiment still needs correct outer amplification primers and a separately checked vector-integration design. Whole-oligo Tm, intended PCA junction Tm and second-round primer 3′ annealing Tm are different quantities. Specified concentrations are a reference model, not an assertion about a proprietary master mix.

Method background: [Xiong et al., Nucleic Acids Research (2004)](https://academic.oup.com/nar/article/32/12/e98/2375836) and [Slinger et al., Frontiers in Bioengineering and Biotechnology (2023)](https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2023.1208784/full). These sources motivate design tradeoffs; this repository does not implement every procedure in either paper or claim a universal optimal oligo length.

## Sequencing

Sanger uses raw AB1 calls, per-base quality and four-channel trace evidence, plus primer-anchored directional alignment against an explicit circular reference. The SEQ sidecar is not an independent replicate. Expected edits are applied only to an in-memory expected reference; source files remain unchanged. Current planned edits are substitutions; observed insertions/deletions are analyzed as abnormalities.

The default policy requires expected high-quality coverage across every required base, sufficient complementary directional coverage, and both directions at each planned edit. Low-quality or unobserved positions cannot be filled from the reference. Trace ratios are heuristic flags, not quantitative estimates of mixture fractions. Exact thresholds and remaining limitations are documented in [SANGER_WORKFLOW.md](SANGER_WORKFLOW.md).

Whole-plasmid functions in `sequencing.py` provide circular normalization, differences and BAM evidence checks. Vendor-specific sample naming, reference selection and reporting from the original laboratory batch are intentionally not exposed as a generic batch CLI.

## Evidence and release boundaries

The software originated in experiments that produced correct sequenced clones. Those experiments do not establish general success rates, superiority over other methods, improved protein expression, or the performance of every later design revision. This public distribution includes the reusable core, the current specialized backbone, and artificial examples; it does not include the original experimental batches.

The template's example POI is replaced with an artificial same-length placeholder and identifying notes are removed; the sequence outside the replaced POI is unchanged. It is not a validated construct for a user's new target. Public tests exclude the original batch-dependent synonymous-search and fallback regressions; those remain in the laboratory workspace. The core algorithms are unchanged in this export; additional CLIs connect existing functions to explicit, portable inputs.

## Dependencies and licensing

Project-authored code is MIT licensed. Python dependencies are installed separately and retain their upstream licenses, including Biopython, ViennaRNA, Primer3/primer3-py, pysam and openpyxl. No dependency license is replaced by this repository's MIT license. Binary redistribution or embedding those dependencies must follow the relevant upstream terms.
