# RAVEN → raven-toolbox migration reference

A function-by-function map from the MATLAB RAVEN Toolbox to raven-toolbox (and cobrapy where
appropriate). For each RAVEN function we record one of four outcomes:

* ✅ **ported** — there's a direct raven-toolbox replacement; see the link.
* 🗒️ **cheatsheet** — cobrapy already covers it; a one-liner or short idiom does the job
  (recorded under each row).
* ⛔ **not ported** — explicit decision not to bring it across, with rationale.
* 🆕 **new in raven-toolbox** — functionality raven-toolbox adds that has no RAVEN counterpart.

For raven-toolbox's deliberate improvements over RAVEN (and which of them are candidates to
upstream into MATLAB RAVEN), see [IMPROVEMENTS.md](improvements.md).

## Design principle

The in-memory object is always a [`cobra.Model`](https://cobrapy.readthedocs.io). There is
no parallel RAVEN struct, no `ravenCobraWrapper`-style adapter. RAVEN fields that cobra
doesn't model natively (`rxnMiriams`, `metDeltaG`, `rxnConfidenceScores`, …) live in
cobra's `annotation` / `notes` dictionaries. raven-toolbox's YAML I/O follows the cobra YAML
standard plus the geckopy enzyme-constrained extension, so ecModels round-trip.

---

## Foundation: utilities, manipulation, I/O

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `addRxns` (equations) | ✅ [`manipulation.add_reactions_from_equations`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/add.py) | The keystone: equation-string → reactions, matching mets by id, name, or `name[comp]`; strict / auto policies for new mets and genes. |
| `addTransport` | ✅ [`manipulation.add_transport_reactions`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/transport.py) | Transport across compartments, matching mets by name, sequential `tr_NNNN` ids. cobra has no transport primitive. |
| `addRxnsGenesMets` | ✅ [`manipulation.add_reactions_from_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/transfer.py) | Copy reactions from a source model, matching mets by `name[comp]` (vs cobra's strict-by-id merge). |
| `replaceMets` | ✅ [`manipulation.replace_metabolite`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/replace.py) | Swap one metabolite for another across every reaction that uses it, merging stoichiometry on collision. |
| `removeBadRxns` | ✅ [`manipulation.remove_bad_reactions`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/bad_reactions.py) | Drops reactions that are unbalanced, duplicated, or otherwise structurally broken. |
| `mergeModels` | ✅ [`manipulation.merge_models`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/merge.py) | N-model merge, unify mets by `name[comp]` (or id), keep all reactions (id collisions renamed). cobra's `merge` is pairwise / strict-by-id. |
| `changeRxns`, `changeGrRules` | ✅ [`manipulation.change_reaction_equations`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/change.py), `change_gene_reaction_rules` | Stoichiometry change in place; batch GPR set/append (`(old) or (new)`). |
| `setParam('var', …)` | ✅ [`manipulation.set_variance_bounds`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/parameters.py) | ±% band around measured values. **Other modes (lb/ub/eq/obj/unc)**: 🗒️ cobra one-liners (`reaction.bounds`, `model.objective`, `Configuration().bounds`). |
| `setExchangeBounds` | ⛔ not ported | cobra's `model.medium = {ex_id: uptake}` covers it. |
| `simplifyModel` (gap modes) | ✅ [`manipulation.remove_dead_end_reactions`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/simplify.py), `remove_duplicate_reactions`, `constrain_reversible_reactions`, `group_linear_reactions` | Cobra-covered modes (no-flux→`find_blocked_reactions`, zero-interval, unconstrained) are 🗒️ cheatsheeted. `group_linear` is lossy (drops genes), per RAVEN. |
| `removeMets`, `removeGenes` | ✅ [`manipulation.remove_metabolites`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/remove.py), `remove_genes` | Delegate to cobra; add `by_name` cross-compartment deletion (mets) and a `blocked_reactions` policy (genes). `removeReactions` itself ⛔ — cobra's `remove_reactions` covers it. |
| `convertToIrrev` | ✅ [`manipulation.convert_to_irreversible`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/irreversible.py) | Splits reversible non-exchange reactions into forward + `_REV`. Adopted from geckopy. |
| `expandModel` | ✅ [`manipulation.expand_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/expand.py) | Splits OR-GPR (isozyme) reactions into one per AND-clause. Adopted from geckopy. |
| `mergeCompartments` | ✅ [`manipulation.merge_compartments`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/compartments.py) | Collapse a multi-compartment model into one; deduplicate identical reactions; optionally drop one-met collapses (`drop_single_metabolite_reactions`). |
| `copyToComps` | ✅ [`manipulation.copy_to_compartment`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/compartments.py) | Duplicate reactions into a target compartment (idempotent; `delete_original=True` makes it a move). |
| `mapCompartments` | ⛔ not ported | Merged predictor compartment columns before `predictLocalization`. Removed from RAVEN; the `compartment_map` argument on the localization scorers covers it. |
| `getElementalBalance` | ✅ [`utils.get_elemental_balance`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/balance.py) | Graded `balanced` / `unbalanced` / `unknown` — `unknown` catches a missing formula that cobra's `check_mass_balance` silently miscounts. |
| `checkModelStruct` (curation subset) | ✅ [`utils.check_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/validate.py) | Structured curation report. RAVEN's struct/type checks are moot in cobra. |
| `is_dnf` / GPR check (from `standardizeGrRules`) | ✅ [`utils.is_dnf`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/gpr.py), `find_non_dnf_grrules` | Lint-only half; cobra auto-normalises GPRs on assignment, so the rewriting half isn't ported. |
| `getIndexes` (`metcomps` sliver) | ✅ [`utils.parse_name_comp`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/parse.py) | The only `getIndexes` bit cobra doesn't already cover. |
| `editMiriam`, `extractMiriam` | ⛔ not ported | cobra's `.annotation` is already a `{namespace: id(s)}` dict — read/write it directly. |
| `getRxnsInComp`, `getMetsInComp` | ⛔ not ported | One-liners over cobra's `reaction.compartments` / `metabolite.compartment`. |
| `constructEquations` | ⛔ not ported | `reaction.build_reaction_string(use_metabolite_names=...)` already does both id and name equations. |
| `sortIdentifiers` | ✅ [`utils.sort_identifiers`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/sort.py) | Model-wide alphabetical sort; also via `sort_ids=` on `write_yaml_model`. |
| `generateNewIds` | ✅ [`utils.generate_new_ids`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/ids.py) | Mints fresh sequential ids after a model's existing numbering for a prefix. |
| `guessComposition` | ✅ [`utils.guess_composition`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/composition.py) | Infers a metabolite's formula from elemental mass balance where every other reaction participant already has one. |

### I/O

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `readYAMLmodel`, `writeYAMLmodel` | ✅ [`io.read_yaml_model`, `write_yaml_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/yaml.py) | Aligned to cobra's `!!omap` writer (RAVEN `fa281a1`). Adds the RAVEN-only top-level per-entry keys (inchis/deltaG/metFrom/notes, confidence_score/references/rxnFrom/deltaG, protein) into `.notes`, plus `version`/`metaData`/GECKO `ec-*`. cobra-readable output verified. |
| `exportToExcelFormat` (export only) | ✅ [`io.export_to_excel`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/excel.py) | RAVEN 5-sheet xlsx (RXNS / METS / COMPS / GENES / MODEL). Excel **import** is intentionally excluded. |
| `exportForGit` | ✅ [`io.export_for_git`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/git.py) | Standard-GEM repo layout (`model/<fmt>/…`). |
| `importYAML/SBML/Mat/Excel` | 🗒️ cobra's standard readers | `cobra.io.read_sbml_model` / `load_json_model` / etc.; Excel import not ported. |

## Reconstruction

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `getModelFromHomology` + `getBlast`/`getDiamond` | ✅ [`reconstruction.homology.get_model_from_homology`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/homology/homology.py), `run_blast`, `run_diamond`, `blast_from_table` | Core homology reconstruction with structured improvements (bidirectional / best-hits-only, AST GPR rewrite, complex policy, bitscore best-hits, DataFrame ortholog map). |
| KEGG download → species model (5 steps) | ✅ [`reconstruction.kegg.*`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/kegg/) | All five steps: `download.fetch_keggdb`, `parse.read_kegg_table` + reference model, `hmm.build_libraries`, `organism.build_kegg_model_for_organism` (no-FASTA), `query.assign_kos` + `run_hmmscan`. |
| `getPhylDist` | ✅ [`reconstruction.kegg.phyl_dist`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/kegg/taxonomy.py) (+ `parse_taxonomy`) | Lineage parsing **and** the distance matrix (RAVEN's `keggPhylDist`), regenerated from the published `taxonomy` artefact so GECKO's organism-distance kcat selection needs no `.mat`. The per-organism HMM-subsampling use is moot here (fixed prok90/euk90 libraries). |
| `getMetaCycModelForOrganism` | ⛔ not ported (and **flagged for removal from MATLAB RAVEN**) | BLAST-to-single-representatives is low-precision at every cutoff. See [IMPROVEMENTS.md](improvements.md) under `R-MetaCyc`. |

## Tasks, gap-filling, INIT, ftINIT

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `parseTaskList`, `checkTasks` | ✅ [`tasks.parse_task_list`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/tasks/tasklist.py), `check_tasks` ([tasks/check.py](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/tasks/check.py)) | `check_tasks` reuses one model across the task list (no per-task model copy) — at genome scale ~12× faster than the copy-per-task implementation it replaced. |
| `fillGaps` (connectivity mode) | ✅ [`gapfilling.connect_blocked_reactions`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/gapfilling/fill.py) | MILP (min penalty-weighted template reactions s.t. blocked reactions carry flux). Targeted mode → 🗒️ `cobra.gapfill`. |
| `findLeakMetabolite` | ✅ [`gapfilling.find_leak_metabolite`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/gapfilling/leak.py) | Finds a metabolite a model can produce or consume for free (mass-balance leak), independent of gap-filling. |
| `gapReport` | ✅ [`gapfilling.gap_report`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/gapfilling/report.py) | Summarises which blocked reactions/metabolites a gap-fill run added, removed, or left unresolved. |
| `scoreModel` | ✅ [`init.score_reactions_from_genes`, `gene_scores_from_expression`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/score.py) | RNA-seq scoring is `5·ln(level/ref)`-clamped. |
| `runINIT`, `getINITModel` (tINIT) | ⛔ not ported | tINIT is MATLAB RAVEN's legacy pre-ftINIT algorithm, kept there for backwards compatibility. raven-toolbox is a new implementation with no installed base to carry forward, so it implements only ftINIT. |
| `ftINIT`, `prepINITModel`, `ftINITInternalAlg`, `getINITSteps` | ✅ [`init.ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py), [`init.prep_init_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/prep.py), [`init.run_ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py), [`init.get_init_steps`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/steps.py) | Staged MILP + linear merge + scaling (`rescaleModelForINIT`). Validated against RAVEN on Human-GEM (Jaccard 0.975–0.980; see [the Human-GEM validation study](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/humangem-validation.md) on raven-docs). Metabolomics-based scoring (`metabolomics`/`prod_weight`) is now included, ported from RAVEN's `develop3` branch. |
| `ftINITFillGaps`, `fitTasks` (`gapFillMode` `'preMerged'`) | ✅ [`init.fill_tasks`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/taskfill.py) | Task-aware gap-filling within ftINIT; in-place `_feasible` check + bounded fill MILP (`mip_gap`, `time_limit`). |
| `mergeLinear` | ✅ [`init.merge_linear`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/merge.py) | Linear merge of unit-stoichiometry chains; bookkeeping (`group_ids`, `reversed_rxns`) to map back to the reference model. |
| `removeLowScoreGenes` | ✅ [`init.remove_low_score_genes`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/genes.py) | Final gene-prune step of ftINIT. |
| `fitTasks` (`gapFillMode` `'merge'`, tINIT's mode, backed by `fillGaps`) | ⛔ not ported | tINIT-only; see the tINIT row above. |

## Curation

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `renameModelGenes` | ✅ [`curation.rename_model_genes`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/curation/rename_genes.py) | Renames genes across the model and every GPR that references them, from a mapping table. Delegates to cobra's `rename_genes` rather than RAVEN's own regex rewrite. |
| `getGeneData`, `downloadGenomeData` | ✅ [`curation.get_gene_data`, `download_genome_data`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/curation/genome_data.py) | Resolves a local GFF3 or an NCBI accession into a gene table; fetches the GFF3 + protein FASTA pair via the NCBI Datasets v2 API. |
| `processProteinFastaFile` | ✅ [`curation.process_protein_fasta_file`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/curation/process_protein_fasta.py) | Renames a protein FASTA's headers via a gene mapping table (e.g. `get_gene_data`'s own output). |

## Omics, analysis, comparison

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `parseHPA`, `parseHPArna`, `scoreModel` | ✅ [`omics.parse_hpa`, `parse_hpa_rna`, `hpa_gene_scores`, `rna_gene_scores`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/omics/hpa.py) | Pandas-tidy DataFrames; scoring adapters reuse `score_reactions_from_genes` (single source of truth for the GPR walk). |
| `reporterMetabolites` | ✅ [`analysis.reporter_metabolites`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/reporter.py) | Exact closed-form background replaces RAVEN's Monte-Carlo (RM1 in IMPROVEMENTS). |
| `FSEOF` | ✅ [`analysis.fseof`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/fseof.py) | Regression slope + correlation, amplify/knockdown/knockout classes, gene aggregation (FS1–FS4 in IMPROVEMENTS). |
| `randomSampling` | ✅ [`analysis.sample`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/sampling.py) | Wraps cobra's flux sampling. |
| `compareFluxes` | ✅ [`analysis.compare_fluxes`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/compare_fluxes.py) | Same selection and classification; returns a tidy DataFrame instead of a printed table, so RAVEN's `nMax` / `verbose` have no counterpart. Replaces the earlier `follow_changed`, ported from `followChanged` before RAVEN removed it in favour of this. |
| `traceFluxPath` | ✅ [`analysis.trace_flux_path`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/trace_flux_path.py) | Best-first search for how flux reaches a target reaction, back through its producing reactions. |
| `walkFluxes` | ✅ [`analysis.walk_fluxes`, `FluxWalker`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/walk.py) | Interactive step-by-step flux explorer around a starting reaction. |
| `getMinNrFluxes` | ✅ [`analysis.get_min_nr_fluxes`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/min_flux_count.py) | Minimises the count of active reactions (MILP) rather than their summed flux (pFBA). |
| `getFluxZ`, `analyzeSampling` | ✅ [`analysis.get_flux_z`, `analyze_sampling`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/differential_sampling.py) | Per-reaction Z-score between two sampled flux sets, combined with an expression t-score into concordant / expression-only / flux-only classes. Takes *samples × reactions*, matching `random_sampling`'s own orientation (RAVEN takes the transpose). |
| `fitParameters` | ✅ [`biomass.fit_parameters`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/biomass/fit_parameters.py) | Fits stoichiometric-matrix coefficients (e.g. a maintenance ATP cost) by nonlinear least squares against measured fluxes. |
| `compareMultipleModels` | ✅ [`comparison.compare_models`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/comparison/compare.py) | Tidy DataFrames (reactions / mets / genes / subsystems presence + pairwise Jaccard + optional `check_tasks` pass/fail). Plotting and tSNE/MDS are 🗒️ one-liners in seaborn/scikit-learn; intentionally not in the function. |
| `runDynamicFBA` | ⛔ not ported | Established Python implementations exist: [`dfba`](https://pypi.org/project/dfba/) (Pinheiro et al.; CVODES-backed), [`reframed`](https://pypi.org/project/reframed/) (Machado lab), [`mewpy`](https://pypi.org/project/mewpy/) (Cunha lab). Cobrapy itself has none, but re-porting would duplicate maintained prior art. |

## Localisation (Phase 7)

| RAVEN | raven-toolbox | Notes |
|---|---|---|
| `predictLocalization` | ✅ [`localization.predict_localization`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py) | Deterministic MILP (not simulated annealing). Caller-passed `reactions_to_relocate` set (everything else pinned). Multi-compartment by default: primary "free", extras pay `multi_compartment_penalty`. Tolerates incomplete models (no silent reaction removal). `apply=False` returns a `LocalizationProposal` diff. Real-data validation against curated yeast-GEM in [yeast_localization_benchmark.md](../studies/yeast_localization_benchmark.md), and against MATLAB's own `predictLocalization` head-to-head in [predictlocalization_comparison.md](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/predictlocalization-comparison.md) (raven-docs). |
| `getWoLFScores`, `parseScores('wolf')` | ➖ Not ported (dropped) | WoLF PSORT is superseded by modern multi-label predictors and the COMPARTMENTS evidence database; use `load_deeploc` / `load_mulocdeep` / `load_compartments` instead. |
| `parseScores('deeploc')` | ✅ [`localization.load_deeploc`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/scores.py) | DeepLoc 2 per-protein CSV (Protein_ID / Localizations / Signals + one column per compartment). `compartment_map` (e.g. `DEFAULT_COMPARTMENT_MAP`) renames labels to model ids. |
| — (new) | ✅ [`localization.load_mulocdeep`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/scores.py) | MULocDeep wide table (id + per-compartment probabilities); auto-detects delimiter. |
| — (new) | ✅ [`localization.load_compartments`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/scores.py) | COMPARTMENTS (jensenlab.org) channel TSV; aggregates per gene×compartment, `min_confidence` filters weak annotations. |
| — (new) | ✅ [`localization.load_uniprot`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/scores.py) | Curated UniProtKB `Subcellular location [CC]` TSV export; qualitative (1.0 per annotated compartment). For yeast-GEM use the ordered-locus id column. |

## Things deliberately not ported

* **`ravenCobraWrapper` / RAVEN struct adapter** — cobra is the canonical object; no parallel struct.
* **`checkModelStruct` struct/type checks** — moot in cobra.
* **`runDynamicFBA`** — see Omics/analysis row.
* **`getMetaCycModelForOrganism`** — see Reconstruction row; flagged for upstream removal.
* **`getPhylDist` per-organism HMM subsampling** — fixed prok90/euk90 libraries make it moot (the distance matrix itself **is** ported, as `reconstruction.kegg.phyl_dist`, for GECKO).
* **`editMiriam`, `extractMiriam`, `getRxnsInComp`, `getMetsInComp`, `constructEquations`, `getIndexes`** (most), **`setExchangeBounds`**, **most `setParam` modes**, **`getBlastFromExcel` Excel branch** — cobra one-liners; recorded above.

## "New in raven-toolbox" entry points

These are not direct RAVEN ports:

* [`comparison.compare_models`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/comparison/compare.py) — already in RAVEN, but the raven-toolbox version returns tidy DataFrames suitable for downstream analysis (RAVEN's version embeds heatmap/tSNE plotting).
* [`reconstruction.homology.run_blast` / `run_diamond` / `blast_from_table`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/homology/blast.py) — generic subprocess wrappers; a DataFrame is the canonical "hits" object (no `blastStructure` struct).
* [`binaries.resolve_binary`, `ensure_binary`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/binaries.py) — version-pinned release-ZIP registry for external tools (BLAST/DIAMOND/HMMER), SHA256-verified.
* [`localization.LocalizationProposal`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py) — the diff-preview mode (`apply=False`).
