# Analysis & comparison

## Analyses — {mod}`raven_toolbox.analysis`

RAVEN analyses that are not in cobrapy's core:

- {func}`raven_toolbox.analysis.reporter_metabolites` — Reporter Metabolites, an
  around-metabolite gene-score test. raven-toolbox uses an exact closed-form background in
  place of RAVEN's Monte-Carlo.
- {func}`raven_toolbox.analysis.fseof` — Flux Scanning based on Enforced Objective Flux:
  regression slope + correlation, with amplify / knockdown / knockout classes and gene
  aggregation.
- {func}`raven_toolbox.analysis.random_sampling` — random-objective flux sampling (wraps
  cobra's samplers); {func}`raven_toolbox.analysis.find_good_reactions` is the companion
  screen.
- {func}`raven_toolbox.analysis.compare_fluxes` — what changed between two flux
  distributions: a DataFrame sorted by the size of the change, plus the reactions that
  were turned on, turned off or reversed. Narrow it to part of the network by naming
  metabolites (`metabolite_list=["ATP", "NADH"]`).

## Comparison — {mod}`raven_toolbox.comparison`

{func}`raven_toolbox.comparison.compare_models` compares any number of models and returns tidy
DataFrames (reaction / metabolite / gene / subsystem presence, pairwise Jaccard, and an
optional `check_tasks` pass/fail matrix). Plotting and tSNE/MDS are deliberately left out —
they are one-liners in seaborn / scikit-learn on the returned frames.
