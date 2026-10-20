# Data dictionary (frozen CSVs)

All CSVs live in `data/frozen/`.

## cv_results.csv
Per-target 5-fold scaffold-split ROC-AUC and bootstrap 95% CI.
Columns: `target`, `fold`, `roc_auc`, `pr_auc`, `mcc`, `balanced_accuracy`, `n_actives`, `n_inactives`.

## coupling_stats.csv
Biosynthesis × target coupling test results.
Columns: `npclassifier_class`, `target`, `n_compounds`, `n_scaffolds`, `n_activity_records`,
`mw_u_stat`, `mw_p`, `rank_biserial`, `rb_lo`, `rb_hi`, `bh_q`, `scaffold_perm_p`.

## coupling_class_mean.csv / coupling_class_log2enrich.csv
Class × target mean predicted probability and log2 enrichment vs panel baseline.

## coupling_pathway_mean.csv / coupling_pathway_log2enrich.csv
Pathway-level aggregates (NPClassifier pathway axis).

## perm_model.csv
Label-permutation null for the multi-task MLP (50 within-fold shuffles).
Columns: `permutation`, `null_macro_roc_auc`.

## perm_coupling.csv
Scaffold-stratified class-label permutation null for the aggregate coupling statistic
Σ −log10(q) across 1,000 permutations.
Columns: `permutation`, `null_aggregate_stat`.

## admet_leads.csv
Pareto-filtered 50-compound set with Lipinski, Veber, QED, SAscore, PAINS and BOILED-Egg
descriptors.
Columns: `compound_id`, `smiles`, `MW`, `logP`, `HBA`, `HBD`, `TPSA`, `rotatable`,
`QED`, `SAscore`, `PAINS_pass`, `BBB`, `GI_absorption`, `in_AD_<target>` (one per target).

## final_leads.csv
Structural-arm leads NP01–NP08 with biosynthetic class and assigned target.

## leads_diverse.csv
Butina-diversified 117-compound set.

## dock_scores.csv / dock_summary.csv / smina_vs_gnina.csv
Per-pose smina kcal/mol score, GNINA CNN score and CNN affinity, grouped per
target–ligand pair.

## gnina_scores.csv
GNINA CNN_score and CNN_affinity for every (target, ligand) pair.

## decoys.csv
Property-matched decoys per target (MW ± 10%, logP ± 1, HBA/HBD band, Tanimoto < 0.4).

## md_selection.csv / md_shortlist.csv / md_inputs.csv
Target–lead complexes advanced to MD, system box sizes, residue counts.

## unigbsa_scores.csv
Single-frame dock-pose Uni-GBSA ΔG per complex.

## evidence_matrix.csv
Cross-arm consistency score per complex (ML prob, docking rank, MD stability, ΔG rank).

## targets.tsv
29-target panel: gene, UniProt, PDB code, chain, binding-site definition.

## biotransform.csv / biotransform_control.csv
SyGMa biotransformation gain for parent + metabolite polypharmacology vs size-matched
random null.

## calibration.csv
Brier score and expected calibration error on held-out OOF predictions.

## compound_key.csv
Mapping between NP01–NP08 shorthand and InChIKey + SMILES.
