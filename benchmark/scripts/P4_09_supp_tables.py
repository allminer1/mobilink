"""Phase 5: supplementary tables of manuscript draft v2 that come from Phase 4 (nothing is recomputed).

Table S10  the survey of 241 public runs: sample, kits, exclusions, per-run and class results, contrasts, jackknife
           replicates, project estimates, within-project slopes, sensitivity analyses, verdicts
Table S11  sub-classes (written by scripts/P4_06c_subclass_summary.py)
Table S7c, S7d  run time and release check of package version 1.0
Table S12  reanalysis of a published comparison, and one sample with several kits (parts that exist are copied)
Output: manuscript/supplement/Table_S1{0,2}*.tsv
usage: envs/linkstat_dev/bin/python -I scripts/P4_09_supp_tables.py
"""
import os
import shutil

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
O = f"{L}/manuscript/supplement"
rd = lambda p: pd.read_csv(f"{L}/{p}", sep="\t")  # noqa: E731

# S10a: the runs, analysed and rejected, with the kit evidence
fin = rd("tables/P4/confirm/p4_final.tsv").assign(status="analysed")
kit = rd("public/screen/p4_kit_probe.tsv")[["run_accession", "kit_from_reads", "kit_stated", "chemistry_stated"]]
rej = rd("public/screen/p4_manual_exclusions.tsv")
rej = rej.assign(status="rejected: " + rej.reason)[["run_accession", "status"]]
s = pd.concat([fin.merge(kit, on="run_accession", how="left"), rej[~rej.run_accession.isin(fin.run_accession)]], ignore_index=True)
s.to_csv(f"{O}/Table_S10a_survey_runs.tsv", sep="\t", index=False)
print("S10a:", int((s.status == "analysed").sum()), "analysed,", int((s.status != "analysed").sum()), "rejected")

copies = [
    ("tables/P4/confirm/per_run.tsv", "Table_S10b_per_run_test.tsv"),
    ("tables/P4/confirm/partB.tsv", "Table_S10c_per_run_cut_and_naive.tsv"),
    ("tables/P4/confirm/classes.tsv", "Table_S10d_class_estimates.tsv"),
    ("tables/P4/confirm/jackknife.tsv", "Table_S10e_jackknife_replicates.tsv"),
    ("tables/P4/confirm/projects.tsv", "Table_S10f_project_estimates.tsv"),
    ("tables/P4/confirm/contrasts.tsv", "Table_S10g_contrasts.tsv"),
    ("tables/P4/confirm/h1c.tsv", "Table_S10h_within_project_slopes.tsv"),
    ("tables/P4/confirm/pooled.tsv", "Table_S10i_pooled_tests.tsv"),
    ("tables/P4/confirm/verdict.tsv", "Table_S10j_verdicts.tsv"),
    ("tables/P4/sens/sens_samples.classes.tsv", "Table_S10k_sensitivity_samples_classes.tsv"),
    ("tables/P4/sens/sens_samples.h3a.tsv", "Table_S10l_sensitivity_samples_contrasts.tsv"),
    ("tables/P4/sens/sens_variants.classes.tsv", "Table_S10m_sensitivity_variants_classes.tsv"),
    ("tables/P4/sens/sens_variants.h3a.tsv", "Table_S10n_sensitivity_variants_contrasts.tsv"),
    ("tables/P4/sens/control_screen.tsv", "Table_S10o_control_read_screen.tsv"),
    ("tables/P4/sens/chimera_counts.tsv", "Table_S10p_chimera_filter_counts.tsv"),
    ("tables/P4/v1check/P4_09_runtime.tsv", "Table_S7c_runtime_version_1.tsv"),
    ("tables/P4/v1check/P4_09_release_check.tsv", "Table_S7d_release_check_version_1.tsv"),
    ("tables/P4/sub/subclass.mge_difference.tsv", "Table_S11b_mge_kind_difference.tsv"),
    ("tables/P4/sub/subclass.composition.tsv", "Table_S11c_arg_composition.tsv"),
    ("public/screen/p4_07_runs.tsv", "Table_S12a_published_comparison_runs.tsv"),
    ("tables/P4/reanalysis/their_runs.tsv", "Table_S12b_published_calls_per_run.tsv"),
    ("tables/P4/reanalysis/their_env.tsv", "Table_S12c_published_calls_per_environment.tsv"),
    ("tables/P4/reanalysis/their_length_classes.tsv", "Table_S12d_published_calls_by_read_length.tsv"),
    ("tables/P4/reanalysis/their_slope.tsv", "Table_S12e_published_calls_length_models.tsv"),
    ("tables/P4/reanalysis/ours.groups.tsv", "Table_S12f_same_runs_this_pipeline_environments.tsv"),
    ("tables/P4/reanalysis/ours.runs.tsv", "Table_S12g_same_runs_this_pipeline_runs.tsv"),
    ("public/screen/p4_07b_runs.tsv", "Table_S12h_kit_sets_runs.tsv"),
    ("tables/P4/kits/kits.groups.tsv", "Table_S12i_kit_sets_cells.tsv"),
    ("tables/P4/kits/kits.slopes.tsv", "Table_S12j_kit_sets_slopes.tsv"),
    ("tables/P4/kits/kits.profiles.tsv", "Table_S12k_kit_sets_profiles.tsv"),
]
for src, dst in copies:
    if os.path.exists(f"{L}/{src}"):
        shutil.copy(f"{L}/{src}", f"{O}/{dst}")
    else:
        print("not yet there:", src)
print(len([x for x in os.listdir(O) if x.startswith(("Table_S10", "Table_S11", "Table_S12"))]), "files of S10-S12")
