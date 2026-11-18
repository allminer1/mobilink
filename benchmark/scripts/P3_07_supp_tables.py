"""P3 step 3.7: supplementary tables of the preprint draft, assembled from the result tables (nothing is recomputed).
Output: manuscript/supplement/Table_S<N>_*.tsv.
usage: envs/linkstat_dev/bin/python -I scripts/P3_07_supp_tables.py
"""
import os
import shutil

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
O = f"{L}/manuscript/supplement"
os.makedirs(O, exist_ok=True)
rd = lambda p: pd.read_csv(f"{L}/{p}", sep="\t")  # noqa: E731

# S1: runs of both stages, accepted and rejected
kits = pd.concat([rd("public/screen/p3_project_kits.tsv"), rd("public/screen/p3x_project_kits.tsv")]).drop_duplicates("project")
a = rd("public/screen/p3_final.tsv").assign(stage=1, status="analysed")
b = rd("public/screen/p3x_final.tsv").assign(stage=2, status="analysed")
pb = pd.concat([rd("tables/P3/confirm/partB.tsv"), rd("tables/P3/confirm_ext/partB.tsv")])[["run", "counts", "meets"]]
s1 = pd.concat([a, b])[["run_accession", "stage", "env", "project", "unit", "subsample_n50", "subsample_bases", "status"]]
s1 = s1.merge(kits[["project", "kit", "chemistry", "basecaller", "source"]], on="project", how="left")
s1 = s1.merge(pb.rename(columns={"run": "run_accession", "counts": "counts_for_H1a_H2", "meets": "meets_H1a_H2"}), on="run_accession", how="left")
rej = pd.concat([rd("public/screen/p3_manual_exclusions.tsv").assign(stage=1), rd("public/screen/p3x_manual_exclusions.tsv").assign(stage=2)])
st2 = rd("public/screen/p3x_status.tsv")
print("stage-2 status columns:", list(st2.columns))
rej = rej.assign(status="rejected: " + rej.reason)[["run_accession", "stage", "status"]]
pd.concat([s1, rej], ignore_index=True).to_csv(f"{O}/Table_S1_public_runs.tsv", sep="\t", index=False)
shutil.copy(f"{L}/public/screen/p3x_status.tsv", f"{O}/Table_S1b_stage2_walk.tsv")

# S2: per run; S3: pooled and sensitivity
pd.concat([rd("tables/P3/confirm/per_run.tsv").assign(stage=1), rd("tables/P3/confirm_ext/per_run.tsv").assign(stage=2)]).to_csv(
    f"{O}/Table_S2a_per_run_test.tsv", sep="\t", index=False)
pd.concat([rd("tables/P3/confirm/partB.tsv").assign(stage=1), rd("tables/P3/confirm_ext/partB.tsv").assign(stage=2)]).to_csv(
    f"{O}/Table_S2b_per_run_cut_test.tsv", sep="\t", index=False)
pd.concat([rd("tables/P3/confirm/pooled.tsv").assign(stage=1), rd("tables/P3/confirm_ext/pooled.tsv").assign(stage=2)]).to_csv(
    f"{O}/Table_S3a_pooled.tsv", sep="\t", index=False)
for src, dst in (("tables/P3/P3_04_sensitivity.tsv", "Table_S3b_sensitivity.tsv"),
                 ("tables/P2/P2_05_mocks.tsv", "Table_S4a_mock_communities.tsv"),
                 ("tables/P2/P2_06_annot_sens.tsv", "Table_S4b_annotation_variants.tsv"),
                 ("tables/G3/G3a_settings.tsv", "Table_S4c_simulation_G3a.tsv"),
                 ("tables/G3/G3b_settings.tsv", "Table_S4d_simulation_G3b.tsv"),
                 ("tables/P3/P3_06b_cut_bias.tsv", "Table_S5a_cut_whole_gene_retention.tsv"),
                 ("tables/P3/P3_06d_netw_cut.tsv", "Table_S5b_cut_unit_level_estimator.tsv"),
                 ("tables/P3/P3_06c_control_strand.tsv", "Table_S6a_control_strand_screen.tsv"),
                 ("tables/P3/P3_06c_control_effect.tsv", "Table_S6b_control_strand_effect.tsv"),
                 ("tables/P3/P3_06_runtime.tsv", "Table_S7a_runtime.tsv"),
                 ("tables/P3/P3_06_release_check.tsv", "Table_S7b_release_check.tsv"),
                 ("tables/P3/P3_05_pilot.nulls.tsv", "Table_S8a_contigs_published_loci_nulls.tsv"),
                 ("tables/P3/P3_05_pilot.linkstat.tsv", "Table_S8b_contigs_published_loci_mobilink.tsv"),
                 ("tables/P3/P3_05_null_sim.summary.tsv", "Table_S8c_contigs_no_linkage_simulation.tsv"),
                 ("tables/P3/P3_05_soil_contigs.linkstat.tsv", "Table_S8d_contigs_eight_soils.tsv"),
                 ("tables/P3/P3_04_soil_reads.tsv", "Table_S9a_soil_reads_deep.tsv"),
                 ("tables/P3/P3_04d_s1.linkstat.tsv", "Table_S9b_soil_reads_low_coverage.tsv")):
    if os.path.exists(f"{L}/{src}"):
        shutil.copy(f"{L}/{src}", f"{O}/{dst}")
    else:
        print("missing:", src)
print(sorted(os.listdir(O)))
