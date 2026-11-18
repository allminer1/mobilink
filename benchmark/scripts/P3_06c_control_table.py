"""P3 step 3.5, exploratory: the registered per-run estimates with and without the suspected control-strand reads
(scripts/P3_06c_control_exclude.sh). Output: tables/P3/P3_06c_control_effect.tsv (one row per run and pair)."""
import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
RUNS = ["SRR23278716", "SRR23278685", "SRR23278728", "SRR23278753", "ERR14825711", "ERR14825617", "SRR33225186"]
reg = pd.concat([pd.read_csv(f"{L}/tables/P3/{d}/per_run.tsv", sep="\t") for d in ("confirm", "confirm_ext")])
rows = []
for run in RUNS:
    a = pd.read_csv(f"{L}/tmp/p35c/{run}/results.tsv", sep="\t")
    n = sum(1 for _ in open(f"{L}/tmp/p35c/{run}/control_suspects.txt"))
    for r in a.itertuples():
        b = reg[(reg.run == run) & (reg.class_a == r.class_a) & (reg.class_b == r.class_b)].iloc[0]
        rows.append(dict(run=run, env=b.env, pair=f"{r.class_a}:{r.class_b}", suspected_control_reads=n, molecules=b.molecules,
                         net_w_registered=b.net_w, net_w_without=r.net_w, net_w_lo_registered=b.net_w_lo, net_w_hi_registered=b.net_w_hi,
                         net_w_p_registered=b.net_w_p, net_w_p_without=r.net_w_p,
                         net_w_E_registered=b.net_w_E, net_w_E_without=r.net_w_E, net_w_O=r.net_w_O,
                         net_mol_registered=b.net_mol, net_mol_without=r.net_mol,
                         net_mol_p_registered=b.net_mol_p, net_mol_p_without=r.net_mol_p,
                         hosts_effective_registered=b.hosts_effective, hosts_effective_without=r.hosts_effective))
t = pd.DataFrame(rows)
t["net_w_diff"] = t.net_w_without - t.net_w_registered
t["net_mol_diff"] = t.net_mol_without - t.net_mol_registered
t.to_csv(f"{L}/tables/P3/P3_06c_control_effect.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
print(t[["run", "env", "pair", "suspected_control_reads", "net_w_registered", "net_w_without", "net_w_diff", "net_w_p_registered",
         "net_w_p_without", "net_mol_registered", "net_mol_without", "net_mol_diff", "net_mol_p_registered", "net_mol_p_without"]].round(4).to_string(index=False))
m = t[t.pair == "ARG:MGE"]
print("ARG:MGE largest |change| net_w", np.nanmax(np.abs(m.net_w_diff)), "net_mol", np.nanmax(np.abs(m.net_mol_diff)))
