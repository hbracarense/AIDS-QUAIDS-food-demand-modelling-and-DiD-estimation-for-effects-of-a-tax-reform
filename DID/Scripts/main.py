from pathlib import Path
import numpy as np
import pandas as pd
from utils import read_elas_matrix_with_ci 
from functions import simulate_post_shares, build_panel, run_did_set, parametric_bootstrap, placebos_randomization, run_quaids_checks

def main(base_path, elas_path, pre_path, s1_path, s2_path, outdir, boot_B=0, placebo_runs=0, seed=42):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Base e objetos
    base = pd.read_csv(base_path)
    if "id" not in base.columns: base["id"] = np.arange(len(base))

    groups = list(range(1,7))
    Wpre = pd.DataFrame({f"w{g}": base[f"w_despesahat{g}"].astype(float).values for g in groups})
    Ppre = pd.DataFrame({f"p{g}": base[f"preco_por_kg{g}"].astype(float).values for g in groups})
    P1   = pd.DataFrame({f"p{g}": base[f"preco_com_reforma1{g}"].astype(float).values for g in groups})
    P2   = pd.DataFrame({f"p{g}": base[f"preco_com_reforma2{g}"].astype(float).values for g in groups})
    Xtot = base["gasto_total_atualhat"].astype(float).values

    # Elasticidades Marshall status quo (com SD aprox p/ bootstrap)
    m_pre, sd_pre = read_elas_matrix_with_ci(pd.read_excel(elas_path, "Marshall - status quo"))

    # DiD baseline (pontual)
    W1 = simulate_post_shares(Wpre, Ppre, P1, Xtot, m_pre)
    W2 = simulate_post_shares(Wpre, Ppre, P2, Xtot, m_pre)
    pan1 = build_panel(base, Wpre, W1, Ppre, P1, Xtot)
    pan2 = build_panel(base, Wpre, W2, Ppre, P2, Xtot)
    A1, BG1, BGadj1, C1 = run_did_set(pan1, "s1")
    A2, BG2, BGadj2, C2 = run_did_set(pan2, "s2")

    did_main = pd.concat([A1.assign(scenario=1), A2.assign(scenario=2)], ignore_index=True)
    did_main.to_csv(outdir/"did_main_overall.csv", index=False)
    BGadj1.to_csv(outdir/"did_bygroup_adj_s1.csv", index=False)
    BGadj2.to_csv(outdir/"did_bygroup_adj_s2.csv", index=False)
    C1.to_csv(outdir/"did_continuous_s1.csv", index=False)
    C2.to_csv(outdir/"did_continuous_s2.csv", index=False)

    # Bootstrap paramétrico
    if boot_B and boot_B > 0:
        bootA, bootBG, bootC = parametric_bootstrap(base, Wpre, Ppre, P1, P2, Xtot, m_pre, sd_pre, boot_B, seed)
        bootA.to_csv(outdir/"boot_overall_post.csv", index=False)
        bootBG.to_csv(outdir/"boot_bygroup_postg.csv", index=False)
        bootC.to_csv(outdir/"boot_continuous_posttax.csv", index=False)

    # Placebos (randomization inference)
    if placebo_runs and placebo_runs > 0:
        plA, plBG, plC = placebos_randomization(base, Wpre, Ppre, P1, P2, Xtot, m_pre, placebo_runs, seed)
        plA.to_csv(outdir/"placebo_overall_post.csv", index=False)
        plBG.to_csv(outdir/"placebo_bygroup_postg.csv", index=False)
        plC.to_csv(outdir/"placebo_continuous_posttax.csv", index=False)

    # QUAIDS checks (pre, s1, s2) + índices Stone vs Translog
    overall_list=[]; blocks_list=[]; idx_list=[]
    ov, bl, idxc = run_quaids_checks("pre", pre_path, base, [f"preco_por_kg{k}" for k in range(1,7)], outdir)
    overall_list.append(ov); blocks_list.append(bl); idx_list.append(idxc)
    ov, bl, idxc = run_quaids_checks("1", s1_path,  base, [f"preco_com_reforma1{k}" for k in range(1,7)], outdir)
    overall_list.append(ov); blocks_list.append(bl); idx_list.append(idxc)
    ov, bl, idxc = run_quaids_checks("2", s2_path,  base, [f"preco_com_reforma2{k}" for k in range(1,7)], outdir)
    overall_list.append(ov); blocks_list.append(bl); idx_list.append(idxc)
    pd.concat(overall_list, ignore_index=True).to_csv(outdir/"quaids_wald_overall_all.csv", index=False)
    pd.concat(blocks_list,  ignore_index=True).to_csv(outdir/"quaids_wald_blocks_all.csv", index=False)
    pd.concat(idx_list,     ignore_index=True).to_csv(outdir/"quaids_price_index_compare_all.csv", index=False)

    # Pequeno sumário no console
    print("Concluído. Principais arquivos em", outdir)
    print("- did_main_overall.csv, did_bygroup_adj_s*.csv, did_continuous_s*.csv")
    if boot_B and boot_B > 0:
        print("- boot_overall_post.csv, boot_bygroup_postg.csv, boot_continuous_posttax.csv")
    if placebo_runs and placebo_runs > 0:
        print("- placebo_overall_post.csv, placebo_bygroup_postg.csv, placebo_continuous_posttax.csv")
    print("- quaids_wald_overall_all.csv, quaids_wald_blocks_all.csv, quaids_price_index_compare_all.csv")
