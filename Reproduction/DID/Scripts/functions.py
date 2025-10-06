from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2, norm
from statsmodels.stats.multitest import multipletests
from utils import simulate_post_shares, build_panel

# ----------------- TWFE + Cluster duplo (Cameron–Gelbach–Miller) -----------------
def _ols_within_unit(y, X, unit_series):
    ug = unit_series.astype(str)
    y_tilde = y - y.groupby(ug).transform("mean")
    X_tilde = X.copy()
    for c in X.columns:
        X_tilde[c] = X[c] - X[c].groupby(ug).transform("mean")
    X_tilde = sm.add_constant(X_tilde, has_constant="add")
    return y_tilde, X_tilde

def twfe_unit_singlecluster(panel, ycol, xcols):
    df = panel.copy()
    y = df[ycol].astype(float)
    X = df[xcols].astype(float)
    y_t, X_t = _ols_within_unit(y, X, df["unit_group_id"])
    # drop colunas sem variância após within
    keep = [c for c in X_t.columns if c == "const" or (np.nanstd(X_t[c].values) > 1e-12)]
    X_t = X_t[keep]
    res = sm.OLS(y_t, X_t, missing="drop").fit(cov_type="cluster", cov_kwds={"groups": df["id"]})
    ci = res.conf_int()
    return pd.DataFrame({
        "term": res.params.index, "coef": res.params.values, "std_err": res.bse.values,
        "p": res.pvalues.values, "ci_low": ci[0].values, "ci_high": ci[1].values
    })

def _ols_within(y, X, id_series, time_series):
    """Aplica transformação dentro (two-way FE: unit-group e time). Retorna y_tilde, X_tilde, ids, groups."""
    ug = id_series.astype(str)  # id da unidade (aqui id_dom x group já está no painel)
    tt = time_series.astype(int)

    y_tilde = y - y.groupby(ug).transform("mean") - y.groupby(tt).transform("mean") + y.mean()
    X_tilde = X.copy()
    for c in X.columns:
        X_tilde[c] = X[c] - X[c].groupby(ug).transform("mean") - X[c].groupby(tt).transform("mean") + X[c].mean()
    X_tilde = sm.add_constant(X_tilde, has_constant="add")
    return y_tilde, X_tilde

def _cgm_2way_cov(X, resid, g1, g2, ridge=1e-10):
    import numpy as np, pandas as pd
    X = np.asarray(X, dtype=float)
    e = np.asarray(resid, dtype=float).reshape(-1,1)
    p = X.shape[1]

    XT_X = X.T @ X
    lam = ridge * np.trace(XT_X) / max(p, 1)
    XT_X_inv = np.linalg.pinv(XT_X + lam * np.eye(p))

    def meat(groups):
        codes = pd.Series(groups).astype("category").cat.codes.to_numpy()
        G = codes.max() + 1
        S = np.zeros((p, p))
        for g in range(G):
            idx = (codes == g)
            if not np.any(idx):
                continue
            Xg = X[idx,:]; eg = e[idx,:]
            xe = Xg.T @ eg
            S += (xe @ xe.T)
        return S

    S1 = meat(g1)
    S2 = meat(g2)

    pairs = pd.Series(list(zip(g1, g2))).astype("category").cat.codes.to_numpy()
    G12 = pairs.max() + 1
    S12 = np.zeros((p, p))
    for pp in range(G12):
        idx = (pairs == pp)
        if not np.any(idx):
            continue
        Xp = X[idx,:]; ep = e[idx,:]
        xe = Xp.T @ ep
        S12 += (xe @ xe.T)

    V = XT_X_inv @ (S1 + S2 - S12) @ XT_X_inv
    V = 0.5*(V + V.T)  # simetrizar
    # forçar PSD
    w, Q = np.linalg.eigh(V)
    w = np.clip(w, 0.0, None)
    V_psd = (Q * w) @ Q.T
    return V_psd

def twfe_2way_cluster(panel, ycol, xcols, clus1="id", clus2="group"):
    df = panel.copy()
    y = df[ycol].astype(float)
    X = df[xcols].astype(float)

    # (use within de unidade OU unidade+tempo, conforme sua spec;
    # para a Spec A com 'post', use SOMENTE unidade)
    y_t, X_t = _ols_within_unit(y, X, df["unit_group_id"])

    # Drop colunas sem variância após o within (FAÇA AQUI apenas)
    keep_cols = [c for c in X_t.columns if c == "const" or (np.nanstd(X_t[c].values) > 1e-12)]
    X_t = X_t[keep_cols]

    res = sm.OLS(y_t, X_t, missing="drop").fit()
    V_psd = _cgm_2way_cov(X_t.values, res.resid, df[clus1].values, df[clus2].values)

    se = np.sqrt(np.clip(np.diag(V_psd), 0.0, np.inf))
    tvals = res.params.values / np.where(se>0, se, np.nan)
    pvals = 2*(1 - norm.cdf(np.abs(tvals)))
    ci_low = res.params.values - 1.96*se
    ci_high= res.params.values + 1.96*se

    return pd.DataFrame({
        "term": X_t.columns,
        "coef": res.params.values, "std_err": se, "p": pvals,
        "ci_low": ci_low, "ci_high": ci_high
    })
# ----------------- DiD wrappers -----------------

def run_did_set(panel, label=""):
    # --- SPEC A: share ~ post  (FE só de unidade; cluster em id)
    tabA = twfe_unit_singlecluster(panel, "share", ["post"])
    tabA["spec"] = f"A_{label}"

    # --- SPEC BG: share ~ post_gk + log_price + log_x  (2-way cluster)
    pan_bg = panel.copy()
    for g in range(1, 7):
        pan_bg[f"post_g{g}"] = pan_bg["post"] * (pan_bg["group"] == g).astype(int)
    x_bg = [f"post_g{g}" for g in range(1, 7)] + ["log_price", "log_x"]
    tabBG = twfe_2way_cluster(pan_bg, "share", xcols=x_bg, clus1="id", clus2="group")
    tabBG["spec"] = f"BG_{label}"

    # Ajustes múltiplos (Holm 5% e FDR-BH 10%) nos termos post_g*
    sub = tabBG[tabBG["term"].str.contains(r"post_g")].copy()
    rej_holm, p_holm, _, _ = multipletests(sub["p"].to_numpy(), alpha=0.05, method="holm")
    rej_fdr,  p_fdr,  _, _ = multipletests(sub["p"].to_numpy(), alpha=0.10, method="fdr_bh")
    sub["p_holm"] = p_holm; sub["reject_5pct_holm"] = rej_holm
    sub["p_fdr_bh_10pct"] = p_fdr; sub["reject_fdr10"] = rej_fdr

    # --- SPEC contínua: share ~ post + post_tax + log_price + log_x  (2-way cluster)
    # Nota: aqui dá para manter o 'post' como controle (absorve choque comum de tempo) e
    # o efeito diferencial vem de 'post_tax'
    tabC = twfe_2way_cluster(panel, "share", ["post", "post_tax", "log_price", "log_x"],
                             clus1="id", clus2="group")
    tabC["spec"] = f"Ccont_{label}"

    return tabA, tabBG, sub, tabC

# ----------------- QUAIDS checks -----------------

def read_quaids_params_vcov(path_xlsx):
    xls = pd.ExcelFile(path_xlsx)
    params = pd.read_excel(xls, "PARAMS")
    vcov_df = pd.read_excel(xls, "VCOV")
    def make_name(row):
        p = str(row["param"]).strip().lower()
        if p in ("alpha","beta","lambda") and not pd.isna(row.get("i", np.nan)):
            return f"{p}_{int(row['i'])}"
        if p=="gamma" and not pd.isna(row.get("i", np.nan)) and not pd.isna(row.get("j", np.nan)):
            return f"gamma_{int(row['i'])}_{int(row['j'])}"
        raise ValueError("Linha em PARAMS sem rótulo claro")
    params["name"] = params.apply(make_name, axis=1)
    theta = params.set_index("name")["est"].astype(float)
    if "param" in vcov_df.columns:
        vcov_df = vcov_df.set_index("param")
    commons = [n for n in theta.index if n in vcov_df.columns]
    theta = theta.loc[commons]
    V = vcov_df.loc[commons, commons].astype(float)
    return theta, V

def build_constraints(theta):
    G = 6
    names = theta.index.tolist(); p = len(names)
    def idx(n):
        return names.index(n) if n in names else None
    rows = []; qs = []
    # Adding-up
    r = np.zeros(p);  [r.__setitem__(idx(f"alpha_{i}"), 1.0) for i in range(1,G+1) if idx(f"alpha_{i}") is not None]
    rows.append(r); qs.append(1.0)
    r = np.zeros(p);  [r.__setitem__(idx(f"beta_{i}"),  1.0) for i in range(1,G+1) if idx(f"beta_{i}")  is not None]
    rows.append(r); qs.append(0.0)
    r = np.zeros(p);  [r.__setitem__(idx(f"lambda_{i}"),1.0) for i in range(1,G+1) if idx(f"lambda_{i}")is not None]
    rows.append(r); qs.append(0.0)
    # Somas coluna γ (∑i γ_ij = 0)
    for jg in range(1,G+1):
        r = np.zeros(p)
        for ig in range(1,G+1):
            j = idx(f"gamma_{ig}_{jg}")
            if j is not None: r[j]=1.0
        rows.append(r); qs.append(0.0)
    # Somas linha γ (homogeneidade: ∑j γ_ij = 0)
    for ig in range(1,G+1):
        r = np.zeros(p)
        for jg in range(1,G+1):
            j = idx(f"gamma_{ig}_{jg}")
            if j is not None: r[j]=1.0
        rows.append(r); qs.append(0.0)
    # Simetria γ_ij − γ_ji = 0
    for i in range(1,G+1):
        for jg in range(i+1, G+1):
            r = np.zeros(p)
            j1 = idx(f"gamma_{i}_{jg}"); j2 = idx(f"gamma_{jg}_{i}")
            if j1 is not None: r[j1] = 1.0
            if j2 is not None: r[j2] -= 1.0
            rows.append(r); qs.append(0.0)
    R = np.vstack(rows); q = np.array(qs)
    return R, q, names

def wald_test(R, q, theta, V):
    use_cols = np.where(np.abs(R).sum(axis=0) > 1e-12)[0]
    Rsub = R[:, use_cols]
    thetas = theta.values.reshape(-1,1)[use_cols]
    Vsub = V.values[np.ix_(use_cols, use_cols)]
    if np.isnan(Vsub).any():
        raise ValueError("VCOV tem NaN no subconjunto relevante das restrições.")
    diff = Rsub @ thetas - q.reshape(-1,1)
    RVRT = Rsub @ Vsub @ Rsub.T
    try:
        inv = np.linalg.inv(RVRT)
    except np.linalg.LinAlgError:
        inv = np.linalg.pinv(RVRT)
    W = float(diff.T @ inv @ diff); df = Rsub.shape[0]
    pval = float(chi2.sf(W, df))
    return W, df, pval

def run_quaids_checks(label, path_xlsx, base, price_cols, outdir: Path):
    theta, V = read_quaids_params_vcov(path_xlsx)
    R, q, names = build_constraints(theta)
    # Overall
    try:
        W_all, df_all, p_all = wald_test(R, q, theta, V)
        overall = pd.DataFrame([{"model":label, "test":"All (adding-up + homogeneity + symmetry)", "Wald":W_all, "df":df_all, "p_value":p_all}])
    except Exception as e:
        overall = pd.DataFrame([{"model":label, "test":"All (...)", "Wald":np.nan, "df":np.nan, "p_value":np.nan, "note":str(e)}])
    # Blocos
    G=6
    rows_add = list(range(0, 3+G))
    rows_hom = list(range(3+G, 3+2*G))
    rows_sym = list(range(3+2*G, R.shape[0]))
    def block(rows_idx, name):
        try:
            Wb, dfb, pb = wald_test(R[rows_idx,:], q[rows_idx], theta, V)
            return {"model":label, "test":name, "Wald":Wb, "df":dfb, "p_value":pb}
        except Exception as e:
            return {"model":label, "test":name, "Wald":np.nan, "df":np.nan, "p_value":np.nan, "note":str(e)}
    blocks = pd.DataFrame([
        block(rows_add, "Adding-up"),
        block(rows_hom, "Homogeneity"),
        block(rows_sym, "Symmetry")
    ])
    # Índices de preço: Stone vs Translog
    lnp = np.column_stack([np.log(np.where(base[c].astype(float).values>0, base[c].astype(float).values, np.nan)) for c in price_cols])
    lnp = np.nan_to_num(lnp, nan=0.0)
    w_pre = np.column_stack([base[f"w_despesahat{k}"].astype(float).values for k in range(1,7)])
    w_pre = np.nan_to_num(w_pre, nan=0.0)
    lnP_stone = (w_pre * lnp).sum(axis=1)
    alpha = np.array([theta.get(f"alpha_{i}", 0.0) for i in range(1,7)])
    Gamma = np.array([[theta.get(f"gamma_{i}_{j}", 0.0) for j in range(1,7)] for i in range(1,7)])
    lnP_tl = lnp @ alpha + 0.5 * np.sum((lnp @ Gamma) * lnp, axis=1)
    a = lnP_stone - lnP_stone.mean(); b = lnP_tl - lnP_tl.mean()
    rmse = float(np.sqrt(np.mean((a-b)**2)))
    corr = float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
    idx_comp = pd.DataFrame([{"model":label, "idx_a":"Stone", "idx_b":"Translog", "corr":corr, "rmse_demeaned":rmse}])
    # salvar
    overall.to_csv(outdir/f"quaids_wald_overall_{label}.csv", index=False)
    blocks.to_csv(outdir/f"quaids_wald_blocks_{label}.csv", index=False)
    idx_comp.to_csv(outdir/f"quaids_price_index_compare_{label}.csv", index=False)
    return overall, blocks, idx_comp

# ----------------- Bootstrap paramétrico + Placebos -----------------

def parametric_bootstrap(base, Wpre, Ppre, P1, P2, Xtot, M_hat, SD_hat, B, seed):
    rng = np.random.default_rng(seed)
    outA=[]; outBG=[]; outC=[]
    for b in range(1, B+1):
        Mdraw = rng.normal(M_hat, SD_hat)
        # cenário 1
        W1b = simulate_post_shares(Wpre, Ppre, P1, Xtot, Mdraw)
        pan1b = build_panel(base, Wpre, W1b, Ppre, P1, Xtot)
        A1b, BG1b, _, C1b = run_did_set(pan1b, f"s1_b{b}")
        outA.append(A1b.query("term=='post'")[["coef"]].assign(scenario=1, draw=b).rename(columns={"coef":"post_coef"}))
        outC.append(C1b.query("term=='post_tax'")[["coef"]].assign(scenario=1, draw=b).rename(columns={"coef":"post_tax_coef"}))
        gtab = BG1b[BG1b["term"].str.contains("post_g")][["term","coef"]].assign(scenario=1, draw=b)
        outBG.append(gtab)
        # cenário 2
        W2b = simulate_post_shares(Wpre, Ppre, P2, Xtot, Mdraw)
        pan2b = build_panel(base, Wpre, W2b, Ppre, P2, Xtot)
        A2b, BG2b, _, C2b = run_did_set(pan2b, f"s2_b{b}")
        outA.append(A2b.query("term=='post'")[["coef"]].assign(scenario=2, draw=b).rename(columns={"coef":"post_coef"}))
        outC.append(C2b.query("term=='post_tax'")[["coef"]].assign(scenario=2, draw=b).rename(columns={"coef":"post_tax_coef"}))
        gtab = BG2b[BG2b["term"].str.contains("post_g")][["term","coef"]].assign(scenario=2, draw=b)
        outBG.append(gtab)
    return (pd.concat(outA, ignore_index=True),
            pd.concat(outBG, ignore_index=True),
            pd.concat(outC, ignore_index=True))

def placebos_randomization(base, Wpre, Ppre, P1, P2, Xtot, M_hat, R, seed):
    rng = np.random.default_rng(seed + 12345)
    outA=[]; outBG=[]; outC=[]
    n = len(base)
    groups = list(range(1,7))
    for r in range(1, R+1):
        # permutar Δln p entre os 6 grupos dentro de cada domicílio (mantém distribuição, quebra vínculo grupo-choque)
        def permute_prices(Ppost):
            Ppost_perm = Ppost.copy()
            for i in range(n):
                perm = rng.permutation(groups)
                Ppost_perm.iloc[i,:] = Ppost.iloc[i, np.array(perm)-1].values
            return Ppost_perm
        # cenário 1
        P1p = permute_prices(P1)
        W1p = simulate_post_shares(Wpre, Ppre, P1p, Xtot, M_hat)
        pan1p = build_panel(base, Wpre, W1p, Ppre, P1p, Xtot)
        A1p, BG1p, _, C1p = run_did_set(pan1p, f"s1_pl{r}")
        outA.append(A1p.query("term=='post'")[["coef"]].assign(scenario=1, run=r).rename(columns={"coef":"post_coef"}))
        outC.append(C1p.query("term=='post_tax'")[["coef"]].assign(scenario=1, run=r).rename(columns={"coef":"post_tax_coef"}))
        gtab = BG1p[BG1p["term"].str.contains("post_g")][["term","coef"]].assign(scenario=1, run=r)
        outBG.append(gtab)
        # cenário 2
        P2p = permute_prices(P2)
        W2p = simulate_post_shares(Wpre, Ppre, P2p, Xtot, M_hat)
        pan2p = build_panel(base, Wpre, W2p, Ppre, P2p, Xtot)
        A2p, BG2p, _, C2p = run_did_set(pan2p, f"s2_pl{r}")
        outA.append(A2p.query("term=='post'")[["coef"]].assign(scenario=2, run=r).rename(columns={"coef":"post_coef"}))
        outC.append(C2p.query("term=='post_tax'")[["coef"]].assign(scenario=2, run=r).rename(columns={"coef":"post_tax_coef"}))
        gtab = BG2p[BG2p["term"].str.contains("post_g")][["term","coef"]].assign(scenario=2, run=r)
        outBG.append(gtab)
    return (pd.concat(outA, ignore_index=True),
            pd.concat(outBG, ignore_index=True),
            pd.concat(outC, ignore_index=True))