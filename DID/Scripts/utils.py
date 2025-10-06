import numpy as np
import pandas as pd

def make_log(arr):
    arr = np.asarray(arr, dtype=float)
    out = np.log(np.where(arr>0, arr, np.nan))
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)

def read_elas_matrix_with_ci(df):
    """Lê a aba de elasticidades marshallianas (status quo) e retorna matriz 6x6 + desvios-padrão aproximados."""
    cols = {c.strip().lower(): c for c in df.columns}
    est_key = next((cols[k] for k in cols if k in ("est","estimate","coef","beta")), None)
    i_key   = next((cols[k] for k in cols if k in ("i","good_i","row","g_i")), None)
    j_key   = next((cols[k] for k in cols if k in ("j","good_j","col","g_j")), None)
    lo_key  = next((cols[k] for k in cols if k in ("lo","low","lwr","ci_low","lower")), None)
    hi_key  = next((cols[k] for k in cols if k in ("hi","high","upr","ci_high","upper")), None)
    if est_key is None or i_key is None or j_key is None:
        raise ValueError("Aba de elasticidades precisa de colunas i, j, est")
    est = np.zeros((6,6)); sd = np.zeros((6,6))
    for _, row in df.iterrows():
        i = int(row[i_key]) - 1; j = int(row[j_key]) - 1
        est[i,j] = float(row[est_key])
        if lo_key is not None and hi_key is not None and pd.notna(row[lo_key]) and pd.notna(row[hi_key]):
            sd[i,j] = float(row[hi_key] - row[lo_key])/(2*1.96)
        else:
            sd[i,j] = max(abs(est[i,j])*0.10, 1e-6)  # fallback conservador
    return est, sd

def simulate_post_shares(Wpre, Ppre, Ppost, Xtot, M_pre, dlogx=0.0, inc_elast=None):
    """Simula shares pós com elasticidades Marshall do status quo (QUAIDS), opcional dln x (via inc_elast)."""
    n, G = len(Wpre), 6
    p0 = Ppre.to_numpy();  w0 = Wpre.to_numpy();  x0 = Xtot.reshape(-1,1)
    q0 = (w0 * x0) / np.where(p0==0, np.nan, p0)
    q0 = np.nan_to_num(q0, nan=0.0, posinf=0.0, neginf=0.0)

    p1 = Ppost.to_numpy()
    dlp = np.log(np.where(p1>0, p1, np.nan)) - np.log(np.where(p0>0, p0, np.nan))
    dlp = np.nan_to_num(dlp, nan=0.0, posinf=0.0, neginf=0.0)

    dlnq = dlp @ M_pre.T
    if inc_elast is not None and dlogx != 0.0:
        dlnq += dlogx * inc_elast.reshape(1,-1)

    q1 = q0 * np.exp(dlnq)
    num = p1 * q1; den = num.sum(axis=1).reshape(-1,1)
    w1 = np.divide(num, np.where(den==0, np.nan, den))
    w1 = np.nan_to_num(w1, nan=1/6, posinf=1/6, neginf=1/6)
    w1 = w1 / w1.sum(axis=1, keepdims=True)
    return pd.DataFrame(w1, columns=[f"w{g}" for g in range(1,7)])

def build_panel(base, Wpre, Wpost, Ppre, Ppost, Xtot):
    ids = base["id"].to_numpy(); rows = []
    for g in range(1,7):
        rows.append(pd.DataFrame({
            "id": ids, "group": g, "time": 0, "post": 0,
            "share": Wpre[f"w{g}"].to_numpy(),
            "log_price": make_log(Ppre[f"p{g}"].to_numpy()),
            "log_x": make_log(Xtot)}))
        rows.append(pd.DataFrame({
            "id": ids, "group": g, "time": 1, "post": 1,
            "share": Wpost[f"w{g}"].to_numpy(),
            "log_price": make_log(Ppost[f"p{g}"].to_numpy()),
            "log_x": make_log(Xtot)}))
    pan = pd.concat(rows, ignore_index=True)
    pan["unit_group_id"] = pan["id"].astype(str) + "_" + pan["group"].astype(int).astype(str)
    # Δln p por id×grupo (achata 2D→1D com ordem p{1..6} por id; repetimos nas 2 linhas (pre e post) para conveniência)
    dlp = (np.log(np.where(Ppost.to_numpy()>0, Ppost.to_numpy(), np.nan)) -
           np.log(np.where(Ppre.to_numpy()>0, Ppre.to_numpy(), np.nan)))
    dlp = np.nan_to_num(dlp, nan=0.0)
    pan["tax_shock"] = np.tile(dlp.flatten(), 2)
    pan["post_tax"] = pan["post"] * pan["tax_shock"]
    return pan