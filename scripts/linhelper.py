#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""控制系统分析辅助工具（运动学/动力学 skill 专用）

仅依赖 numpy；安装 scipy 后自动加速（ARE 方程与 expm 用官方实现）。
用法：
    python3 linhelper.py demo          # 倒立摆端到端演示
    或在脚本中 import:
        from linhelper import rank_report, lqr, dlqr, c2d, place, step_metrics, freq_margins
"""
import numpy as np

try:
    import scipy.linalg as _sla
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


# ---------------------------------------------------------------- 矩阵指数
def expm(M):
    """矩阵指数。scipy 可用则用 scipy，否则 scaling-and-squaring + 泰勒级数。"""
    if _HAVE_SCIPY:
        return _sla.expm(M)
    M = np.asarray(M, dtype=float)
    n = M.shape[0]
    norm = np.max(np.sum(np.abs(M), axis=1)) if n else 0.0
    s = max(0, int(np.ceil(np.log2(norm + 1e-300)))) + 1
    A = M / (2.0 ** s)
    E = np.eye(n)
    term = np.eye(n)
    for k in range(1, 40):
        term = term @ A / k
        E = E + term
        if np.max(np.abs(term)) < 1e-18:
            break
    for _ in range(s):
        E = E @ E
    return E


# ---------------------------------------------------------------- 能控能观
def ctrb(A, B):
    A = np.atleast_2d(np.asarray(A, float))
    B = np.asarray(B, float)
    if B.ndim == 1:
        B = B.reshape(-1, 1)
    n = A.shape[0]
    return np.hstack([np.linalg.matrix_power(A, i) @ B for i in range(n)])


def obsv(A, C):
    A = np.atleast_2d(np.asarray(A, float))
    C = np.asarray(C, float)
    if C.ndim == 1:
        C = C.reshape(1, -1)
    n = A.shape[0]
    return np.vstack([C @ np.linalg.matrix_power(A, i) for i in range(n)])


def pbh_uncontrollable_modes(A, B):
    """返回不可控特征值列表（PBH 判据）。"""
    n = A.shape[0]
    Bm = np.asarray(B, float).reshape(n, -1)
    bad = []
    for lam in np.linalg.eigvals(A):
        M = np.hstack([lam * np.eye(n) - A, Bm])
        if np.linalg.matrix_rank(M, tol=1e-8) < n:
            if not any(np.abs(lam - u) < 1e-6 for u in bad):
                bad.append(lam)
    return bad


def pbh_unobservable_modes(A, C):
    eig = np.linalg.eigvals(A)
    n = A.shape[0]
    Cm = np.asarray(C, float).reshape(-1, n)
    bad = []
    for lam in eig:
        M = np.vstack([lam * np.eye(n) - A, Cm])
        if np.linalg.matrix_rank(M, tol=1e-8) < n:
            if not any(np.abs(lam - u) < 1e-6 for u in bad):
                bad.append(lam)
    return bad


def rank_report(A, B, C, name="系统"):
    """一次性给出：极点、能控性、能观性、PBH 坏模态、能稳/能检结论。"""
    A = np.atleast_2d(np.asarray(A, float))
    n = A.shape[0]
    poles = np.linalg.eigvals(A)
    rc = np.linalg.matrix_rank(ctrb(A, B))
    ro = np.linalg.matrix_rank(obsv(A, C))
    unctl = pbh_uncontrollable_modes(A, B)
    unobs = pbh_unobservable_modes(A, C)

    print(f"===== {name} 分析报告 =====")
    print("极点:")
    for p in sorted(poles, key=lambda z: (-z.real, z.imag)):
        print(f"   {p:.4g}")
    print(f"能控性矩阵秩: {rc}/{n}  -> {'完全能控' if rc == n else '不完全能控'}")
    if unctl:
        print(f"   不可控模态(PBH): {[f'{p:.4g}' for p in unctl]}")
        stab = all(u.real < 0 for u in unctl)
        print(f"   {'可稳（不可控模态均稳定）' if stab else '不能稳！物理上该执行器无法镇定系统'}")
    print(f"能观性矩阵秩: {ro}/{n}  -> {'完全能观' if ro == n else '不完全能观'}")
    if unobs:
        print(f"   不可观模态(PBH): {[f'{p:.4g}' for p in unobs]}")
        det = all(u.real < 0 for u in unobs)
        print(f"   {'能检（不可观模态均稳定）' if det else '不能检，无法构造收敛观测器'}")
    print("=" * 30)
    return {"poles": poles, "rank_ctlb": rc, "rank_obsv": ro,
            "uncontrollable": unctl, "unobservable": unobs}


# ---------------------------------------------------------------- LQR / DLQR
def _lyapunov_solve(Aa, M):
    """解 AaᵀP + P·Aa = −M（vec + kron，行主序约定）。"""
    n = Aa.shape[0]
    coef = np.kron(Aa.T, np.eye(n)) + np.kron(np.eye(n), Aa.T)
    p = np.linalg.solve(coef, -M.reshape(-1))
    return p.reshape(n, n)


def _initial_stabilizing_gain(A, B):
    """Kleinman 迭代的初始镇定增益：Ackermann 配置全 −α 重极点，α 逐档加大。"""
    n, m = A.shape[0], B.shape[1]
    b1 = B[:, :1]
    for alpha in (0.5, 1.0, 2.0, 4.0, 8.0):
        try:
            k1 = np.asarray(place(A, b1, np.full(n, -alpha, dtype=complex))).reshape(1, -1)
        except Exception:
            continue
        if np.max(np.real(np.linalg.eigvals(A - b1 @ k1))) < 0:
            K0 = np.zeros((m, n))
            K0[0] = k1[0]
            return K0
    raise RuntimeError("找不到初始镇定增益（第一输入列可能不完全能控），请安装 scipy: pip3 install scipy")


def lqr(A, B, Q, R):
    """连续代数Riccati方程。scipy 可用用官方实现，否则 Kleinman 迭代。
    返回 (K, P, closed_loop_eigs)。"""
    A, B, Q, R = map(lambda x: np.atleast_2d(np.asarray(x, float)), (A, B, Q, R))
    if _HAVE_SCIPY:
        P = _sla.solve_continuous_are(A, B, Q, R)
    else:
        K0 = _initial_stabilizing_gain(A, B)
        P = _lyapunov_solve(A - B @ K0, Q + K0.T @ R @ K0)
        for _ in range(100):
            K = np.linalg.solve(R, B.T @ P)
            Pn = _lyapunov_solve(A - B @ K, Q + K.T @ R @ K)
            if np.max(np.abs(Pn - P)) < 1e-10:
                P = Pn
                break
            P = Pn
    K = np.linalg.solve(R, B.T @ P)
    return K, P, np.linalg.eigvals(A - B @ K)


def dlqr(Ad, Bd, Q, R):
    """离散代数Riccati方程（值迭代，稳定条件弱依赖，scipy 可用则官方实现）。
    返回 (Kd, P, closed_loop_eigs)。"""
    Ad, Bd, Q, R = map(lambda x: np.atleast_2d(np.asarray(x, float)), (Ad, Bd, Q, R))
    if _HAVE_SCIPY:
        P = _sla.solve_discrete_are(Ad, Bd, Q, R)
    else:
        P = Q.copy()
        for _ in range(2000):
            BtPB = Bd.T @ P @ Bd
            K = np.linalg.solve(R + BtPB, Bd.T @ P @ Ad)
            Pn = Q + Ad.T @ P @ (Ad - Bd @ K)
            if np.max(np.abs(Pn - P)) < 1e-11:
                P = Pn
                break
            P = Pn
    K = np.linalg.solve(R + Bd.T @ P @ Bd, Bd.T @ P @ Ad)
    return K, P, np.linalg.eigvals(Ad - Bd @ K)


# ---------------------------------------------------------------- 极点配置
def place(A, B, poles):
    """Ackermann 公式极点配置（SISO；MIMO 请取一列能控输入）。
    返回 K 使 eig(A - B K) = poles。"""
    A = np.atleast_2d(np.asarray(A, float))
    b = np.asarray(B, float).reshape(-1, 1) if np.asarray(B).ndim == 1 else np.asarray(B, float)[:, :1]
    n = A.shape[0]
    if np.linalg.matrix_rank(ctrb(A, b)) != n:
        raise ValueError("系统不完全能控，无法任意配置极点（见 rank_report）")
    coeffs = np.poly(np.asarray(poles, complex).real if np.max(np.abs(np.asarray(poles, complex).imag)) < 1e-9
                     else np.asarray(poles, complex))
    # phi(A) = A^n + c1 A^{n-1} + ... + cn I（Horner）
    phiA = np.zeros((n, n))
    for c in coeffs:
        phiA = phiA @ A + c * np.eye(n)
    e_n = np.zeros(n)
    e_n[-1] = 1.0
    K = e_n @ np.linalg.inv(ctrb(A, b)) @ phiA
    # 返回形状与 B 的输入形状兼容，使 A - B @ K 直接可算
    return K if np.asarray(B, float).ndim == 1 else K.reshape(1, -1)


# ---------------------------------------------------------------- 离散化
def c2d(A, B, Ts):
    """ZOH 精确离散化（增广矩阵指数法）。返回 (Ad, Bd)。"""
    A = np.atleast_2d(np.asarray(A, float))
    B = np.asarray(B, float)
    if B.ndim == 1:
        B = B.reshape(-1, 1)
        squeeze = True
    else:
        squeeze = False
    n, m = A.shape[0], B.shape[1]
    Z = np.zeros((n + m, n + m))
    Z[:n, :n] = A
    Z[:n, n:] = B
    E = expm(Z * Ts)
    Ad, Bd = E[:n, :n], E[:n, n:]
    return (Ad, Bd.reshape(-1)) if squeeze else (Ad, Bd)


# ---------------------------------------------------------------- 指标计算
def step_metrics(t, y, tol=0.02):
    """阶跃响应指标：超调%、峰值时间、上升时间(10-90%)、调节时间(默认2%)、稳态值。"""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    yf = float(np.mean(y[-max(1, len(y) // 20):]))
    scale = abs(yf) if abs(yf) > 1e-9 else 1.0
    i_pk = int(np.argmax(y))
    Mp = max(0.0, (y[i_pk] - yf) / scale) * 100.0
    tp = float(t[i_pk])

    def _cross(level):
        for i in range(1, len(y)):
            if (y[i-1] - level) * (y[i] - level) <= 0 and y[i] != y[i-1]:
                frac = (level - y[i-1]) / (y[i] - y[i-1])
                return float(t[i-1] + frac * (t[i] - t[i-1]))
        return None

    tr = None
    if abs(yf) > 1e-9:
        t10, t90 = _cross(0.1 * yf), _cross(0.9 * yf)
        if t10 is not None and t90 is not None:
            tr = abs(t90 - t10)
    outside = np.abs(y - yf) > tol * scale
    ts = float(t[np.max(np.nonzero(outside)[0]) + 1]) if outside.any() else float(t[0])
    return {"overshoot_pct": Mp, "peak_time": tp, "rise_time": tr,
            "settling_time": ts, "steady_state": yf}


def freq_margins(A, B, C, D=None, wmin=1e-3, wmax=1e4, npts=20001):
    """SISO 频域裕度扫描。返回 dict(OneGainCross_wc, phase_margin_deg,
    PhaseCross_wg, gain_margin_db)。"""
    A = np.atleast_2d(np.asarray(A, float))
    B = np.asarray(B, float).reshape(-1, 1)
    C = np.asarray(C, float).reshape(1, -1)
    D = np.zeros((1, 1)) if D is None else np.asarray(D, float).reshape(1, 1)
    n = A.shape[0]
    w = np.logspace(np.log10(wmin), np.log10(wmax), npts)
    mag = np.empty(npts)
    ph = np.empty(npts)
    for i, wi in enumerate(w):
        G = C @ np.linalg.solve(1j * wi * np.eye(n) - A, B) + D
        mag[i] = np.abs(G[0, 0])
        ph[i] = np.degrees(np.angle(G[0, 0]))
    mag_db = 20 * np.log10(mag)
    res = {"wc": None, "phase_margin_deg": None, "wg": None, "gain_margin_db": None}
    # 相位裕度：|G|=1 处
    idx = np.where(np.diff(np.sign(mag_db)) < 0)[0]
    if len(idx):
        i = idx[0]
        wc = np.interp(0.0, mag_db[i:i+2][::-1], w[i:i+2][::-1])
        pm = 180.0 + np.interp(np.log(wc), np.log(w), ph)
        res["wc"], res["phase_margin_deg"] = float(wc), float(pm)
    # 增益裕度：相位 -180° 处
    phw = (ph + 180) % 360 - 180
    idx = np.where(np.diff(np.sign(phw + 180.0)) != 0)[0]
    for i in idx:
        wg = float(np.sqrt(w[i] * w[i + 1]))
        gmdB = float(-20 * np.log10(np.interp(wg, w, mag)))
        if gmdB > 0:
            res["wg"], res["gain_margin_db"] = wg, gmdB
            break
    return res


# ---------------------------------------------------------------- 演示
def _demo():
    """小车倒立摆端到端（modeling.md 算例 C 数值）。"""
    M, m, l, g = 1.0, 0.1, 0.5, 9.81
    I = m * l * l / 3.0
    D = (M + m) * (I + m * l * l) - (m * l) ** 2
    A = np.array([
        [0, 1, 0, 0],
        [0, 0, -m * m * g * l * l / D, (I + m * l * l) / D],
        [0, 0, 0, 1],
        [0, 0, (M + m) * m * g * l / D, 0],
    ])
    B = np.array([[0.0], [(I + m * l * l) / D], [0.0], [-m * l / D]])
    C = np.array([[1.0, 0, 0, 0]])
    rep = rank_report(A, B, C, "倒立摆")
    assert rep["rank_ctlb"] == 4 and rep["rank_obsv"] == 4, "能控能观检验失败"

    Q = np.diag([1.0, 1.0, 100.0, 10.0])   # Bryson 后手动加权摆角
    R = np.array([[1.0]])
    K, P, eig_cl = lqr(A, B, Q, R)
    print("LQR K =", np.round(K, 4))
    print("闭环极点 =", np.round(eig_cl, 4))
    assert np.max(np.real(eig_cl)) < 0, "LQR 闭环应稳定"

    pd = np.array([-4, -8, -8, -16], dtype=complex)
    Kp = place(A, B, pd)
    print("Ackermann K =", np.round(Kp, 4),
          " 实际闭环极点 =", np.sort_complex(np.linalg.eigvals(A - B @ Kp)))

    Ts = 0.01
    Ad, Bd = c2d(A, B, Ts)
    # ZOH 一致性校验：eig(Ad) 应等于 e^{λT}（开环不稳定对象模态在单位圆外是正常的）
    lam = np.sort_complex(np.linalg.eigvals(A))
    lam_d = np.sort_complex(np.linalg.eigvals(Ad))
    err = np.max(np.abs(lam_d - np.exp(lam * Ts)))
    print("ZOH 校验 max|eig(Ad) − e^(λT)| =", err)
    assert err < 1e-6, "ZOH 离散化与连续特征值不一致"
    Kd, _, eig_dd = dlqr(Ad, Bd, Q, R)
    print("DLQR Kd =", np.round(Kd, 4), " 闭环模态 =", np.round(np.abs(eig_dd), 4))

    t = np.linspace(0, 5, 2001)
    y = 1 - np.exp(-2 * t) * (np.cos(3 * t) + 0.5 * np.sin(3 * t))
    print("step_metrics 示例:", {k: (round(v, 4) if v is not None else None)
                                  for k, v in step_metrics(t, y).items()})
    print("scipy:", "可用（ARE/expm 用官方实现）" if _HAVE_SCIPY else "不可用（Kleinman 回退）")


if __name__ == "__main__":
    _demo()
