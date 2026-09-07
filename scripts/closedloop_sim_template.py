#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""非线性闭环仿真骨架（skill Stage 8 标准工具）

结构：RK4 积分连续非线性被控对象 + 每个采样周期 Ts 调用一次离散控制器（ZOH）。
用法：复制本文件改名，替换标有 [REPLACE] 的三处（对象动力学、控制器、参考），
运行后打印 step 指标表并可导出对拍 CSV。仅依赖 numpy。
"""
import numpy as np

from linhelper import step_metrics

# ================================================================ [REPLACE 1/3] 被控对象
# 以小车倒立摆为例（modeling.md 算例 C，θ 为与竖直向上夹角）。
PARAMS = dict(M=1.0, m=0.1, l=0.5, g=9.81, I=0.1 * 0.5**2 / 3.0)


def plant_f(x, u):
    """非线性动力学 ẋ = f(x, u)。x = [x, ẋ, θ, θ̇]。"""
    M, m, l, g, I = (PARAMS[k] for k in ("M", "m", "l", "g", "I"))
    px, vx, th, dth = x
    s, c = np.sin(th), np.cos(th)
    D = (M + m) * (I + m * l * l) - (m * l) ** 2
    # 由拉格朗日方程联立消元的显式非线性解：
    # [(M+m), ml c; ml c, (I+ml²)] [ẍ; θ̈] = [F + ml θ̇² s; mgl s]
    a11, a12 = M + m, m * l * c
    a21, a22 = m * l * c, I + m * l * l
    b1 = u + m * l * dth * dth * s
    b2 = m * g * l * s
    det = a11 * a22 - a12 * a21
    xdd = (a22 * b1 - a12 * b2) / det
    thdd = (-a21 * b1 + a11 * b2) / det
    return np.array([vx, xdd, dth, thdd])


PLANT_X0 = np.array([0.0, 0.0, np.radians(10.0), 0.0])  # 初始条件
STATE_NAMES = ["x [m]", "xdot [m/s]", "theta [rad]", "thetadot [rad/s]"]
OUTPUT_IDX = [0, 1, 2, 3]    # 测量：本演示全状态可测；输出反馈时改为此处子集并替换观测器
CTRL_OUT_IDX = 2             # 报告指标所关注的主输出分量（摆角）


def rk4(x, u, h):
    k1 = plant_f(x, u)
    k2 = plant_f(x + 0.5 * h * k1, u)
    k3 = plant_f(x + 0.5 * h * k2, u)
    k4 = plant_f(x + h * k3, u)
    return x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


# ================================================================ [REPLACE 2/3] 离散控制器
# 以 LQR 增益为例（由 linhelper 离线算出；生成代码后此处换成对拍过的差分方程实现）。
from linhelper import lqr  # noqa: E402

_A = np.array([[0, 1, 0, 0],
               [0, 0, -0.7178, 0.9756],
               [0, 0, 0, 1],
               [0, 0, 15.7917, 0]])
_B = np.array([[0.0], [0.9756], [0.0], [-1.4634]])
K, _, _ = lqr(_A, _B, np.diag([1.0, 1.0, 100.0, 10.0]), np.array([[1.0]]))


class Controller:
    """离散控制器实例：每个 Ts 调一次 update。u = −Kx（演示版，全状态反馈）。"""

    def __init__(self):
        self.K = K.ravel().copy()
        self.u_min, self.u_max = -20.0, 20.0     # 执行器限幅 [N]
        self.u = 0.0

    def update(self, y_meas, r=0.0):
        x_est = np.array(y_meas)                  # [REPLACE] 输出不可全测时换成观测器差分方程
        u = -self.K @ x_est + 0.0 * r
        self.u = float(np.clip(u, self.u_min, self.u_max))
        return self.u


# ================================================================ 主循环
def simulate(Ts=0.01, T_end=8.0):
    ctrl = Controller()
    x = PLANT_X0.copy()
    n_sub = 20                                    # 每个控制拍内 RK4 子步（ plant 精积分）
    h = Ts / n_sub
    t_arr, x_arr, u_arr = [0.0], [x.copy()], [0.0]
    n_steps = int(round(T_end / Ts))
    for k in range(n_steps):
        y_meas = x[OUTPUT_IDX]                    # 传感器采样
        u = ctrl.update(y_meas)                   # 控制器一拍
        for _ in range(n_sub):                    # ZOH 下积分对象
            x = rk4(x, u, h)
        t_arr.append((k + 1) * Ts)
        x_arr.append(x.copy())
        u_arr.append(u)                           # u_arr[k] 为 [k-1,k) 保持的控制量
    return np.array(t_arr), np.array(x_arr), np.array(u_arr)


def main():
    t, X, U = simulate()
    # 镇定问题：初值→0 是负向过程，取"初值−当前"正向化后再统计指标
    y = X[0, CTRL_OUT_IDX] - X[:, CTRL_OUT_IDX]
    m = step_metrics(t, y)
    print("===== 非线性闭环仿真结果 =====")
    print(f"{'指标':<14}{'数值':>12}")
    for key in ("overshoot_pct", "peak_time", "rise_time", "settling_time", "steady_state"):
        v = m[key]
        print(f"{key:<14}{(round(v, 4) if v is not None else 'N/A'):>12}")
    print(f"控制量峰值     {np.max(np.abs(U)):>10.4f}  (限幅 ±{Controller().u_max})")
    print(f"末态各状态: " + ", ".join(f"{n}={v:.4g}" for n, v in zip(STATE_NAMES, X[-1])))
    # 对拍 CSV：每行 t, y..., u  —— 供 C harness 逐拍回放（codegen.md §6）
    np.savetxt("sim_golden.csv", np.hstack([t[:, None], X[:, OUTPUT_IDX], U[:, None]]),
               delimiter=",", header="t," + ",".join(f"y{i}" for i in OUTPUT_IDX) + ",u", comments="")
    print("已导出 sim_golden.csv")
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax[0].plot(t, np.degrees(X[:, 2]), label="theta [deg]")
        ax[0].plot(t, X[:, 0], label="x [m]")
        ax[0].legend(), ax[0].grid(True), ax[0].set_ylabel("states")
        ax[1].plot(t, U, label="u [N]")
        ax[1].legend(), ax[1].grid(True), ax[1].set_xlabel("t [s]")
        fig.savefig("closedloop_sim.png", dpi=120)
        print("已保存 closedloop_sim.png")
    except ImportError:
        print("（无 matplotlib，跳过画图）")


if __name__ == "__main__":
    main()
