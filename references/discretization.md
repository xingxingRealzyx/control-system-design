# 离散化参考：采样率 / ZOH / Tustin / 差分方程

## 1. 采样率选择（先论证再取值）

| 依据 | 规则 |
|---|---|
| 闭环带宽 ω_b | 采样频率 ω_s = 2π/Ts ≥ 10×ω_b，保守取 20~30× |
| 主导极点/上升时间 | Ts ≤ t_r/10 ~ t_r/20 |
| 常用工程值 | 机械系统 100–1000 Hz；姿态控制 200–500 Hz；温度等慢过程 1–10 Hz |
| 硬件约束 | 受 ADC/总线/算法耗时限制时取能满足上述下限的最大可行 Ts |

报告中必须写出：由指标得 ω_b（或 t_r）→ 取 ω_s ≥ 20ω_b → Ts = xxx s（取整到硬件可行值）。

## 2. ZOH 精确离散化（被控对象标准方法）

$$A_d = e^{AT},\qquad B_d = \int_0^T e^{A\tau}d\tau\, B = A^{-1}(A_d - I)B$$

A 不可逆时**不要**用公式第二行，用增广矩阵指数一次得到（推荐，永远成立）：

$$\exp\!\left(\begin{bmatrix}A & B\\ 0 & 0\end{bmatrix}T\right) = \begin{bmatrix}A_d & B_d\\ 0 & I\end{bmatrix}$$

D 与 C 离散化不变（Cd = C, Dd = D，ZOH 假设 u[t_k, t_k+T) 恒定）。数值计算用 `linhelper.py c2d(A,B,T)`（scaling-and-squaring）。校验：Ad 特征值应等于 e^{λᵢT}。

## 3. 近似离散化方法对比（何时用哪种）

| 方法 | 映射 | 优点 | 缺点 | 用途 |
|---|---|---|---|---|
| 前向欧拉 | Ad ≈ I + AT | 最简单 | Ts 大时可能把稳定系统离散成不稳定 | 仅粗估 |
| 后向欧拉 | Ad ≈ (I−AT)⁻¹ | 无条件稳定 | 精度低 | 快速原型 |
| Tustin/双线性 | s ← (2/Ts)·(z−1)/(z+1) | 频率轴畸变有规律、稳态增益准 | 需防频率混叠 | **滤波器/补偿器离散首选** |
| Tustin 预畸变 | s ← (ω₀/Ts)·cot(ω₀Ts/2)·(z−1)/(z+1) | 在频率 ω₀ 处精确匹配 | 只匹配单点 | 谐振/陷波器 |

混叠提醒：ω_s ≥ 20ω_b 基本免疫；若对象有显著高于 ω_s/2 的谐振，需加抗混叠模拟滤波（二阶 Butterworth，截止 ≈ ω_s/4）。

## 4. 矩阵指数数值算法（linhelper.py 实现，勿另造）

Scaling-and-squaring：取 s = max(0, ceil(log2‖A‖T))，A ← A·T/2^s，用 Padé(13) 或泰勒级数（至项 < 1e-16）算 e^A，再平方 s 次还原。误差受浮点精度限制 ~1e-15。

## 5. 控制器离散实现（Stage 7 核心产出：差分方程）

### 5.1 PID 位置式（默认用这个，配 clamping 抗饱和）

$$u_I[k] = u_I[k{-}1] + K_i T_s e[k]\quad(\text{含饱和禁止条件})$$
$$u_D[k] = \frac{K_d}{T_s + K_d/N}\Big(u_D[k{-}1] + T_s N\,(-\dot y[k])\Big)\ \text{，或直接}\ u_D[k] = \frac{K_d}{N T_s}\cdot\frac{(N T_s)(-y[k]) - \text{(滤波历史)}}{1 + \cdot}$$

标准一阶滤波微分（推荐写法，令 τ_d = Kd/N）：

$$d_f[k] = \frac{\tau_d}{\tau_d + T_s}d_f[k{-}1] + \frac{K_d}{\tau_d + T_s}\big(y[k] - y[k{-}1]\big)$$

微分取 −y（见 controllers.md §1.3），输出 u = Kp·e + u_I + d_f，再限幅。

### 5.2 PID 增量式（执行器接受增量/步进电机时用）

$$\Delta u[k] = K_p(e[k]-e[k-1]) + K_i T_s\,e[k] + K_d\frac{e[k]-2e[k-1]+e[k-2]}{T_s}$$
$$u[k] = \mathrm{sat}\big(u[k-1] + \Delta u[k]\big)$$

抗饱和天然内建（限幅在累加处）；缺点是积分基准不显式，切换手动/自动需预置 u[k−1]。

### 5.3 状态反馈 + 观测器

$$u[k] = -K\hat x[k]\ (+\ N_r r[k])$$
$$\hat x^-[k] = A_d\hat x[k{-}1] + B_d u[k{-}1],\qquad \hat x[k] = \hat x^- + L_d\big(y[k] - C_d\hat x^-[k]\big)$$

Ld 由离散观测器极点（|λ|≈e^{λ_c T}）或稳态 Kalman 增益离线求出。执行顺序（每拍）：采样 y → 更新观测器 → 算 u → 输出保持。

LQR/极点配置设计出的 K 直接用于离散（增益不变）；若连续 LQR 的 Ts 较大（>1/20 带宽），改用 DLQR 重算 K_d 保证最优性。

### 5.4 MPC

天然离散（§5.1 controllers.md 的 Δu 模型即离散域）。每拍：滚动 ξ[k] = [x̂[k]; u[k−1]] → 解 QP → 施加 Δu[0]。

### 5.5 滤波与补偿器离散化

连续传函 G(s) 先转成可控标准型状态空间再 ZOH，或 Tustin 逐环节变换展开为差分方程。陷波/低通直接用双线性模板。

## 6. 抗饱和离散伪码（与 5.1 配套）

```
u_unsat = Kp*e + I + D_f
if (u_unsat > u_max and e > 0) or (u_unsat < u_min and e < 0):
    pass                    # 本拍不积分
else:
    I = I + Ki*Ts*e
I  = clamp(I, I_min, I_max)  # 双保险
u  = clamp(Kp*e + I + D_f, u_min, u_max)
```

## 7. 常见坑清单（离散化阶段逐条对照）

- [ ] 代数环：D ≠ 0（或 K 作用于含当前拍输出的量）且一拍内反馈 → 检查信号时序，控制器输出只用本拍及以前的数据。
- [ ] 微分噪声：未加滤波（N≥8）或对 r 微分。
- [ ] 积分饱和：无 clamping/限幅。
- [ ] 单位不一致：Ts 单位 s，带宽 rad/s，别混 Hz。
- [ ] 离散稳定性未校验：eig(Ad) 全部 |λ|<1。
- [ ] ZOH 假设破坏：控制器输出在拍间被改写（如 PWM 中断重写占空比）——保持行为写入代码注释。
- [ ] Ts 与实际定时漂移：报告注明 Ts 必须由硬件定时器保证。
