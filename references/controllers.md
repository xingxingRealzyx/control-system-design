# 控制器设计参考：PID / 极点配置 / LQR / 观测器 / MPC / 进阶

## 0. 选型决策矩阵（详细版）

| 方法 | 系统类型 | 需要的模型信息 | 优势 | 局限 | 计算量(在线) |
|---|---|---|---|---|---|
| PID | SISO | 只需整定 | 结构简单、工业标准、免模型 | 多变量/强耦合无力；性能靠试 | 极低 |
| 极点配置 | MIMO 全状态可测 | A,B | 指标→极点直接映射 | 无最优性；对饱和敏感 | 低 |
| LQR | MIMO 全状态可测 | A,B,Q,R | 最优折中、固有裕度(≥60°PM/∞GM)、调参只有 Q,R | 无显式超调约束；需全状态 | 低（K 离线算好） |
| LQG = LQR+KF | 部分状态可测 | +噪声统计 | 输出反馈最优 | 裕度保证失效，需验证 | 中 |
| MPC | MIMO + 约束 | 离散模型 | 显式处理约束、可预见参考、每步重优化 | 需在线解 QP、调参复杂 | 高 |

选择顺序：SISO 简单对象 → PID；指标用"能量/控制量折中"表述或 MIMO → LQR；用户明确提约束/预测/滚动优化 → MPC；输出不可全测 → 加观测器（任何上述方法都可配观测器，分离原理）。

## 1. PID（经典）

### 1.1 形式与参数意义

$$u(t) = K_p e(t) + K_i\int_0^t e\,d\tau + K_d \dot e(t)$$

Kp 提刚度（加快、升超调）、Ki 消稳态误差（降低稳定裕度）、Kd 抑制超调（放大噪声）。规范型 Ti = Ki/Kp、Td = Kd/Kp。

### 1.2 整定路线

**A. 临界比例法（Ziegler–Nichols）**：纯 P 闭环，增 Kp 至等幅振荡，记 Ku 与振荡周期 Pu：

| 类型 | Kp | Ti | Td |
|---|---|---|---|
| P | 0.5Ku | — | — |
| PI | 0.45Ku | Pu/2 | — |
| PID | 0.6Ku | Pu/2 | Pu/8 |

Z-N 结果偏激进（超调大），取 0.5×Kp 作起点微调。

**B. 极点配置法（推荐，和模型挂钩）**：把 (Kd s²+Kp s+Ki) 与对象构成的开环闭环特征多项式，与期望 (s²+2ζωₙs+ωₙ²)(s+p₀) 比较系数，解出 Kp,Ki,Kd。期望极点由指标反推：ωₙ ≈ 4/(ζ·t_s)，ζ 由超调反查。

**C. 频域法**：调 Kp 使穿越频率到目标带宽，补偿相位使 PM≥50°；Kd 提供超前相位（最大 +90° 单个、双领先最多 ~+120°），Ki 兜底低频增益。

### 1.3 工程必备结构（缺一不可，写入设计与代码）

1. **微分对测量不对参考**：微分项作用于 −y 而非 r，避免设定值阶跃微分冲击。
2. **微分滤波**：D(s) = Kd s/(1 + s·Td/N)，N = 8~20，噪声大的传感器取小 N。
3. **抗积分饱和（clamping）**：

```
u_unsat = u_P + u_I + u_D
if 饱和 且 sign(e) 使积分继续加深:
    本拍不积分（跳过 u_I += Ki·e·Ts）
else:
    u_I += Ki·e·Ts
```

4. **输出限幅与速率限幅**（对执行器而言），限幅值写入配置。
5. 积分限幅（双保险）：|u_I| ≤ u_max。

### 1.4 串级 PID（机械/飞行控制常用）

外环(位置/姿态角)输出作为内环(速度/角速度)参考，带宽比 ≥ 3~5×（外环慢、内环快），逐环由内向外整定。四旋翼典型三级串级：位置→姿态→角速度。

## 2. 极点配置 + Ackermann（现代经典）

**Ackermann 公式**（完全能控前提）：

$$u = -Kx,\qquad K = \begin{bmatrix}0&\cdots&0&1\end{bmatrix}\mathcal{C}^{-1}\,\varphi_d(A)$$

φ_d(s) 为期望特征多项式，φ_d(A) 为矩阵多项式。期望极点选取：主导二阶对 (ζ,ωₙ) 按指标反推，其余实极点放在主导对的 3~5 倍远处；n≥4 时用 Butterworth 模板（单位 ωₙ）：

n=2: s²+1.414s+1；n=3: s³+2s²+2s+1；n=4: s⁴+2.613s³+3.414s²+2.613s+1。整体系数乘 ωₙ 的对应次幂。

**参考跟踪两条路**（任选其一，推荐增广）：

1. **前馈增益**（快但不消扰动）：u = −Kx + N̄r，SISO 时 N̄ = −1/[C(A−BK)⁻¹B]。
2. **积分器增广**（无静差，抗常值扰动）：

$$\tilde x = \begin{bmatrix}x\\ x_a\end{bmatrix},\quad \dot x_a = r - Cx,\quad
\tilde A = \begin{bmatrix}A & 0\\ -C & 0\end{bmatrix},\quad \tilde B = \begin{bmatrix}B\\ 0\end{bmatrix}$$

对 (Ã,B̃) 配极点/LQR 得 u = −Kx − kₐxₐ。增广系统可控 ⟺ 原系统 (A,B) 可控且 C∉left-null(A,B)（一般成立）。

## 3. LQR（线性二次型调节器）

**连续**：min J = ∫₀^∞ (xᵀQx + uᵀRu)dt，Q=Qᵀ≥0, R=Rᵀ>0，(A,B) 能控（或能稳）：

$$u = -Kx,\quad K = R^{-1}B^{\mathsf T}P,\quad A^{\mathsf T}P + PA - PBR^{-1}B^{\mathsf T}P + Q = 0\ (\text{CARE})$$

性质（设计合理性论据）：全状态 LQR 有无限增益裕度和 ≥60° 相位裕度；K 使闭环 A−BK 稳定。

**Q/R 取法（Bryson 法则，给出初值的标准依据）**：

$$Q_{ii} = \frac{1}{(x_i^{\max})^2},\qquad R_{jj} = \frac{1}{(u_j^{\max})^2}$$

x_i^max 为该状态允许的最大偏差，u_j^max 为最大允许控制量。调参方向：加大 Q 中某项 → 该状态收敛更快但控制量增大；加大 R → 更省力、更平缓。过程：Bryson 初值 → 仿真 → 按指标对角微调 2~3 轮。

**离散（DLQR）**：min Σ xᵀQx + uᵀRu（每步），

$$u[k] = -K_d x[k],\quad K_d = (R + B_d^{\mathsf T}P B_d)^{-1}B_d^{\mathsf T}P A_d$$
$$P = A_d^{\mathsf T}P A_d - A_d^{\mathsf T}P B_d (R+B_d^{\mathsf T}P B_d)^{-1}B_d^{\mathsf T}P A_d + Q\ (\text{DARE})$$

连续设计后离散实现，小 Ts 下 K≈K_d；要求严格最优直接用 DLQR 于 (Ad,Bd)。

**LQR 跟踪**：同 §2 的积分增广，对 (Ã,B̃) 解 LQR（对 x_a 也给权重 Q_a，q_a 小→慢消差、q_a 大→快但易超调）。

**数值**：`linhelper.py lqr()/dlqr()`（scipy 可用时用 solved_*_are，否则 Kleinman/值迭代）。

## 4. 观测器与 Kalman 滤波

### 4.1 Luenberger 观测器（确定性）

$$\dot{\hat x} = A\hat x + Bu + L(y - C\hat x),\qquad \dot e = (A - LC)e$$

极点配置在对偶系统：place(Aᵀ, Cᵀ, p_obs) 得 L = Pᵀ。观测极点取控制极点 3~5 倍远（更快收敛）；传感器噪声明显时放缓至 2~3 倍。**分离原理**：反馈用 x̂ 代替 x，闭环极点 = 控制极点 ∪ 观测极点，二者独立设计。

### 4.2 离散 Kalman 滤波（随机，嵌入式首选）

模型 x[k+1] = Ad x[k] + Bd u[k] + w，y[k] = Cd x[k] + v，w~N(0,Qₙ), v~N(0,Rₙ)：

```
预测:  x̂⁻[k] = Ad·x̂[k−1] + Bd·u[k−1]
       P⁻    = Ad·P·Adᵀ + Qn
增益:  Kf    = P⁻·Cdᵀ·(Cd·P⁻·Cdᵀ + Rn)⁻¹
更新:  x̂[k]  = x̂⁻ + Kf·(y[k] − Cd·x̂⁻)
       P     = (I − Kf·Cd)·P⁻        # 数值敏感时用 Joseph 形式:
                                      # P = (I−Kf·Cd)P⁻(I−Kf·Cd)ᵀ + Kf·Rn·Kfᵀ
```

- Qₙ/Rₙ 调法：Rₙ 取传感器方差（由 datasheet 或静止数据算）；Qₙ 从小往大调——大→跟踪快但噪、小→平滑但滞后。
- 初值 x̂=测量值, P=I×10⁴；协方差 P 离线验证收敛（跑 1000 拍看 P 是否趋稳）。
- 稳态 Kalman 增益可离线算好存成常数（嵌入式推荐）：迭代 DARE 对偶式直至 Kf 收敛。

### 4.3 LQG

LQR（§3）+ Kalman（§4.2）组合：u = −Kx̂。注意：分离设计后**必须回 Stage 8 仿真**验证，LQG 不继承 LQR 裕度。

## 5. MPC（模型预测控制）

### 5.1 标准构造（速度形式，自带积分作用）

用 ZOH 离散模型 (Ad,Bd)，增广状态 ξ = [x; u]，控制量改为增量 Δu：

$$\xi[k+1] = \underbrace{\begin{bmatrix}A_d & B_d\\ 0 & I\end{bmatrix}}_{\tilde A}\xi[k] + \underbrace{\begin{bmatrix}B_d\\ I\end{bmatrix}}_{\tilde B}\Delta u[k]$$

预测 N 步（堆叠）：Y = Φξ[k] + ΓΔU，

$$\Phi = \begin{bmatrix}\tilde C\tilde A\\ \vdots\\ \tilde C\tilde A^{N}\end{bmatrix},\quad
\Gamma = \begin{bmatrix}\tilde C\tilde B & & \\ \vdots & \ddots & \\ \tilde C\tilde A^{N-1}\tilde B & \cdots & \tilde C\tilde B\end{bmatrix}$$

### 5.2 代价与求解

$$J = \sum_{i=1}^{N} (y_{k+i}-r_{k+i})^{\mathsf T}\bar Q (y_{k+i}-r_{k+i}) + \sum_{i=0}^{N_c-1}\Delta u_{k+i}^{\mathsf T}\bar R\,\Delta u_{k+i} + x_{k+N}^{\mathsf T}P_t x_{k+N}$$

展开为 QP：J = ΔUᵀHΔU + 2fᵀΔU + const，H = ΓᵀQ̄Γ + R̄ +（数值加 εI 保正定），f = ΓᵀQ̄(Φξ − R_seq)。

- **无约束解析解** ΔU* = −H⁻¹f，只取第一个分量 Δu[k] 施加，滚动进行（低维系统嵌入式可用）。
- **有约束** min ΔU s.t. 线性不等式（|u|≤u_max、|Δu|≤Δu_max、|x|≤x_max、终端集）→ 在线 QP。实现选项：①解析解 + 饱和投影（工程近似，必须仿真验证不失稳）；②小型活跃集/ADMM 求解器（自写 ≤200 行或移植 OSQP-light）；③上位机算轨迹下发。
- 稳定性保障：终端代价 P_t 取对应 DLQR 的 DARE 解，或加终端约束集。
- 参数经验：Ts ≈ 开环主导时间常数/10~20；N 覆盖 ~2×主导响应时间；控制时域 Nc = 3~10 并配合块移动(move blocking)降维。
- **实时性红线**：一拍内必须解完（先离线测 QP 平均/最坏耗时再定 Ts），嵌入式上首选"无约束解析解+约束投影"或短 horizon。

## 6. 进阶方法（非线性/特殊场景速查）

| 方法 | 场景 | 一句话做法 |
|---|---|---|
| 增益调度 | 对象随工况大变（飞行全包线） | 分工作点各自线性化设计 LQR/PID，插值增益，切换加限速滤波 |
| 反馈线性化 | 模型准、全状态可测 | 用 u = β(x)⁻¹[−α(x) + v] 把动态变成线性积分链，v 再用 LQR/PID |
| 滑模控制 SMC | 鲁棒抗扰、匹配不确定 | s = (d/dt+λ)e，u = u_eq − k·sat(s/φ)，边界层抗抖振 |
| 输入整形 | 吊车/柔性结构残余振荡 | 用整形卷积序列(如 ZV: [1,1]/2 间隔 T_d/2)滤波参考 |
| 自抗扰 ADRC | 模型不准、扰动大 | 用 LESO 估计总扰动并补偿，带宽法整定 |

以上每种都遵循同一报告模板；选了进阶方法时，Stage 4 之后仍需完成离散化与仿真验证门槛。

## 7. 设计完成自检清单（Stage 5 出口检查）

- [ ] 选型理由写明（对照 §0 矩阵）
- [ ] 所有增益有数值（代入参数算过，不是符号）
- [ ] 闭环特征值全部左半平面（打印列表）
- [ ] 跟踪方案含积分增广或前馈（若任务含跟踪）
- [ ] 执行器饱和下的方案（抗饱和结构/控制量峰值预估）
- [ ] mermaid 闭环框图已画
