# 分析参考：s域 / 时频域 / 能控能观

## 1. 两条建模路线

- **现代（物理→状态空间，默认）**：由 Stage 1–2 直接得 (A,B,C,D)。
- **经典（物理→微分方程→拉氏变换→G(s)）**：由 G(s) 反推状态空间用**能控标准型**（companion）：

  对 $G(s) = \dfrac{b_{n-1}s^{n-1}+\cdots+b_1s+b_0}{s^n+a_{n-1}s^{n-1}+\cdots+a_0}$：

  $$A = \begin{bmatrix}0&1&&\\&\ddots&\ddots&\\&&0&1\\-a_0&-a_1&\cdots&-a_{n-1}\end{bmatrix},\quad B=\begin{bmatrix}0\\\vdots\\0\\1\end{bmatrix},\quad C=\begin{bmatrix}b_0&b_1&\cdots&b_{n-1}\end{bmatrix}$$

  能观标准型为其转置对偶。约当型用于含重极点系统。报告中默认用物理状态空间，经典路线作为对照。

## 2. 传递函数与极零点

$$G(s) = C(sI-A)^{-1}B + D = \frac{C\,\mathrm{adj}(sI-A)\,B + D\det(sI-A)}{\det(sI-A)}$$

- 极点 = det(sI−A)=0 的根 = 特征值（相似变换不变）。
- n 阶 m 零点系统极零点计数：分母 n 次、分子 m 次；m<n 时有 n−m 个无穷远零点。
- 主导极点：实部最靠近虚轴且未被零点对消的一对，决定动态特性；远极点（>5×）可忽略。
- 手算技巧：n≤3 展开 det(sI−A)；n≥4 用 `linhelper.py`（numpy 特征值）。

## 3. 二阶近似指标公式（背下来，处处用）

对主导极点对 $s^2 + 2\zeta\omega_n s + \omega_n^2$（欠阻尼 0<ζ<1，单位阶跃）：

| 指标 | 公式 |
|---|---|
| 超调量 | $M_p = e^{-\pi\zeta/\sqrt{1-\zeta^2}}\times 100\%$ |
| 峰值时间 | $t_p = \pi/(\omega_n\sqrt{1-\zeta^2})$ |
| 调节时间(2%) | $t_s \approx 4/(\zeta\omega_n)$；(5%) ≈ 3/(ζωₙ) |
| 上升时间(0→100%) | $t_r \approx (\pi - \arccos\zeta)/(\omega_n\sqrt{1-\zeta^2})$ |
| 阻尼角 | 极点位置 $s = -\zeta\omega_n \pm j\omega_n\sqrt{1-\zeta^2}$，与负实轴夹角 θ=arccos ζ |

反查表：ζ=0.5→Mp≈16.3%；ζ=0.707→4.3%；ζ=0.8→1.5%；ζ=1→0%。

## 4. 稳定性判据

- **极点判据（一票否决）**：全部极点实部 < 0 稳定；任一 > 0 不稳定；实轴重根 =0 临界。
- **Routh 判据**（无特征值数值时）：建 Routh 表，第一列全正 ⟺ 稳定；第一列变号次数 = 右半平面极点数。
  - 速查：一阶 a₁s+a₀：a₁,a₀>0；二阶 a₂s²+a₁s+a₀：三者全正；三阶 a₃s³+a₂s²+a₁s+a₀：全正且 a₂a₁ > a₃a₀。
- **Lyapunov（非线性/线性通用）**：线性定常系统稳定 ⟺ 任意 Q=Qᵀ>0，李雅普诺夫方程 AᵀP + PA = −Q 有唯一解 P=Pᵀ>0。非线性系统找 V(x) 正定且 V̇ 负定可证局部/全局稳定（非线性验证用，报告 Stage 8 用数值仿真佐证）。

## 5. 状态转移矩阵

$$\Phi(t) = e^{At} = \mathcal{L}^{-1}\{(sI-A)^{-1}\}$$

计算方法（按优先级）：①A 可对角化 → e^{At} = Ve^{Λt}V⁻¹；②Cayley-Hamilton：e^{At} = Σαᵢ(t)Aⁱ（n 项，αᵢ 由特征值代入解线性方程组）；③数值：`linhelper.py` 的 expm（scaling-and-squaring）。

关键恒等式：x(t) = e^{A(t−t₀)}x(t₀) + ∫ₜ₀ᵗ e^{A(t−τ)}Bu(τ)dτ。

## 6. 能控性与能观性（Stage 4 核心）

**能控性**：能控性矩阵
$$\mathcal{C} = \begin{bmatrix}B & AB & A^2B & \cdots & A^{n-1}B\end{bmatrix}$$
完全能控 ⟺ rank 𝒞 = n。SISO 时 det 𝒞 ≠ 0。

**能观性**：能观性矩阵
$$\mathcal{O} = \begin{bmatrix}C\\CA\\ \vdots \\ CA^{n-1}\end{bmatrix}$$
完全能观 ⟺ rank 𝒪 = n。

**PBH 判据**（找出坏模态）：对每个特征值 λᵢ，
$$\mathrm{rank}\begin{bmatrix}\lambda_i I - A & B\end{bmatrix} = n$$
成立则 λᵢ 可控；不成立的 λᵢ 即不可控模态。对偶地 rank[(λI−A); C] = n 检验可观。报告里必须用 PBH 指名"哪个模态不可控/不可观"。

**对偶性**：系统 (A,B,C) 的能控性 ⟺ 对偶系统 (Aᵀ,Cᵀ,Bᵀ) 的能观性。这就是观测器设计转成对偶系统配极点的根据。

**能稳/能检（部分可控时的判据）**：
- 能稳（stabilizable）：不可控模态全部稳定 → 仍可用状态反馈镇定，但极点不能任意配；
- 能检（detectable）：不可观模态全部稳定 → 仍能构造收敛观测器。
- 结论句模板："系统不完全能控，但不可控模态 {λ} 位于左半平面，系统可稳，状态反馈可实现镇定；极点配置仅能在可控子空间内进行。"

**Kalman 标准分解**（要写过程）：用能控性矩阵 𝒞 做 QR，取其前 rank 列构成 T₁，补正交基 T₂，变换 x̄ = Tᵀx 得分块：
$$\bar A = T^{-1}AT = \begin{bmatrix}A_c & A_{12}\\ 0 & A_{\bar c}\end{bmatrix},\quad \bar B = T^{-1}B = \begin{bmatrix}B_c\\ 0\end{bmatrix}$$
可控子空间 (A_c,B_c) 维数 = rank 𝒞；A_̄c 特征值即不可控模态。能观分解对偶同理。串接做完全分解（可控∩可观 / 可控不可观 / …四块）仅当需要时做。

**结论 → 设计映射**：

| 检验结果 | 极点配置 | LQR | 观测器/Kalman |
|---|---|---|---|
| 完全能控 | 任意配置 | 有唯一正定解 | — |
| 能稳 | 只能镇定可控部分 | P 存在但非唯一镇定意义 | — |
| 完全能观 | — | — | 观测器极点任意配 |
| 能检 | — | — | 可构造收敛观测器 |
| 不能稳 | **失败，物理不可行** | — | — |

## 7. 频域分析（经典补充）

- 频率响应 G(jω) = C(jωI−A)⁻¹B + D；幅值 |G(jω)|、相位 ∠G(jω)。数值扫描：ω 取对数栅格 10⁻²~10⁴ rad/s（每十倍程 50 点）。
- **增益裕度** GM = 1/|G(jω₁)|，ω₁ 为相位 = −180° 之频率；**相位裕度** PM = 180° + ∠G(jω_c)，ω_c 为 |G|=1 之穿越频率。工程要求：PM ≥ 45°（希望 45–60°），GM ≥ 6dB。
- Bode 渐近线：一阶环节 ±20dB/dec、转折处 −3dB、−45°；二阶振荡环节在 ζ<0.5 时谐振峰 Mᵣ ≈ 1/(2ζ√(1−ζ²))。
- 闭环带宽 ≈ ω_c（单位反馈下）；带宽与调节时间关系 ωₙ ≈ 4/(ζ·t_s)。
- 一票否决式检查：开环稳定 + 负反馈下 Nyquist 不包 −1 点 ⟺ 闭环稳定（数值上用闭环特征多项式检验代替，更可靠）。

## 8. 数值实践（配 scripts/linhelper.py）

- 秩判定永远用 `matrix_rank`（基于 SVD，默认容差），**禁止**用 det ≠ 0 判断高阶矩阵。
- SVD 奇异值在 1e-8~1e-10 量级时报告"数值上接近奇异，结论附条件数"。
- 特征值复数用 numpy 直接复数输出，报告里写 a±bj。
- 调用示例：

```bash
python3 scripts/linhelper.py demo        # 自带倒立摆演示
```

```python
from linhelper import ctrb, obsv, rank_report, lqr, c2d
r = rank_report(A, B, C)   # 一次性给出能控/能观/PBH坏模态/极点
K, P, eig_cl = lqr(A, B, Q, R)
Ad, Bd = c2d(A, B, Ts)     # ZOH 精确离散
```
