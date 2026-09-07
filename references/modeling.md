# 建模参考：运动学 / 动力学 / 线性化

## 1. 建模路线选择

| 路线 | 适用 | 特点 |
|---|---|---|
| 拉格朗日 | 多自由度机械系统（摆、机械臂、小车系统） | 纯能量法，不画受力图，约束力自动消去，**默认首选** |
| 牛顿-欧拉 | 刚体链、已知各部件受力关系 | 需画自由体图，含内力，适合简单 2 体系统 |
| D-H 参数 | 开链机械臂运动学/动力学 | 系统化建 4×4 变换矩阵链 |
| 能量法(功率平衡) | 电机、传动系统 | 简单回路系统最快 |

## 2. 拉格朗日五步法（机械系统标准流程）

1. **选广义坐标** q = [q₁…qₙ]：取独立自由度，角度取"静止方向到杆的夹角"，明确正方向并画出坐标图（mermaid 或文字描述）。
2. **写各质点位置与速度**：用 q, q̇ 表示每个刚体质心的位置矢量，求出速率平方 vᵢ²。
3. **写能量**：动能 T = Σ(½mᵢvᵢ² + ½Iᵢωᵢ²)，势能 V = Σmᵢg·hᵢ（+ 弹簧 ½kx²）。
4. **拉格朗日方程**：L = T − V，对每个 qᵢ：
   $$\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} - \frac{\partial L}{\partial q_i} = Q_i$$
   Qᵢ 为广义力（输入力/力矩 + 非保守力）。
5. **整理成标准形** M(q)q̈ + C(q,q̇)q̇ + G(q) = Bu（+ 摩擦项 D·q̇ 等）。

注意：θ̇² 项是向心力项（cos θ 乘积项保留，sin/cos 的平方项展开时勿丢 sinθ·θ̇²）。

## 3. 运动学速查

- 平面旋转：点 p 绕原点转 θ → $p' = \begin{bmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{bmatrix}p$。
- 刚体上任一点速度：v_P = v_O + ω × r_{OP}。
- 欧拉角 ZYX（yaw-pitch-roll，航空/四旋翼惯例）：角速度映射
  $$\begin{bmatrix}p\\q\\r\end{bmatrix} = \begin{bmatrix}1&0&-\sin\theta\\0&\cos\phi&\sin\phi\cos\theta\\0&-\sin\phi&\cos\phi\cos\theta\end{bmatrix}\begin{bmatrix}\dot\phi\\\dot\theta\\\dot\psi\end{bmatrix}$$
  小角度线性化时直接 p≈φ̇, q≈θ̇, r≈ψ̇。
- 机械臂正运动学：ᵢ₋₁Tᵢ = Rot(z,θᵢ)Trans(z,dᵢ)Trans(x,aᵢ)Rot(x,αᵢ)，逐项相乘；逆运动学一般用数值迭代（牛顿法 on 位置误差），控制中常用笛卡尔空间 PID 于关节空间映射。

## 4. 工作点线性化（Stage 2 标准操作）

非线性系统 ẋ = f(x, u)，平衡点 (x₀, u₀) 满足 f(x₀,u₀) = 0。

1. 解平衡点方程（倒立摆类：顶点 θ₀=0 配合求支撑输入 u₀；悬垂摆 θ₀=π 或直接以偏差角 θ=π−φ 为变量取 θ₀=0）。
2. 雅可比：
   $$A = \left.\frac{\partial f}{\partial x}\right|_{(x_0,u_0)},\qquad B = \left.\frac{\partial f}{\partial u}\right|_{(x_0,u_0)}$$
   手工求导：对每个状态方程分别对每个状态/输入求偏导，逐元素写矩阵；无法手工时用数值差分（中心差分步长 1e-5·量程）。
3. 输出方程 y = Cx + Du：C 由"用户能测哪些物理量"决定，D 一般为 0。
4. 偏差变量 δx = x − x₀，δu = u − u₀；模型写 δẋ = Aδx + Bδu 后统一去掉 δ（报告中说明"此后 x 均指偏差量"）。

## 5. 常见非线性项处理

| 项 | 线性模型 | 非线性仿真 |
|---|---|---|
| 库仑摩擦 F = μN·sign(v) | 丢弃 | 保留（注意 sign(0) 取 0） |
| 粘性摩擦 bv | 并入 A | 保留 |
| 执行器饱和 | 不建模 | 必须 limsat(u, u_max) |
| 量化/死区 | 丢弃 | 可选保留 |

## 6. 算例 A：弹簧-质量-阻尼（SISO 二阶）

m ÿ + c ẏ + k y = F。

状态 x = [y, ẏ]ᵀ，u = F：

$$A = \begin{bmatrix}0 & 1\\ -k/m & -c/m\end{bmatrix},\quad B = \begin{bmatrix}0\\ 1/m\end{bmatrix},\quad C = \begin{bmatrix}1 & 0\end{bmatrix},\quad D=0$$

G(s) = 1/(ms² + cs + k)，ωₙ = √(k/m)，ζ = c/(2√(km))。

## 7. 算例 B：电枢控制直流电机（SISO，机电耦合）

电气方程 L di/dt + R i = V − Kₑω；机械方程 J ω̇ + b ω = Kₜ i。状态 x = [ω, i]ᵀ：

$$A = \begin{bmatrix}-b/J & K_t/J\\ -K_e/L & -R/L\end{bmatrix},\quad B = \begin{bmatrix}0\\ 1/L\end{bmatrix},\quad C = \begin{bmatrix}1 & 0\end{bmatrix}$$

$$\frac{\Omega(s)}{V(s)} = \frac{K_t}{(Js+b)(Ls+R) + K_tK_e}$$

速度环经典结论：反电动势项 Kₑ 提供"天然速度反馈"；位置控制时状态取 [θ, ω, i]。

## 8. 算例 C：小车倒立摆（4 阶欠驱动，拉格朗日全过程示范）

小车质量 M，摆杆质量 m，转轴到质心长 l，绕质心转动惯量 I = ml²/3（均质杆），输入水平力 F，θ 为摆杆与**竖直向上**的偏角（θ=0 为倒立点，正方向与 x 正向同侧倾倒）。

**质心位置**：x_p = x + l sinθ，y_p = l cosθ；v_p² = ẋ² + l²θ̇² + 2lẋθ̇cosθ。

**能量**：T = ½Mẋ² + ½m v_p² + ½Iθ̇²，V = mgl cosθ。

**代入拉格朗日方程**整理得非线性方程：

$$(M+m)\ddot x + ml\ddot\theta\cos\theta - ml\dot\theta^2\sin\theta = F$$
$$(I+ml^2)\ddot\theta + ml\ddot x\cos\theta - mgl\sin\theta = 0$$

**在 θ=0 线性化**（cosθ≈1, sinθ≈θ, θ̇²≈0），令 D = (M+m)(I+ml²) − m²l²：

$$\begin{bmatrix}\ddot x\\ \ddot\theta\end{bmatrix} = \frac{1}{D}\begin{bmatrix}(I+ml^2)F - m^2gl^2\,\theta\\ -ml\,F + (M+m)mgl\,\theta\end{bmatrix}$$

取 x = [x, ẋ, θ, θ̇]ᵀ：

$$A = \begin{bmatrix}0&1&0&0\\ 0&0&-\dfrac{m^2gl^2}{D}&\dfrac{I+ml^2}{D}\\ 0&0&0&1\\ 0&0&\dfrac{(M+m)mgl}{D}&0\end{bmatrix},\quad B = \begin{bmatrix}0\\ \dfrac{I+ml^2}{D}\\ 0\\ -\dfrac{ml}{D}\end{bmatrix}$$

典型取值 M=1.0kg, m=0.1kg, l=0.5m, g=9.81：可算得 A 有两个右半平面极点（倒立点不稳定），且**摆角不可控当且仅当 B 中摆通道为零**——本结构下完全能控（用 linhelper 验证）。

## 9. 其他常用对象速查

- **四旋翼姿态**（滚转通道简化）：Iₓ φ̈ = ℓ·u₂ − ... 小角度下双积分 + 电机一阶滞后 G(s) ≈ K/(s²(T_m s+1))；位置控制用串级：位置环 → 期望加速度 → 期望姿态角 → 姿态环。
- **球杆系统**：小球加速度 ≈ (5/7)g sinα（实心球），α 由电机转角决定，二阶欠驱动同倒立摆套路。
- **一级/二级倒立摆、龙门吊（吊车防摆）**：均为拉格朗日多坐标推广；吊车注意摆角在加速结束后是**不能控**欠驱动模态需轨迹整形或输入整形。
- **船舶/车辆航向**：Nomoto 一阶 G(s) = K/(s(Ts+1))。

## 10. 参数表规范

报告 Stage 1 末尾必须给：

| 符号 | 含义 | 数值 | 单位 | 来源 |
|---|---|---|---|---|
| M | 小车质量 | 1.0 | kg | 假设 |
| … | … | … | … | 给定/假设 |

"假设"项在 Stage 10 假设回顾中再次列出，提醒用户替换实测值后需重跑 Stage 3–8。
