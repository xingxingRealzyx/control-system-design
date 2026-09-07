# 代码生成参考：C（默认）/ 其他语言 / 定点化 / 对拍测试

## 1. 通用硬性原则（任何语言都适用）

1. **确定性**：无动态内存、无递归、无隐藏全局状态；一切状态在控制器 struct/对象里。
2. **参数与状态分离**：增益、限幅、Ts 在配置结构体，运行期状态在实例结构体；支持多实例。
3. **接口形状统一**：`init(cfg)` 一次 + `update(measurement, setpoint) -> output` 每拍一次。
4. **防护**：所有除法先查零；饱和函数统一 `limsat(v, lo, hi)`；禁止 NaN 传播（异常输入打错误标志并输出安全值，如上次输出或 0）。
5. **单位注释**：每个物理量标注单位（rad、rad/s、N·m、s…），每个增益标注由哪一节设计公式得出。
6. **数值类型显式**：`config.h` 提供 `CTRL_REAL`（float/double 切换），代码一律用 `CTRL_REAL` 不裸写 float。
7. 每拍计算量固定（无数据依赖分支长度差异），保证实时性可测。

## 2. C 代码组织（默认输出）

```
ctrl_<object>.h      # 类型定义 + 接口声明（参数结构体、状态结构体、init/update）
ctrl_<object>.c      # 实现
ctrl_config.h        # CTRL_REAL、限幅、采样率等编译期配置
test_ctrl.c          # 对拍 harness（§6）
```

接口模板：

```c
/* 控制器参数（设计期常量） */
typedef struct {
    CTRL_REAL K[4];        /* 状态反馈增益, Stage5 式(x) */
    CTRL_REAL L[4];        /* 观测器增益, Stage6 */
    CTRL_REAL Ts;          /* 采样周期 [s] */
    CTRL_REAL u_min, u_max; /* 执行器限幅 [V] */
} ctrl_cfg_t;

/* 控制器实例（运行期状态） */
typedef struct {
    const ctrl_cfg_t *cfg;
    CTRL_REAL x_hat[4];    /* 观测状态 [rad, rad/s, ...] */
    CTRL_REAL u_prev;      /* 上一拍输出 */
    int   err;             /* 0=正常, 1=输入异常 */
} ctrl_t;

void ctrl_init(ctrl_t *c, const ctrl_cfg_t *cfg);
CTRL_REAL ctrl_update(ctrl_t *c, CTRL_REAL y, CTRL_REAL r);
```

## 3. PID 参考 C 实现（抗饱和 + 微分滤波，含头文件可直接用）

```c
#include "ctrl_pid.h"

/* 结构体见 .h:
typedef struct {
    CTRL_REAL Kp, Ki, Kd, Tf, Ts, u_min, u_max;
} pid_cfg_t;
typedef struct {
    const pid_cfg_t *cfg;
    CTRL_REAL i_term;   /* 积分累积 */
    CTRL_REAL d_state;  /* 微分滤波器状态 */
    CTRL_REAL y_prev;   /* 上一拍测量 */
    int       first;    /* 首拍标志 */
    CTRL_REAL u_prev;
} pid_t; */

static inline CTRL_REAL limsat(CTRL_REAL v, CTRL_REAL lo, CTRL_REAL hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

void pid_init(pid_t *p, const pid_cfg_t *cfg)
{
    p->cfg = cfg;
    p->i_term = (CTRL_REAL)0;
    p->d_state = (CTRL_REAL)0;
    p->first = 1;
    p->u_prev = (CTRL_REAL)0;
}

CTRL_REAL pid_update(pid_t *p, CTRL_REAL y, CTRL_REAL r)
{
    const pid_cfg_t *c = p->cfg;
    CTRL_REAL e = r - y;

    /* P */
    CTRL_REAL u_p = c->Kp * e;

    /* I + clamping 抗饱和：饱和且误差同向时不积分 */
    CTRL_REAL u_unsat = p->u_prev + c->Kp * e;   /* 预估不含积分项的走向 */
    if (!((u_unsat > c->u_max && e > (CTRL_REAL)0) ||
          (u_unsat < c->u_min && e < (CTRL_REAL)0)))
        p->i_term += c->Ki * c->Ts * e;
    p->i_term = limsat(p->i_term, c->u_min, c->u_max);

    /* D：作用于测量 -y（避免设定值冲击），一阶滤波 tau = Kd/N, N=1/Tf·Kd 关系已折算进 Tf */
    CTRL_REAL d_f = p->d_state;
    if (p->first) { p->y_prev = y; p->first = 0; }
    CTRL_REAL alpha = c->Tf / (c->Tf + c->Ts);
    d_f = alpha * p->d_state
        + ((CTRL_REAL)1 - alpha) * c->Kd * (-(y - p->y_prev)) / c->Ts;
    p->d_state = d_f;
    p->y_prev  = y;

    CTRL_REAL u = limsat(u_p + p->i_term + d_f, c->u_min, c->u_max);
    p->u_prev = u;
    return u;
}
```

（以上为骨架参考；生成时按具体对象重命名、补单位注释与头文件，并保持与设计公式一致。）

## 4. 状态反馈 + 观测器 C 实现要点

```c
CTRL_REAL ctrl_update(ctrl_t *c, const CTRL_REAL y_meas[], CTRL_REAL r)
{
    /* 1. 观测器更新: x̂⁻ = Ad x̂ + Bd u_prev; x̂ = x̂⁻ + L(y − C x̂⁻) */
    /* 2. 积分增广: xa += Ts*(r − C x̂) */
    /* 3. u = −K x̂ − ka xa + Nr r;  u = limsat(u, cfg->u_min, cfg->u_max) */
    /* 矩阵乘法全部展开写常维 for 循环（n 为编译期常量），不做通用矩阵库 */
}
```

要点：矩阵元素用 `static const CTRL_REAL Ad[4][4] = {...}` 写数值，注释标"由 Stage7 ZOH, T=0.01s"；u_prev 存实例；Kalman 稳态增益离线算好后与普通观测器同构。

## 5. 其他语言要点

| 语言 | 要点 |
|---|---|
| C++ | class 封装同结构；禁异常/RTTI（嵌入式）；数组用 std::array 固定维 |
| Python | 交付两份：仿真验证版（numpy 向量化）+ 部署版（纯标量 class，与 C 逐行对应便于对拍）；禁全局可变状态 |
| Rust | struct + impl，状态字段不 pub；用 f32/f64 显式；无 unsafe；饱和用 clamp |
| JavaScript/TS | Float64Array 存状态；类封装；注意没有 int 截断问题但有精度——默认 double |

生成非 C 语言时同样遵守 §1 原则与 §6 对拍。

## 6. 对拍测试协议（Stage 9 门槛，必做）

1. Python 侧：把仿真中控制器"每拍输入(y, r) → 输出(u)"序列写成 CSV（≥2000 拍，覆盖阶跃+反向+饱和段）。
2. 目标语言 harness：读同一 CSV 逐拍喂给控制器，输出 u_c[k] 与期望 u_py[k] 比对。
3. 容差：double — 相对误差 1e-9（或绝对 1e-12）；float — 相对误差 1e-5。超差说明两侧实现不一致（最常见：运算次序、float/double 混用、抗饱和分支差异），必须修到过。
4. harness 模板（C）：

```c
/* 用法: ./test_ctrl golden.csv ; 每行: y,r,u_expect */
while (fscanf(f, "%lf,%lf,%lf", &y, &r, &u_exp) == 3) {
    double u = ctrl_update(&c, (CTRL_REAL)y, (CTRL_REAL)r);
    if (fabs(u - u_exp) > TOL * (fabs(u_exp) + 1.0)) { fails++; }
}
printf("FAIL %d/%d\n", fails, n) / printf("PASS\n");
```

## 7. 定点化指引（仅在用户要求/硬件无 FPU 时）

- 定标：Qm.n，m=整数位含符号，n=小数位；转换 x_q = round(x·2^n)，范围 ±2^m。
- 采样/增益等设计量全部定标表列出（变量 | Q 格式 | 范围 | 分辨率）。
- 乘法后右移 n 位，注意中间结果位增长（(m1+n)+(m2+n) 位需 32/64 位累加器）。
- 积分器、滤波器状态用更高精度（如 Q15 状态配 Q30 累加）。
- 每处乘加后显式饱和（DSP 的 SAT 位或手写 limsat）。
- 定点版必须单独对拍（容差放宽为量化步长量级）并报告最大误差拍。

## 8. 交付物清单（Stage 9 完成定义）

- [ ] `<lang>` 源文件 + 头文件/接口 + config
- [ ] 使用示例（10 行内：init → 周期调用 update）
- [ ] 对拍 harness + 通过结果（贴出容差与结果行）
- [ ] README 段落：参数含义、如何替换实测参数、已知局限（线性化工作点、假设参数）
