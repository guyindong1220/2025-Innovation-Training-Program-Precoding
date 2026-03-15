# [cite_start]频谱效率优化中的动态预编码设计方法研究 

![MATLAB](https://img.shields.io/badge/MATLAB-R2023a+-blue.svg)
![QuaDRiGa](https://img.shields.io/badge/QuaDRiGa-Channel%20Model-brightgreen.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

## 📖 1. 项目简介 (About The Project)

本项目研究频谱效率优化中的动态预编码设计方法，重点分析混合预编码技术在毫米波大规模MIMO系统中的应用 。在传统的全数字预编码结构中，海量天线对应的模数/数模转换器及射频功放会带来不可估量的硬件成本和功耗瓶颈 。为此，联合数字与模拟域的混合预编码架构应运而生，实现了硬件复杂度与系统性能的良好折中 。

然而，在动态移动通信中，终端处于高速移动状态时信道矩阵发生快速演进，不可避免地导致预编码矩阵产生误差（即信道老化现象）。

为实现动态场景下的精确预编码追踪，本项目跳出了传统的欧氏空间线性预测框架，将预编码矩阵的正交属性映射至 Stiefel 流形。通过建立严格的常微分方程初值问题（IVP），提出了一种基于 Stiefel 流形的 $C^1$ 贝塞尔插值与测地线外插预测算法。数值实验表明，该方法在显著降低硬件复杂度的同时，不仅能够逼近全数字预编码的系统容量，更在动态追踪场景下展现出极高的逼近精度与鲁棒性，有效提升了系统的频谱效率。
## 🧠 3. 理论与核心算法 (Mathematical Theory & Algorithms)

本项目的核心算法架构分为静态多用户干扰抑制与动态流形追踪两个阶段：

### 3.1 混合预编码系统模型
在毫米波大规模 MIMO 系统中，为平衡极高的硬件功耗，基站端采用两阶段混合预编码架构。发射信号经过基带数字预编码器 $\mathbf{F}_{\mathrm{BB}}$ 与射频模拟预编码器 $\mathbf{F}_{\mathrm{RF}}$ 的联合处理，其中射频层受限于恒模约束。

### 3.2 WMMSE 多用户优化算法
针对单小区多用户环境下的严重干扰问题，本项目利用加权最小均方误差 (WMMSE) 算法进行纯数字域的理想干扰压制。
* 算法通过巧妙的代数变换，将高度非凸的加权和速率最大化问题，等价转化为矩阵加权均方误差 (MSE) 最小化问题。
* 采用块坐标下降法，通过交替更新预编码矩阵 $\mathbf{F}_i$、接收器矩阵 $\mathbf{W}_i$ 和权重矩阵 $\mathbf{U}_i$，在满足总功率约束的前提下逼近系统极限容量。

### 3.3 Stiefel 流形上的测地线外插 (Geodesic Extrapolation IVP)
在终端连续移动的场景下，为了克服信道老化带来的预编码失效，本算法不再依赖传统的欧氏空间线性外推。由于数字预编码矩阵的列正交性，其解空间并非平坦的欧氏空间，而是构成了一个特定的黎曼流形——Stiefel 流形 $\operatorname{St}(n,p)$。
* **初值问题 (IVP) 建模：** 我们将未来的预编码演进物理等效为流形上的无外力惯性运动，即建立满足 $\nabla_{\dot{\gamma}}\dot{\gamma}(t) = 0$ 的测地线方程。
* **边界条件确定：** 提取 $C^1$ 贝塞尔插值曲线末端的终点位置 $\mathbf{P}_n$ 与末端切向量 $\mathbf{V}_{end}$ 作为常微分方程的初始状态。
* **指数映射求解：** 利用黎曼几何理论，该 IVP 的唯一解可通过流形的指数映射显式求出：$\hat{\mathbf{P}}(t) = \operatorname{Exp}_{\mathbf{P}_n} \left( (t - n) \cdot \mathbf{V}_{end} \right)$。这种保约束的几何外推在极大降低实时计算负担的同时，维持了优异的通信容量。

## ⚙️ 4. 环境与依赖 (Prerequisites)

本项目基于 MATLAB 开发，在运行仿真脚本前，请确保您的计算环境满足以下配置要求：

* **核心计算环境:** MATLAB (推荐 R2023a 或更高版本)
* **信道生成引擎:** [QuaDRiGa 工具包](https://quadriga-channel-model.de/)。本项目使用其生成具有空间一致性的 3GPP 38.901 UMa 真实感信道数据。
* **依赖的 MATLAB 工具箱:**
  * Communications Toolbox
  * Optimization Toolbox
  * Phased Array System Toolbox。
## 📁 5. 目录结构 (Repository Structure)

本项目按照论文的章节逻辑进行严密组织，核心算法均存放在 `code` 目录下。完整的真实目录树如下所示：

```text
├── README.md
├── code/
│   ├── 第二节/
│   │   ├── 数值实验/
│   │   │   ├── OUP.m               # 最优无约束预编码 (SVD/MMSE)
│   │   │   ├── OMP.m               # 正交匹配追踪混合预编码算法
│   │   │   └── BS.m                # 波束导向 (Beam Steering) 预编码
│   │   └── code.m                  # 第二节辅助测试与可视化主程序
│   ├── 第三节/
│   │   └── WMMSE_MMSE.m            # WMMSE 与 MMSE 算法性能及和速率分布对比
│   ├── 第四节/
│   │   ├── Achievable_rates...m    # 单径理想信道下单小区多用户可达速率分析
│   │   └── Coverage_probability...m# 蜂窝网络覆盖概率与性能分析
│   ├── 第五节/
│   │   └── run_quadriga_submit.m   # 动态信道生成、用户移动轨迹及混合预编码动态追踪对比
│   ├── 第六节/
│   │   └── InterpolationOnManifolds.m # 球面上的 C^1 贝塞尔曲线验证实验
│   ├── 第七节/
│   │   ├── LoS_5users_data...m     # 5用户视距 (LoS) 场景流形追踪数据生成
│   │   ├── stiefel_utils.m         # 指数映射、对数映射等核心黎曼几何算子
│   │   ├── verify_stiefel_map...m  # Stiefel 流形映射互逆性数值验证
│   │   └── stiefel_interpolati...m # Stiefel 流形贝塞尔插值验证实验
│   ├── 第八节/
│   │   ├── line_waicha.m           # 直线轨迹下的测地线外推算法与误差分析
│   │   └── sin_waicha.m            # 正弦曲线轨迹下的测地线外推鲁棒性验证
│   └── 第九节/
│       ├── H_data.mlx              # 单小区多用户系统级动态信道数据生成与处理 (实时脚本)
│       ├── my_v.mlx                # 综合预测实验：流形内插重构与外推和速率评估 (实时脚本)
│       └── narrowband_channels.mat # 预生成的窄带信道核心数据文件
└── figures/                        # 仿真结果图表保存目录
```
## 📊 6. 核心复现指南 (Reproducibility & Figure Mapping)

本项目高度重视学术研究的可复现性。各章节代码采用了“主执行脚本 + 底层函数库”的结构。您只需在 MATLAB 中运行各章节对应的主脚本，即可一键生成论文中对应的核心图表与结论：

| 执行主脚本 (Main Script) | 调用的核心函数 / 数据 | 对应论文图表 / 章节 | 生成内容与学术结论描述 |
| :--- | :--- | :--- | :--- |
| `第二节/code.m` | `OUP.m`, `OMP.m`, `BS.m` | 第二节 / 图 4 | 预编码基础基准： 生成 64×16 与 256×64 MIMO 系统下的全数字（OUP）、混合预编码（OMP）与波束导向（BS）性能对比图。 |
| `第三节/WMMSE_MMSE.m` | (内置算子) | 第三节 / 节内散图 | 纯数字域优化： 生成 WMMSE 与 MMSE 算法在单小区多用户场景下的和速率对比曲线及箱线图分布。 |
| `第四节/Achievable_rates...m` <br> `Coverage_probability...m` | (内置算子) | 第四节 / 图 5 | 混合架构理论性能： 验证单径理想信道下的单小区多用户可达速率，及蜂窝网络覆盖概率。 |
| `第五节/run_quadriga_submit.m` | (QuaDRiGa 接口) | 第五节 / 图 6, 7, 8 | 动态追踪基准测试： 生成多用户 LoS 二维移动轨迹，输出相邻快照间的 Chordal 距离漂移；对比 WMMSE 优化解与近似解的动态鲁棒性。 |
| `第六节/InterpolationOnManifolds.m`| (内置几何算子) | 第六节 / 节内散图 | 流形几何验证： 生成单位球面上的多节点（如5节点、10节点）$C^1$ 贝塞尔曲线插值三维可视化图。 |
| `第七节/verify_stiefel_map...m` <br> `stiefel_interpolati...m`| (内置算子) | 第七节 / 节内散图 | Stiefel 流形算子测试： 验证指数/对数映射互逆性（精度可达极高量级），及流形插值重建误差分析。 |
| `第八节/line_waicha.m` <br> `sin_waicha.m` | (内置算子)| 第八节 / 图 9, 10 | 流形外推核心实验： 执行测地线外推预测。对比不同预测步长下 Chordal 误差，及直线与正弦曲线轨迹下的频谱效率损失性能。 |
| `第九节/my_v.mlx` <br> `H_data.mlx` <br> `narrowband_channels.mat`| (内置算子) | 第九节 / 图 11, 12, 13, 14 | 系统级应用验证： 压轴综合仿真。生成偶数节点内插重构与未来节点测地线外推的系统总和速率评估，验证实际动态信道下的预测表现。 |
## 🚀 7. 快速开始 (Quick Start)

**步骤 1：克隆仓库**
首先，将本项目克隆到您的本地计算机上：
```bash
git clone https://github.com/guyindong1220/2025-Innovation-Training-Program-Precoding.git
cd 2025-Innovation-Training-Program-Precoding
```
**步骤 2：配置 MATLAB 路径**
打开 MATLAB，将本项目根目录及所有子文件夹添加到工作路径中，以确保各脚本能顺利读取跨章节保存的数据集（如 `.mat` 文件）：
```matlab
addpath(genpath(pwd));
```
**步骤 3：配置 QuaDRiGa 信道模型 (可选但强烈推荐)**
本项目第五节与第九节的动态信道数据生成高度依赖 QuaDRiGa 工具包。如果您希望完整复现信道生成过程，请提前下载 [QuaDRiGa](https://quadriga-channel-model.de/) 并将其添加至 MATLAB 路径：
```matlab
addpath(genpath('您的_QuaDRiGa_安装路径'));
```
*(注：如果您仅想验证第八节等核心流形预测算法，可直接使用仓库中已预生成的 `.mat` 数据，跳过此步)*

**步骤 4：运行核心实验 (以第八节为例)**
路径配置完成后，您就可以进入对应的章节目录运行主脚本了。例如，一键复现正弦曲线轨迹下的测地线外推实验：
1. 在 MATLAB 当前文件夹中导航至 `code/第八节/` 目录。
2. 打开 `sin_waicha.m` 并点击运行。
3. 稍等片刻，即可直接生成论文中对应的用户轨迹、Chordal 距离误差与频谱效率对比图。

## 🤝 8. 致谢与参考文献 (Acknowledgements & References)

* **核心几何理论:** 本项目的流形优化算子基于 P.-A. Absil 等人提出的矩阵流形优化框架，以及 P. Y. Gousenbourger 等人提出的黎曼流形分段贝塞尔插值理论。
* **混合预编码架构:** 参考了 Ahmed Alkhateeb 等人提出的多用户毫米波系统混合预编码可达速率分析框架。

**核心参考文献:**
1. P-A Absil, Robert Mahony, and Rodolphe Sepulchre. Optimization algorithms on matrix manifolds. Princeton University Press, 2009.
2. Ahmed Alkhateeb, Robert W. Heath Jr., and Geert Leus. Achievable rates of multi-user millimeter wave systems with hybrid precoding, 2015.
3. Pierre Yves Gousenbourger, Chafik Samir, and P. A. Absil. Piecewise-bézier $C^1$ interpolation on riemannian manifolds with application to 2d shape morphing. IEEE Computer Society, 2014.