# CST／CRP 追踪与审核

从 [分步处理 Notebook](../notebooks/Tractography.ipynb) 开始。初学时先运行 [合成数据教程](../notebooks/Tractography_Tutorial.ipynb)，练习审核流线起点、走向、门顺序和覆盖区域。两份 Notebook 都调用同一个 Python 后端。报告由后端生成。

本分支接收已经完成质控的校正 DWI 和批准的结构配准链。已验证保存的父流线筛选、导出、坐标链合成检查及 HTML。原生 MRtrix 拟合与追踪在命令入口中显式开启，新增入口仍待实际 MRtrix 环境验收。新个体使用前检查 DWI 梯度、畸变、脑干覆盖、b0 到个体 T1 的对齐和图谱位置。

## 安装与快速试用

沿用仓库 Python 3.12 环境与锁定分析依赖。另为该环境安装 JupyterLab 和 ipywidgets，然后使用同一个 Python 环境开启 Notebook。真实图谱映射需要 README 里列出的 SimpleITK。只有运行原生拟合和追踪时，才需要本地 PATH 上的 MRtrix 3.0.4。

```powershell
python run_tractography.py demo --output work/tractography_demo
python run_tractography.py report --results work/tractography_demo/screened --output work/tractography_demo/review.html
python -m pip install jupyterlab ipywidgets
python -m jupyterlab
```

演示输入由代码生成。输出目录每次取新名字。正式数据的各阶段输出放在该访视的 results/DWI_CST。正式配置与文件使用本地路径和 SHA-256，配置里的相对路径从配置文件自身位置解析。发布仓库保留演示代码与通用配置模板。

## 分阶段操作

1. **核对 DWI、结构空间和图谱。** 以校正 b0 为输出网格，使用已审的 DWI→T1 刚体坐标变换和同访视 T1↔MNI。检查脑干采集范围、脑域边界及图谱侧别。
2. **准备种子、追踪支持区、门和皮层带。** 种子决定从哪里启动，追踪支持区限制模型可行走的区域，普通解剖门与冻结共同脑域相交。不同模型可各自声明追踪支持区。记录皮层带来源和皮层距离阈值。
3. **按需拟合和追踪。** 用 fit-plan、track-plan 审核配置和命令。fit 生成张量与 FA、MD、AD、RD，并从单壳估计单组织白质 FOD。track 根据配置生成父 TCK 和播种日志。每个模型、起点及侧别的输出都能单独追溯。
4. **筛选。** 用 screen 处理实际父 TCK。它核对文件、记录的流线数及种子坐标，逐条识别皮层起点、向下走向、CST 和 CRP 候选。已有父 TCK 时直接进入此阶段。
5. **人工审核和导出。** 用 report 生成含计数和内嵌 b0 覆盖图的 HTML，用 Notebook 对照三维流线和原生 b0 切片上的覆盖图。用 export 导出单个门的 NIfTI、候选索引及非空的候选 TCK。按同一父集和同一门比较结果。

```powershell
python run_tractography.py fit-plan --config fit.local.json --output work/fit
python run_tractography.py track-plan --config track.local.json --output work/tracking
python run_tractography.py screen --config replay.local.json --output work/screened
python run_tractography.py report --results work/screened --output work/CST_CRP_report.html
python run_tractography.py export --results work/screened --run ST_iFOD2_L_M1 --gate CST_PLIC_plus2_CP_Tang_pons --output work/selected
```

输入字段见 [筛选配置](../configs/examples/tractography_replay.example.json)、[拟合配置](../configs/examples/tractography_fit.example.json)和[图谱映射配置](../configs/examples/tractography_label.example.json)。模板内的 SHA-256 占位符用相应文件的实际哈希替换。完整演示配置可由 demo 生成。

**拟合和追踪的额外字段。** fit 需要在相同 b0 网格准备两个已校正、带梯度表的单壳 MIF，并分别声明目标 b 值。track 使用同一份筛选配置的参考、脑域和 runs 元数据，每个 run 额外登记带哈希的 tracker_input 和 rng_seed。Tensor_Det 的 tracker_input 是校正单壳 DWI，另填 shell_b。ST_iFOD2 的 tracker_input 是 lmax=6 的白质 FOD。父 TCK 和种子日志由追踪生成。原生入口成功后输出 replay.local.json 供 screen 使用。

**原生空间图谱。** map-label 输入已阈值化的二值标签。图谱在 MNI152NLin2009cAsym 时，从目标 DWI 点依次回查个体 T1 和 MNI2009c。图谱在 MNI152NLin6Asym 时，再接锁定的 MNI2009c → Gen6 点变换。配置里的 H5 按影像方向命名。刚体矩阵以 RAS 给出，执行时转换为 LPS。组合后进行一次最近邻采样。图谱和标准空间变换由使用者按各自许可准备。正式图谱映射及真实 `screen` 的侧别判定都读取绑定实际影像、变换和解剖地标的审核记录。用 `audit-label` 在新目录生成待审核候选，见[空间审核回执](SPATIAL_INPUT_REVIEW.md)。

## 这些参数回答什么问题

- **播种与传播。** 默认每组 200,000 次播种尝试。起步和停止阈值均为 0.10，长度范围 5–250 mm。记录 RNG 和线程。不同调度下多线程追踪可能生成不同父流线。查看实际日志中的尝试、生成和保存数。
- **候选门。** 皮层起点距种子掩膜最多 2 mm，皮层带容差由输入声明。净向下位移至少 10 mm，检查沿流线的门顺序。Tang 脑桥和延髓门看同侧候选，Navigator／Boyne 网状结构门依次检查 mRt 到各下游结构。严格 CST 路径用同侧 PLIC、脑脚、原始 Tang 脑桥门。ST-CSD 可对比仅将 PLIC 外扩 2 mm 的候选。M2 浏览入口分别展示 SMA 与 preSMA。
- **覆盖与核心。** 访问图在原生 DWI 网格，每条候选在每个体素最多计一次。按 RAS 物理高度分层，在已访问体素中保留达到该层峰值 15% 的位置。两种图的候选索引一致。核心掩膜用于显示候选分布，候选数本身受到采集、模型、种子和门定义影响。
- **范围。** DICOM 转换、DWI 畸变校正、DWI→T1 注册优化、病灶负荷指标和个体临床解释另行验收。以 mRt 为独立起点和降低局部传播阈值的试验另列。

MRtrix 参数、梯度处理和拟合命令可查官方说明。[tckgen](https://mrtrix.readthedocs.io/en/3.0.4/reference/commands/tckgen.html) · [mrinfo](https://mrtrix.readthedocs.io/en/3.0.4/reference/commands/mrinfo.html) · [dwi2response](https://mrtrix.readthedocs.io/en/3.0.4/reference/commands/dwi2response.html) · [dwi2fod](https://mrtrix.readthedocs.io/en/3.0.4/reference/commands/dwi2fod.html)
