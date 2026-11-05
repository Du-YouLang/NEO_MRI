# NEO_MRI · T1/T2结构像发布候选

版本：**0.2.0rc1**。本轮发布只包含已经锚定的T1/T2结构流程与审核工具；DWI/CST、BOLD、GLM、SDC、DoubleScanner和Mv-SPFM不在本安装包内。没有Notebook、患者数据、模型权重、许可证或个人路径配置。

T1默认HD-BET 2.0.1 CPU/noTTA；SynthStrip和ANTs两个脑提取模板显式可选，不自动后备。原生脑掩膜人工批准后，导出审核T1/mask，再进行组织分割与T1→MNI，保存实际mask采用证据，结合闪烁/切片、往返与形变审核。工具固定4线程。

T2使用已冻结的`T2_structural_v1`。默认HD-BET脑域为主结果，T1脑域为补充/显式可选；原始强度/N4均保留，不擅自指定唯一强度方案。已有批准T1/MNI变换复用，单次组合重采样与两次插值对照保留。详见[T2说明](docs/T2_STABLE_V1.md)。

## 安装

使用Python 3.12，在独立环境中执行：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-lock.txt
.venv\Scripts\python -m pip install '.[anatomy]' --no-deps
.venv\Scripts\python -m pip install 'SimpleITK==2.5.6' --no-deps
.venv\Scripts\python -m pytest
```

核心分析依赖保持原锁，不升级。HD-BET、SynthStrip、ANTs和fMRIPrep各自单独准备，不随pip自动下载模型或启动服务。HD-BET推荐独立CPU环境；fMRIPrep固定25.2.7且必须核对本地不可变镜像ID，许可证自备。

## 命令与处理顺序

```powershell
neo-mri --help
neo-mri-extract policy
neo-mri-extract plan --config extraction.local.json
neo-mri-mask --help
neo-mri t1-plan --config T1.local.json
neo-mri-t2-v1 policy
neo-mri-t2-v1 plan --config T2.local.json
```

plan不运行影像计算。run需要明确的execution_enabled、有效输入与本地运行环境；Linux容器需要自己预先启动并完成资源准入，代码不启动/修复WSL或Docker。

1. 准备明确身份的三维NIfTI/BIDS和N4全头T1。本版接入已准备的N4图；不提供从任意医院DICOM到完整T1输出的一键通用转换/N4入口。
2. 用HD-BET等生成候选，在原生红线页逐层审核；记录明确批准，必要时ITK-SNAP手修并重新审核。
3. 导出与原始BIDS T1身份绑定的审核解剖输入，再计划anat-only fMRIPrep。例子：configs/examples/t1_case.example.json与t1_runtime.example.json；所有路径相对配置文件所在目录。
4. 运行结束检查实际mask采用及变换输入链。通用容器启动器仍是发布候选，未在这一轮重跑真实容器验收；不能把既往专用控制器的九分支结果写成该新入口的整体验收。后端未重复导出原生mask时，需要在固定容器内采集anatomy_runtime_evidence并在T1配置提供runtime_evidence，不能复制输入冒充输出。
5. 检查MNI闪烁/逐层切片、往返一致性及非病灶/形变；批准后才给T2配置使用。跨模态T2配准用已审T1作为个体中心，不重新拟合标准化。

## 验证与边界

本地证据：T1三人九分支mask/标准化对照；T2两例卒中原始/N4与两种脑域对照；另外一例健康FLAIR仅为技术演示。人工批准绑定当时的文件，不传给新个体。脑域Dice、病灶往返Dice都不是临床定位准确率。

见[发布与验收范围](docs/RELEASE_SCOPE.md)、[T1使用说明](docs/T1_MASK_REVIEW.md)、[T2稳定版本](docs/T2_STABLE_V1.md)。Jupyter分步交互封装后续再做。

## 病灶分割接口（研究预览）

另提供 `neo-mri-lesion` 独立命令，支持 3D Slicer 人工勾画、外部 SynthStroke 候选、版本审核与标准 NIfTI 掩膜导出。模型、权重和推理环境由使用者自行准备；不更改现有 T1/T2 入口、冻结引擎或脑掩膜审核。T2/FLAIR 背景复用现有配准结果，不新增另一套配准。

见[病灶接口](docs/LESION.md)、[Slicer 操作指南](docs/SLICER_GUIDE.md)、[SynthStroke 接口](docs/SYNTHSTROKE.md)及[本次验证范围](docs/LESION_VALIDATION.md)。CFM 自动接入及功能影像工具本次不发布。
