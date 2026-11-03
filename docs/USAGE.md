# 独立工作流使用说明

## 文件与空间

代码目录不存放原始影像。建议外部工作目录按人/访视组织：

```text
/analysis/data/sub-demo/sub-demo_s01/
  raw/                    原始 DICOM，保持只读
  preprocessing_v1/
    bids/                 已审阅的 BIDS 输入副本
    derivatives/fmriprep/  输出
    logs/                 本工具的实际命令、退出状态
  results/task_fMRI/v1/    阶段登记、模型版本和报告
```

配置中 `preprocessing.project` 指向含 bids、derivatives、logs 的目录；`workspace` 指向另一个不相交的结果目录。`output_roots` 显式允许个体结果根；不填写时仅允许 `--root/results` 内输出。不要用此版本直接接管本机旧结果状态，不要把本地迁移位置表或旧收据改写成新实验。

文件在指定 `--root` 下解析；`analysis.runs` 的 BOLD、mask、confounds 路径和 `analysis.anatomy` 特别规定为相对 `project/derivatives/fmriprep`。事件和 TR JSON 则相对 `--root` 或使用绝对路径。示例位于 configs/examples，默认 template_only，复制到外部工作目录后填写，不直接执行。

同一个模型中的 BOLD、结构背景和掩膜应处于声明的同一空间；模型掩膜必须同网格。原生 T1 与 MNI 不能因数组尺寸相同就视为同空间。标签重采样用最近邻；连续图像用适当连续插值。配准质量另看非病灶解剖对应与形变，不能只看往返 Dice。

## 新预处理：Linux Docker 后端

这是独立版新增的移植后端，尚未进行真实容器验收，区别于本地主流程已测过的专用 WSL 调度。

1. 由使用者启动自己的 Linux/WSL、Docker；本工具不启动或修复它们。
2. 自行准备固定 fMRIPrep 25.2.7 镜像、TemplateFlow 缓存、FreeSurfer license。运行 `docker image inspect nipreps/fmriprep:25.2.7 --format '{{.Id}}'` 获取实际本地镜像 ID，填入阶段配置 `expected_image_id`。不要用示例全零值代替核验。
3. 将原始输入转换成审阅后的 BIDS staging。DICOM 转换和医院元数据建账仍由外部流程完成，首版不包含自动 DICOM intake。
4. 在 staging 中准备与 `bids_filter` 完全一致的 `neo_bids_filter.json`；在 `.bidsignore` 中列出该辅助文件。工具不编辑原始 BIDS 或原始 DICOM。
5. 填写外部 `*.local.json` 和 runtime。work 是每个新版本专属的空目录，与输入/输出项目分离，建议位于 Linux 原生 ext4；输出也宜在快速原生盘完成后自行校验导出。不要把 Windows 的路径字符串直接塞给 Linux Docker，也不要连接远端 Docker daemon。
6. 先执行 `neo-mri container-plan ...` 检查命令；它只读配置，不探测 Docker。人工核对后启用 runtime 与 scientific_protocol 两个 execution_enabled，再明确执行 preprocess。

固定工具线程4，工作流预算与容器配额分别配置；不使用 sloppy、low-mem，不自动更换输出空间，不自动忽略 BIDS 验证。容器名 neo-mri-exclusive 避免同名重叠启动；冲突时保留失败，人工协调。模板 cache 目前仅声明目录，没有提供跨机器完整内容锁，这是首版复现限制。

失败退出、服务不可达、缺失输入和不匹配的科学条件均停止；容器和日志保留。此预览不提供断点缓存恢复控制器。终端中断也不意味着容器已停止，须查看日志保存的容器名并自行核实，不能盲目再启动。

只做结构像时 runtime 设置 anat_only=true，任务时序设 paused，analysis.runs=[]；登记后可 report 和查看结构影像，不能运行有意义的任务 GLM。

## 已有预处理与任务分析

`mode=existing` 需要 project/logs/preprocessing_status.json 中真实 exit_code=0，以及 fmriprep_run_configuration.json 中真实 arguments（用于读取原执行线程）。这两份必须来源于真实运行日志，不应为导入任意结果捏造成功记录。不满足时先讨论专门的第三方导入器，不能绕过证据要求。

执行 preprocess 后生成 SHA256 登记。后续快检大小与修改时间，变化时核实哈希；加 --verify-full 可强制完整核验。不能篡改旧收据来接受变化。统计参数变化产生独立模型目录，未完成的同请求结果保留并停止，须人工检查。

GLM 默认包含24项头动、实际 combined aCompCor、已有 cosine 与异常帧列；任务 HRF 和导数由参数决定。该配置不适合所有数据：常量 confound、秩不足、时序超出扫描长度等问题会停止，不能偷偷删列或截断事件。审计的名义残差自由度不是有效独立样本量，也不是临床质量分。

任务时序 confirmed 需要与预处理指纹、事件 SHA 绑定的 clearance；示例已提供所需字段。共享事件探索需明确授权并保留逐次日志未验证的状态。本预览仍保留历史特定探索门，未将其泛化成通用任务授权系统。纯合成输入只能标为 engineering_subset。

报告展示保存的统计图、设计矩阵和解释；不自动重拟合。多 run 固定效应需要任务和条件相容，不用于群体推断；FDR 只针对本模型共同脑掩膜。静息态可预处理及查看，但首版不进行任务外推或自动静息连接分析。

## CFM、SDC 和配准检查

协议规划 `neo-mri plan --config protocol.json` 不读影像或探测容器。CFM 要求选定 T1 网格的二值病灶 ROI；测量 SDC 要求真实场图或反向编码 EPI、相位编码/读出时间与明确 B0 关联。缺失时停止，不自动换 SyN。协议分支拟合前还需要与实际输出指纹及证据文件哈希绑定的 scientific_protocol_review.json；不得把请求的方法当成已经发生的方法。

registration_audit、registration_geometry 消费已经产生的掩膜、对应点与物理形变，不运行配准，也不自动选赢家。T2/FLAIR 刚性配准、病灶传播和多时相专属交互页尚未移植成通用命令。CST、DoubleScanner 与 Mv-SPFM 后续按独立模块接入。
