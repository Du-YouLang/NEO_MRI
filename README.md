# NEO_MRI · 首版研究预览

结构 MRI 与基础任务态 fMRI 的分阶段研究工具。当前版本 `0.1.0.dev1`。
T1 强度校正、脑提取、组织分割和标准化由固定版本 fMRIPrep/ANTs 等软件完成；本项目负责输入约束、显式调用、结果登记、模型审计、GLM、FDR 和报告。

这份独立共享包未包含患者影像、个人配置、FreeSurfer 许可证、历史报告或 EEG 代码。

## 安装与最小测试

使用 Python 3.12，Windows 可运行统计、审计、报告和查看器；新容器预处理须在用户自己启动的本地 Linux/WSL shell 中运行。

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-lock.txt
.venv\Scripts\python -m pip install --no-deps -e .
.venv\Scripts\python -m pytest
.venv\Scripts\neo-mri-demo --output ..\synthetic_demo
```

Linux 将 `.venv\Scripts\python` 改为 `.venv/bin/python`，命令入口使用 `.venv/bin/neo-mri` 等。
`requirements-lock.txt` 保留原分析环境的精确版本快照，未把它冒充所有平台都可安装的完整锁。若某版本在索引不可得，应记录失败并讨论，不自动替换版本。

合成示例包含 6×6×6 体素、100 帧随机 BOLD、模拟 confounds 和事件表。它演示真实 GLM 数值计算与复用，不验证 fMRIPrep、卒中临床定位或真实任务效应。

## 分开运行

相对路径统一以 `--root` 指定的工作目录为基准；不指定时使用当前目录。代码位置和数据位置可以完全分离，不需要完整 NEO 工作区或 MRI_DataQC。

```text
neo-mri status --root /analysis --config /analysis/case.local.json
neo-mri preprocess --root /analysis --config /analysis/case.local.json
neo-mri audit --root /analysis --config /analysis/case.local.json
neo-mri glm --root /analysis --config /analysis/case.local.json
neo-mri report --root /analysis --config /analysis/case.local.json
```

只有明确的 `preprocess` 才可能启动新容器；`mode=existing` 只登记具有完成证据的既有结果。只改统计参数执行 `glm`、`report`；只改显示执行 `report`。缓存按输入、时序、引擎和参数指纹分版本；更改模型不会覆写旧结果。

## 结构及统计图像查看

```text
neo-mri-view --background /analysis/T1w.nii.gz --space native-T1 --output /analysis/T1_view.html
neo-mri-view --background /analysis/MNI_T1w.nii.gz --overlay /analysis/z.nii.gz --space MNI152NLin2009cAsym --threshold 2.3 --output /analysis/activation_view.html
```

输出使用 Nilearn 的 HTML 交互式正交切片。结构和叠加必须已经位于相同网格；查看器拒绝隐式配准。显示阈值不等于统计显著性，空间名称由用户明确声明，软件不从文件名猜测。本预览没有移植病例定制的四维滑块及病灶联合浏览器。

## 使用与验证边界

- [独立输入、预处理和任务分析说明](docs/USAGE.md)
- [验证范围和首版限制](docs/VALIDATION.md)
- [协作与模块接入](CONTRIBUTING.md)
- [代码来源、依赖及待确认的许可](NOTICE.md)

CFM/SDC 的输入检查、规划和配准审计是可选功能；Dice 只是逆一致性指标，不能独自证明解剖定位准确。不要按激活体素最多挑模型。本项目不自动计算 CST 健康度或临床入组资格，也不包含完整静息态 FC/ICA 去噪流程。

尚未核定代码公开许可；本候选没有擅自赋予 MIT/Apache 等授权。首次公开推送前由维护者确认来源及公开范围。
