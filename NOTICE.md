# 来源及许可状态

本候选提取自 NEO 的 MRI 主流程，已有基础源于团队在 2026-10-27 接收的同事交接代码，并持续扩展。具体作者归属、团队授权及公开许可证尚待仓库维护者确认；未擅自设置开源许可证。公开 GitHub 仓库并不自动等于已授予开源许可。

首版不包含原交接压缩包、历史源码快照或许可证材料。不得删除实际适用的第三方署名；确认来源后应补充正式 LICENSE 及贡献者名单。

外部工具与库各自按其许可证使用，本仓库不再分发 fMRIPrep 镜像、ANTs、FSL、FreeSurfer、Nilearn 等软件。使用者需自行取得适用的 FreeSurfer license，并仅在本地挂载。

追踪入口依据 NEO 团队审核的 CST／CRP 交接协议实现。MRtrix、HMAT、XTRACT、Tang、Navigator、Boyne 和 Jülich 相关工具或图谱由使用者按各自许可准备。仓库提供合成演示。

病灶重叠入口使用外部 SMATT、HCPA 与 Boyne 2025 资产。来源和 SHA-256 见 [资产清单](lesionload/reference_assets.lock.json)。[Boyne 上游 MIT 声明](lesionload/licenses/Boyne2025_MIT.txt)随来源说明保留。

参考：[fMRIPrep 用法](https://fmriprep.org/en/25.2.3/usage.html)、[Nilearn](https://nilearn.github.io/)、[setuptools 项目配置](https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html)。fMRIPrep 链接是同系列文档参考；实际镜像仍固定为 25.2.7，并必须核实本地镜像的命令和身份。
