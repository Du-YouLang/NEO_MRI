# 人工病灶接口（本地预览）

独立于 BOLD、fMRIPrep 和任务时序；没有模型、权重、自动分割、自动配准或自动病灶判读。当前接入 Slicer 窗口，不是网页画笔。

第一次使用请先读[3D Slicer 人工病灶勾画指南](SLICER_GUIDE.md)，含界面截图、三平面浏览、Paint/Draw/Erase、保存草稿、退出及常见问题。

如需模型提供初始候选，另有[SynthStroke 单一模型家族接口](SYNTHSTROKE.md)。它使用外部环境、源码和权重，输出仍需返回本人工接口检查和审核。

T2/FLAIR 背景复用现有 [T2_structural_v1](T2_STABLE_V1.md) 或 [FLAIR 研究入口](T2_FLAIR_EXPERIMENT.md) 的已核查结果，不新增配准。CFM 自动接入与 MNI 功能统计图工具不在本次结构像扩展范围。

## 安装与配置

在 Python 3.12 环境按主 README 安装。3D Slicer 另行安装，在配置中填写可执行文件。代码和患者工作目录分离，结果必须放在代码仓库外。

复制 `configs/examples/lesion.example.json` 到外部工作目录，填写后移除 `template_only`。相对路径默认以该配置文件所在目录为基准；指定 `--root` 时统一以该目录为基准，配置文件参数也以此解析。

- `reference_t1`：最终选定的完整头部 T1 NIfTI（可为原始或已校正 N4 的 T1，明确记录选择）；这是所有输出的参考几何。
- `review_images`：可为空；每项为 `modality`（`T2w` 或 `FLAIR`）、`image`、`space`（`selected-T1w`）。**首版只接收已配准、同网格的背景**。原始T2/FLAIR需在外部完成配准、保存变换并人工核查后再输入；本模块不静默配准/重采样。
- `initial_mask`：空值表示从零勾画，也可指向已有同网格 NIfTI。默认仅接受0/1；多标签须明确填写正整数 `initial_label`。
- `lesion_kind`：`stroke`、`wmh`、`other`；不同类型使用不同工作目录。
- `workspace`：独立版本目录，不可放在代码仓库、原始输入目录内；不能直接接管旧实验结果。
- `editor`：`{"type":"slicer","executable":"完整路径/Slicer.exe"}`。

输入应有有效空间映射，qform/sform不能冲突。没有空间头信息或网格不同将停止，不能仅替换affine绕过检查。初始掩膜和编辑结果不自动与脑组织掩膜求交。

所有输入必须显式声明 NIfTI 长度单位为 `mm`；未知、米或微米单位会停止。需要转换时必须正确处理坐标、affine 与空间头，不能仅改单位字段冒充毫米。只有在已确认数值坐标本来就是毫米时，才可补全缺失的单位元数据。

2026-11-05 单位校验修复对已知前一版准备记录提供精确兼容识别：仍完整验证当前输入、毫米单位和冻结文件哈希，符合条件的既有分割／审核版本可继续使用，不需重画；新准备和新导入记录当前引擎哈希。

Slicer可能更改NIfTI存储轴次序/方向而保留物理采样格。本模块自身保存的草稿在导入时允许无插值的轴置换/翻转，转换后必须严格匹配T1，且记录源affine和转换状态。外部掩膜默认拒绝这类差异，明确确认后可增加`--reorient-to-reference`；该参数不接受平移、旋转配准或分辨率改变，不执行插值。

## 分阶段使用

以下命令使用独立入口 `neo-mri-lesion`，原有 `neo-mri` 结构像入口保持不变。

```text
neo-mri-lesion prepare --root /analysis --config lesion.local.json
neo-mri-lesion edit --root /analysis --config lesion.local.json --dry-run
neo-mri-lesion edit --root /analysis --config lesion.local.json
```

编辑窗口自动加载参考T1和同网格背景，进入Segment Editor。右侧NEO面板可切换背景、保存草稿。保存按钮不表示审核完成。每次保存生成独立快照，包含NIfTI和`.seg.nrrd`；重新打开恢复上次已保存草稿。首版每个会话只有一个标签，不在同一会话混画卒中和WMH。编辑期间如更改分割几何，导出仍以指定T1为参考；保存后必须回本模块检查。

```text
neo-mri-lesion import --root /analysis --config lesion.local.json
```

不指定`--mask`时导入该准备版本最后保存的Slicer草稿。也可导入外部工具输出：

```text
neo-mri-lesion import --root /analysis --config lesion.local.json --mask edited_mask.nii.gz
neo-mri-lesion import --root /analysis --config lesion.local.json --mask labels.nii.gz --label 2
neo-mri-lesion status --root /analysis --config lesion.local.json
```

导入验证三维、有限值、标签、参考网格，保存独立版本、SHA256、候选体积、多序列三平面PNG和静态`report.html`。PNG只展示代表切片，并非全脑审核证明。体积是当前掩膜描述，不是确诊病灶体积。

```text
neo-mri-lesion review --root /analysis --config lesion.local.json --reviewer "Reviewer" --review-note "Reviewed all slices against T1 and FLAIR"
```

审核是用户明确作出的记录，不是身份认证或临床有效性认证。空掩膜需额外`--confirm-empty`，区分尚未勾画和审核后未标出病灶。编辑器关窗、保存草稿或导入检查通过都不会自动审核。掩膜、T1或背景改变后，旧审核不能用于新版本。`--version`可指定先前已导入的完整64位版本键，仍须匹配当前准备输入。

```text
neo-mri-lesion export --root /analysis --config lesion.local.json --purpose analysis --export-dir exports/lesion_v1
```

只导出已审核版本；目录必须全新。`analysis`输出标准二值掩膜、审核记录与来源。空结果可导出供记录，但不能把未知损伤负荷当成0。

## 复用 T2/FLAIR 结果

以结构配置中 `t1` 指定的同一幅完整头部 T1 作为 `reference_t1`。将已核查分支的 `raw_independent_HDBET_T1.nii.gz` 或 `n4_independent_HDBET_T1.nii.gz` 填入 `review_images[].image`；显式选用 T1 脑域时使用对应 `*_T1_only_T1.nii.gz`。原始强度与 N4 的选择须记录，不自动选优。FLAIR 仍按其独立研究路线核查。

若参考是另一幅原始 T1，先核实它与上述结构参考的物理空间和网格一致；不一致时停止，不能替换 affine 冒充配准。病灶模块只接收同网格背景，不修改冻结 T2 引擎、审批或结果。

如病灶已在源 T2 上勾画，现有结构流程的 `lesion` 字段接收源网格二值掩膜，并随同一变换用最近邻生成 `*_lesion_T1.nii.gz`；核查后可作为本接口 `initial_mask`。不要把已在 T1 上勾画的掩膜填入要求源 T2 空间的字段。已有勾画只需核查、必要时修订，无需重新从零画。

导出的 `lesion_T1w_mask.nii.gz` 可供后续同一 T1 空间的分析读取，附带审核与来源记录；不会自动作为脑提取掩膜，也不会自动修改已批准 T1/MNI 或启动重算。CFM、DWI/CST 负荷和功能统计流程待相应模块恢复后另行接入。

## 输出与恢复

```text
workspace/
  prepared/<key>/     冻结参考影像、初始掩膜、session与收据
  editor/<key>/       可编辑草稿快照、last_draft.json
  versions/<key>/     标准掩膜、来源、QC、HTML与收据
  reviews/<key>/      绑定掩膜/参考身份的审核记录
  current.json        最近导入版本
```

原始输入只读。所有注册文件使用完整哈希检查；相同导入复用，新掩膜新建版本。失败的未完成目录保留，检查后换独立工作目录，不覆盖/篡改收据。

## 验证边界

见[本次验证范围](LESION_VALIDATION.md)。合成影像和文件检查不能证明真实病灶边界准确；全体积解剖复核仍由使用者负责。运行结果与患者数据均保留在代码仓库外。
