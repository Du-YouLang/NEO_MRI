# SMATT、HCPA Max-WLL 与 Boyne 2025 病灶重叠

这个研究入口从同访视的原生 T2 二值病灶和已经审核的配准开始。它测量病灶与**参考图谱**的空间重叠，输出整束、轴位、视野和可选 FA。各方法使用自己的模板网格与分母。先用合成例练习，再接入已审核的个体。

## 练习

用仓库 Python 3.12 依赖和 SimpleITK 2.5.6，从检出目录运行。也可以打开[分步教程](../notebooks/Lesionload_Tutorial.ipynb)。

```powershell
python run_lesionload.py demo --output work/lesionload_demo
```

打开输出的 `results/HCPA/report.html`，在左侧 CST 查看两个有束体素的切片。第一层 5 个参考体素、5 个有覆盖、命中 1 个。第二层 10 个参考体素、8 个有覆盖、命中 3 个。合计 15 个参考体素、13 个有覆盖、4 个命中、2 个未知。覆盖内比例为 4/13，全参考束的区间为 4/15 到 6/15。最大轴位面积是 10，所以 HCPA Max-WLL 区间为 3 到 5 个体素等效。

切换侧别可比较右侧完整覆盖。左侧还含一个空束，空分母显示为空。查看 D、F 的合成 FA，注意第一层命中点有一个因拟合掩膜无效而排除。拖动轴位滑块，三维轮廓、轴位、C 和 D 同步变化。点击 F 的轴位位置可反向定位。颜色和线条用于显示，表格保留原始整数分母。

## 接入真实输入

1. 使用同访视的[冻结 T2 流程](T2_LESION_PROPAGATION_V1.md)，准备与原生 T2 同网格的二值标注。明确选择已完成的 HD-BET 分支及 raw/N4。填写 [配置示例](../configs/examples/lesionload.example.json)，其中的相对路径以配置所在目录为准。已有归档链可登记 `source`、`lesion`、`t1`、`mni`、`t1_to_t2_pull`、`mni_to_t1_pull` 和 `review`，每项形式为 `{"path":"…","sha256":"…"}`。审核材料必须对应实际链。
2. 按 [资产来源与哈希](../lesionload/reference_assets.lock.json)准备本地资产。可使用 `asset_dir` 指向同一目录，或用 `asset_paths` 按清单文件名登记现存路径。SMATT 从官方 ZIP 读取 12 个二值成员。HCPA 读取固定版本的左右二值模板。Gen6 需要固定背景及用于数值求逆的 TemplateFlow H5。Boyne 需要左右派生计数图，可从已核验资产复用，或使用下方命令重建。
3. 运行时选择一个或多个方法和一个新输出目录。每个方法生成病灶、源 T2 视野掩膜、整束和轴位 CSV，以及独立 HTML。总回执包含网格、分支、输入和输出哈希。正式报告和病例配置存到个体访视结果下。

在新计算前登记同一访视的 [空间审核回执](SPATIAL_INPUT_REVIEW.md)。启用 FA 时另提供配准和模型来源审核。已有结构批准可以在证据确实覆盖所选分支、原始病灶与变换时复用。

```powershell
python run_lesionload.py build-boyne --assets local/reference_assets --output work/boyne_count_maps
python run_lesionload.py run --config lesionload.local.json --output work/lesionload_subject
```

重建 Boyne 计数图需要两份原始 TRK、18 张 ROI 与固定作者脚本文件供溯源。计算按 TRK 几何解析每条流线，并按线段穿越体素计数。同一参考参与者对一个通道的同一个体素最多贡献一次。各通道阈值为至少 1/50 位参与者。计数图供重叠入口读取时应与锁定图像一致。通道顺序见 [Boyne 源码](../lesionload/boyne.py)。

## 统计和空间口径

- **空间链。** HCPA 直接使用 MNI152NLin2009cAsym 1 mm。SMATT 和 Boyne 使用 FSL MNI152NLin6Asym 1 mm，先在目标网格求解指定 TemplateFlow 变换的数值逆并核查点残差和图谱支持域的 Jacobian，再接上 MNI2009c→T1→原生 T2 的 LPS pull。原始标签与同网格全一 FOV 分别使用同一组合链做一次最近邻采样。
- **视野头信息。** 验证原始 T2 与标注同网格后，SMATT/Boyne 沿用 T2 强度像的视野头，HCPA 沿用原标注的视野头。这个选择保留三例归档的边界体素行为。
- **覆盖。** 设参考束体素数 N、覆盖数 C、病灶命中数 H、未覆盖数 U=N−C。完整覆盖报告 H/N。覆盖不全同时报告 H/C 和全参考区间 [H/N, (H+U)/N]。轴位峰值的下界和上界分别跨层取最大值，并保留正下界并列层及分母。并集用逻辑 OR。
- **HCPA Max-WLL。** 以侧别各自最大束截面积乘以各层命中比例，再跨层取最大值。部分覆盖分别计算上下界。单位为 1 mm 网格的体素等效。SMATT 与 Boyne 使用各分束面积和体积口径。
- **可选 FA。** 在配置的 `fa` 中登记 `b0`、`fa`、`md`、`fit_mask`、`flirt_dwi_to_t1`、`flirt_inverse`、`reviewed_mask_in_t1`、`review` 文件记录，另填 `t1_sha256`、`valid_native_voxels`、`world_forward_ras` 及根级 `affected_side`（Left 或 Right）。可登记 `model_receipt`。仅在原生 DWI 做一次三线性取样，要求完整插值支持、有限且在 [0,1] 的 FA、正且有限的 MD 和有效拟合域。SMATT 分别报告病灶内、病灶外、对侧同名束。Boyne 的左右是皮层起源，分别报告各通道的覆盖内、病灶内和病灶外。

## DLPFC 待解决

现有 Boyne 方法 C 皮层分束是按作者规则和可用掩膜做的本地重放。等待作者原版左右三维体积掩膜，或经作者确认的 Gen6 白质面及生成链。补齐后重新核对 FEFarea、Other、Unassigned 的排他归属和整束计数。整束当前使用作者公布的全部流线。分束输出是研究预览。

[SMATT 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5907352/) · [HCPA 上游实现](https://github.com/willi3by/pyLesionLoad/tree/a4563f88bf310f4b2fb139d507061030de2328c5) · [Boyne 方法论文](https://pubmed.ncbi.nlm.nih.gov/40816536/) · [Boyne 公开资料](https://doi.org/10.17605/OSF.IO/Q3GW7)
