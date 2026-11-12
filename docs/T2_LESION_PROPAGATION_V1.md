# T2 原生病灶传播 v1

该入口复用 [T2_structural_v1](T2_STABLE_V1.md) 已完成分支的刚性 T2／T1 配准与同访视审核过的 T1／MNI 标准化。输入是原生 T2 上的三维二值病灶标注。运行只读既有影像和变换，另建传播结果目录。原生病灶分割、病灶排除后重新配准和 SMATT 负荷各自独立。

## 操作

复制 [配置示例](../configs/examples/t2_lesion_v1.example.json)，填写路径和要复用的 `branch`。HD-BET 是 T2 流程的默认脑域，`raw_independent_HDBET` 和 `n4_independent_HDBET` 按已审结果显式选择。已有 T1 脑域对照也可显式指定。已迁移的结构输入可以在 `input_paths` 中按原 `request.json` 的字段重绑路径，文件内容的 SHA-256 必须一致。原结构 `request.json`、`receipt.json` 和所选分支的变换及配准强度图必须齐全。

```powershell
neo-mri-t2-lesion-v1 plan --config t2_lesion.local.json
neo-mri-t2-lesion-v1 run --config t2_lesion.local.json
```

未安装项目时也可以从检出目录运行 `python code/t2_lesion_propagation.py` 加相同参数。使用 Python 3.12 和 SimpleITK 2.5.6。plan 核对上游 T2 引擎、原请求指纹及相关输入和分支产物的哈希。run 还检查原生 T2 与标注的物理网格、二值体素，并检查每个产物读回的网格与体素。输出使用新目录。

## 空间约定与审核

- T1 网格上的点沿 LPS pull 变换查询原生 T2。MNI 网格上的点先回查 T1，再回查原生 T2。先组合两段变换，然后在每个目标网格对原生标签做一次最近邻采样。参考网格采用原 T2 运行实际绑定的 T1 和 MNI 影像。[SimpleITK 重采样约定](https://simpleitk.readthedocs.io/en/main/fundamentalConcepts.html#resampling)
- 用原始 T2 网格内的全一掩膜，沿同一条变换、相同插值方式生成源视野覆盖。`source_fov=1` 表示目标体素能从 T2 文件的网格取样，网格中原有的零强度也包含在内。`source_fov=0` 表示超出原始文件视野，不能计为已观察到的无病灶背景。脑域仅用于既有配准拟合。
- `T1_lesion.nii.gz`、`MNI_lesion.nii.gz` 是病灶掩膜。`T1_source_fov.nii.gz`、`MNI_source_fov.nii.gz` 是对应的源视野掩膜。`request.json` 记录变换方向、输入、SHA-256 与分支。`receipt.json` 记录输出哈希、视野覆盖率和病灶体素数。原生标注有前景且目标标签为空时会在对应空间标记。
- 打开 `overlay.html` 在原生 T2、T1 和 MNI 逐层看病灶。可在同一空间切换配准后的 T2 与背景，核对病灶与局部解剖。若用结果计算重叠，先在相同物理网格打开对应 `source_fov`，区分无覆盖与覆盖内无病灶。合成检查验证变换与输出，真实病例仍需审核。
- MNI 输出沿用上游审批绑定的模板及网格。与外部束图计算重叠前要分别核对模板版本、分辨率和仿射。需使用其他网格时，应从原生标签和原生 FOV 出发，明确对应关系后再组合采样。
