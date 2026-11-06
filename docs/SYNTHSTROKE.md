# SynthStroke 候选分割接口（本地预览）

当前只接入 SynthStroke 一个模型家族；通过配置选择 `baseline` 或 `synth`。模型源码、权重和推理依赖都放在仓库外，由用户手动提供；没有自动下载。原有纯人工流程仍可独立使用。

## 配置

复制 `configs/examples/synthstroke.example.json` 到外部工作目录，填好后移除 `template_only`。相对路径默认以配置文件目录为基准；指定 `--root` 时统一以该目录为基准，包括嵌套病灶配置内的路径。

- `lesion_config`：已有 `neo-lesion-1` 配置，`lesion_kind=stroke`、`initial_mask=null`。参考必须为选定的完整头部 T1。已有草稿使用原配置继续编辑，本命令不覆盖它。
- `workspace`：候选推理目录，必须与人工工作目录、模型、源码、运行环境和输入路径分开，也必须位于代码仓库外。
- `python`：外部 Python 可执行文件；该环境需要兼容的 PyTorch、MONAI、nibabel、NumPy、safetensors、huggingface_hub。主工作流不增加这些依赖。
- `python_paths`：通常为空；仅在使用独立依赖目录时填写完整路径。不是 pip 包名。
- `model_source`：包含官方 `synthstroke_model.py` 的外部源码目录。
- `model_dir`：包含匹配的 `config.json` 与 `model.safetensors`，不接受 `.pth`、ONNX、Hub ID 或任意其他模型。
- `variant`：`baseline` 对应 1 通道、2 输出、卒中标签 1；`synth` 对应 1 通道、6 输出、卒中标签 5。本版不接收 4 通道模型。
- `device`：`cpu` 或 `cuda`，不可用时报错，不自动切换；`threads`、`patch_size` 默认 6、128。

T2/FLAIR 当前只作同 T1 网格的审核背景，不会作为额外输入通道送进网络。

T1、背景和输出均须使用显式 `mm` 空间单位；未知或其他长度单位会停止，避免把数值相同但实际尺度不同的影像当成同一空间。

## 运行与人工流程衔接

```text
neo-mri-lesion synthstroke --root /analysis --config synthstroke.local.json --dry-run
neo-mri-lesion synthstroke --root /analysis --config synthstroke.local.json
```

dry-run 校验路径、影像几何、模型通道和文件哈希，不启动外部 Python，也不写结果；不保证外部依赖、GPU 和权重能成功加载。实际运行每次创建独立目录，不复用或覆盖旧预测。

成功返回 `CANDIDATE_REVIEW_PENDING`，输出：

- `candidate_T1w_mask.nii.gz`：原始 T1 网格的 uint8 二值候选。
- `lesion_softmax_score.nii.gz`：原始 T1 网格的类别 softmax 分数，不是经过校准的置信度。
- `request.json`、`worker_result.json`、`receipt.json`、`inference.log`：输入、模型及源码哈希、实际解释器、依赖版本和执行记录。
- `lesion.with-candidate.json`：指向候选掩膜的人工流程配置，沿用原人工 workspace；原配置文件不会改写。

将下面 `RUN` 替换为命令返回的实际输出目录：

```text
neo-mri-lesion prepare --root /analysis --config RUN/lesion.with-candidate.json
neo-mri-lesion edit --root /analysis --config RUN/lesion.with-candidate.json
neo-mri-lesion import --root /analysis --config RUN/lesion.with-candidate.json
neo-mri-lesion review --root /analysis --config RUN/lesion.with-candidate.json --reviewer "Reviewer" --review-note "Reviewed all slices against T1 and FLAIR"
```

先在 Slicer 检查、修改并保存草稿，再 import。推理成功、空预测、导入校验通过都不代表审核完成。后续使用[人工病灶接口](LESION.md)导出审核后的掩膜。CFM 自动接入留待功能流程恢复后另行接入。

## 方法及边界

沿用本地已有方案：完整头部 T1 → 保留原 affine 的 MONAI 载入 → RAS → 1 mm 重采样 → 直方图归一化 → 全图 z 归一化 → 滑窗、50% 重叠、Gaussian logits 融合，无 TTA → logits 逆变换至原 T1 → softmax/argmax。输出保存原始 qform/sform。没有隐式脑掩膜裁剪、连通域去除、T2/FLAIR 配准或自动审核。

默认 patch 128 沿用已有本地推理设置，不宣称复现论文训练／测试指标。接入可运行与患者分割准确是两个验证目标；患者结果仍需逐层检查。

来源：[SynthStroke 官方代码](https://github.com/liamchalcroft/SynthStroke)。此仓库仅提供适配脚本，不分发外部模型代码及权重。
