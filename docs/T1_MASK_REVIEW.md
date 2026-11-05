# T1脑提取、人工关口与MNI闪烁核验

## 为什么加入人工关口

脑提取是生成包含脑组织的mask，不等于灰质/白质/CSF组织分割，更不等于已经重建皮层表面。现有ANTs/OASIS在少数病例失败不代表ANTs在所有病例都不适用；我们以同一个体、同一完整头部T1的对照和人工观察选择方法。

新流程：N4全头T1 → 脑掩膜候选 → 原生网格红线逐层看 → 明确批准 → 导出已审核解剖输入 → 新组织分割/标准化 → MNI闪烁与解剖检查。候选不通过可换方法，最后手工修正。没有独立参考时不报告自动准确率，不按体积最大、与旧mask最相似来选赢家。

新入口在独立NEO_MRI；旧主工作区新preprocess启动器已暂停，已有输出登记、GLM和报告仍可读取。未修改原GLM数值引擎或历史结果。现有结果不会被新审核记录追溯改写为通过。

## 候选算法与预算

2026-11-04用户指定：未来首轮默认 **HD-BET 2.0.1 CPU noTTA**（保持本轮实测配置），其他显式可选方案为SynthStrip、ANTs/OASIS30ANTs和ANTs/MNI152NLin2009cAsym。`neo-mri-extract policy`可读取当前选择规则；配置不填method时选择HD-BET，已有显式method继续生效，不追溯改写旧结果。查看器在存在HD-BET候选时优先显示它，但不会自动批准。

不通过人工核验时，再明确选择其他方法，不因缺模型、超时或mask体积/Dice等指标自动换算法。其他三项没有自动后备顺序。新受试者仍须原生脑掩膜审核再进入MNI；本轮九份MNI的用户批准不转移给新数据。

ANTs的两个命名选项要求runtime.template_id与template/prior/registration_mask三份资源身份一致。原始OASIS结果由fMRIPrep/NiWorkflows生成，通用antsBrainExtraction.sh＋OASIS适配器不冒充原封不动复现该后端；配置/命令中会保留实际脚本和资源SHA。提取模板选择与最终MNI标准空间是两个独立参数。

配置的method可取`hd_bet_cpu_no_tta`（省略时默认）、`synthstrip`、`ants_oasis`、`ants_mni2009c`。旧的显式`ants_template`仍可读取。runtime只填写选中方法的资源；也可在runtimes对象内按method分别配置四套资源。SynthStrip需要python/script/model；HD-BET需要python/model_dir；ANTs需要script/template/prior/registration_mask，以及两个命名分支对应的template_id。资源文件保留TemplateFlow的`tpl-<template_id>_...`命名；修改method不能代替修改对应资源。任何示例路径都需在本地配置为真实路径，不会自动安装或启动环境。

|方法|用途与成本依据|本流程策略|
|---|---|---|
|SynthStrip|三例CPU实测10–14秒，安装下载另计|可选，默认边界/含CSF策略不擅改；记录模型/脚本或不可变镜像身份|
|HD-BET 2.0.1 CPU noTTA|三例CPU实测30–144秒，安装下载另计|首轮默认，不混同GPU默认TTA；推理4线程，预处理/导出各1进程|
|ANTs+替代模板|模板/先验可能改善特定病例，但需图像、概率先验和registration mask相互一致|30分钟/例上限；不使用debug/sloppy求速度，不把模板mask仅改文件名当作有效替代|
|ITK-SNAP手修|人操作时间差异很大，不能保证30分钟完成|最后关口，原生网格二值mask，修后仍须审核|

上述是官方参考量级而非本机保证。一次环境安装、模型/镜像下载应单列；实际算法运行单例上限1800秒，超时记录未完成，不自动放宽预算或阈值。算法切换应由审阅决定，不能把网页上的推荐顺序冒充三例比较已证明的最优顺序。

SynthStrip提供--no-csf和边界距离参数，但收紧mask可能再伤皮层；本轮默认模型不自动加这些参数。SynthSeg用于组织标注而非直接替代本轮脑掩膜对照，暂不引入更多模型增加解释负担。

## 1．生成同图候选与比较页

`neo-mri-extract plan --config extraction.local.json` 检查算法、源/模型身份和命令，不启动环境。准备好独立Linux Python/工具与模型后才用 `run`。模型安装不合并进核心requirements-lock；HD-BET官方依赖单独环境管理，不能让pip顺手替换核心统计环境。

候选CLI支持SynthStrip官方Python脚本、HD-BET库适配器、antsBrainExtraction.sh。它不启动WSL/Docker或下载依赖；当前本机容器执行控制器另存本地实验工作区，不进入公共仓库的个人配置。run只向新目录写结果，不自动覆盖、重试或批准。外部调度仍需协调CPU/内存和共享锁。

```text
neo-mri-extract plan --config extraction.local.json
neo-mri-extract run --config extraction.local.json
neo-mri-mask compare --config mask_review.local.json
```

comparison.json将每个候选绑定到完整头部T1与mask的SHA；review.html在同一灰度窗内切换候选，红色轮廓/填充，完整头部/保留/排除，以及相对第一项的新增绿、移除紫差异。保持原生空间分辨率，仅显示强度作8bit编码，不能用已经置零的skull-stripped图像检查被误删的皮层。

认真检查额颞极、顶部硬膜/静脉窦、枕极、颅底、小脑后下缘、脑干、病灶邻近及M1。差异图中的“新增”不等于正确；可能是补回皮层，也可能多留硬膜。红线本身不是病灶。

## 2．明确记录人工判断

如果未通过：

```text
neo-mri-mask decision --manifest /analysis/review/comparison.json --candidate synthstrip --verdict NEEDS_MANUAL_EDIT --reviewer your-name --notes "Posterior border needs correction" --output /analysis/review/decision-v1.json
```

确实逐层确认可接受后，才明确执行：

```text
neo-mri-mask decision --manifest /analysis/review/comparison.json --candidate synthstrip --verdict APPROVED --reviewer your-name --notes "Reviewed all planes; cortical coverage and posterior boundary acceptable" --reviewed-all-slices --output /analysis/review/approved-v1.json
```

这不是自动生成用户批准的命令模板；reviewer/notes要填写实际判断。不要求软件通过“滑动了多少层”猜测是否认真审核。批准绑定T1、mask和比较清单；任何字节变化使旧批准失效，不覆写批准JSON来接受新mask。

如果用户在实际会话中直接明确批准候选，管理端可通过`decision(..., user_statement=原文)`记录`EXPLICIT_USER_APPROVAL`。这时可以如实保留`reviewed_all_slices=false`，而不是代替用户声明已逐层看完；同样绑定文件SHA，仅批准该掩膜用于下一步计算，不表示MNI配准或临床结论已经通过。没有明确批准原文时仍执行原有完整逐层确认关口。

比较页JSON路径必须在运行该命令的操作系统中可访问。Windows G:路径与Linux /mnt/g路径不能直接互用；同一次审核/导出/执行使用同一个明确路径环境。本版不提供自动跨系统重写冻结审核记录。

## 3．手工修正首选ITK-SNAP

对本需求的单个二值mask，ITK-SNAP的三面逐层画笔、擦除及3D检查足够直接；这是操作适配性的选择，不声称它在所有分割任务都优于3D Slicer。日后多标签/复杂表面编辑可再用Slicer。

```text
neo-mri-mask manual-export --manifest /analysis/review/comparison.json --candidate synthstrip --output /analysis/manual-v1
```

导出独立编辑副本mask_to_edit.nii.gz，原候选不动。manual_edit.json记录主T1路径及父mask。

1. 在ITK-SNAP把完整头部原生T1作为Main Image打开，不加载MNI T1。
2. 加载mask_to_edit.nii.gz为Segmentation；标签1是保留脑，0是背景。
3. 三个方向确认位置；用画笔补回误删组织、以0擦除硬膜/脂肪等；不要改变影像网格或重采样T1。
4. 用三面和3D检查连续性，另存mask_corrected.nii.gz，不能覆盖原始mask。
5. 回读检查：

```text
neo-mri-mask manual-import --bundle /analysis/manual-v1/manual_edit.json --edited /analysis/manual-v1/mask_corrected.nii.gz --output /analysis/manual-review-v1
```

回读拒绝错网格、非二值、空/全掩膜；通过这些只是格式有效，仍生成新的PENDING比较页。再次人工批准后才能用于标准化。不要用自动膨胀/闭操作代替尚未确认的修正。

## 4．把已批准mask交给下一轮fMRIPrep

```text
neo-mri-mask export-reviewed --review /analysis/review/approved-v1.json --raw-t1 /analysis/project/bids/sub-demo/anat/sub-demo_T1w.nii.gz --output /analysis/reviewed_anatomy-v1
```

仅导出对应的desc-preproc_T1w和desc-brain_mask及审核绑定，不带旧dseg/probseg或旧变换，以免把已换掩膜的结果与旧组织分割混用。只支持审核的单T1与原BIDS T1同网格；多T1合成或不同网格需要独立评估，不能隐式缩放修复。

阶段配置增加anatomy_review、preprocessing.reviewed_anatomy_dir。新容器命令以只读--derivatives导入；在本轮没有重新运行正式fMRIPrep/MNI，所以还需实际确认固定25.2.7拾取了正确的预计算mask。不能仅凭命令含--derivatives就宣称采用成功；运行方法报告须显示预计算脑掩膜，并核对输出mask身份/几何。旧标准化结果不能直接改标为使用了新mask。

## 5．MNI闪烁页：先代码测试，后真实验收

```text
neo-mri-blink --config blink.local.json
```

配置需绑定已批准的原生T1/mask、标准化T1、准确模板和实际变换证据。normalization_provenance保存t1_sha256、mask_sha256、normalized_t1_sha256、template_sha256及transforms文件记录。缺批准、沿用旧mask的标准化或网格不同均停止；查看器不重新配准、不缩放脑组织来凑对齐。

输出blink.html和registration_blink.gif。同一切面交替显示模板与已归一化T1，能暂停、继续、单看模板/患者、调间隔。灰度各自作1–99%显示归一化，空间不变。先看脑室/胼胝体/丘脑等内部对应，再看中央沟、内囊邻近及病灶边缘；结合往返、自洽与Jacobian，不单凭外轮廓。

屏幕上的切片会放大或缩小以便观看，并保留体素物理间距决定的纵横比。这只是显示尺寸，不是改变三维空间配准或把患者脑缩放到模板。

本轮真实T1掩膜尚未获正式批准，因此只给合成例子，不启动真实MNI验收；合成APPROVED记录明确为SYNTHETIC_TEST_ONLY，不得用于患者。

## 参考资料

- SynthStrip官方方法、CPU耗时与参数：https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/
- HD-BET官方CPU/noTTA说明：https://github.com/MIC-DKFZ/HD-BET
- ITK-SNAP手工分割教程：https://www.itksnap.org/pmwiki/pmwiki.php?n=Documentation.TutorialSectionManualSegmentation
- fMRIPrep结构与配准查看教程：https://andysbrainbook.readthedocs.io/en/stable/OpenScience/OS/fMRIPrep_Demo_3_ExaminingPreprocData.html
