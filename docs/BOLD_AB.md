# BOLD v2 与 DoubleScanner A/B

这是代码检出目录的独立任务 BOLD 入口。结构与病灶 CLI 保持原接口。BOLD 处理读取经核验的 fMRIPrep 衍生文件，报告和比较只读取已保存的拟合。

## 按阶段使用

从 configs/examples 复制 task_bold_v2.example.json 或 task_bold_ab.example.json 到自己的配置目录。相对路径从配置文件的目录解析。将 data_root 指向受试者数据根目录，workspace 指向该访视 results/task_fMRI 下的新子目录。病例身份、任务标签、实际采样时间与事件应逐条确认。一个采集只选择一份 NIfTI，其他空间、回波或 part 不另算 run。每条 run 可独立暂停，只有确认同任务、同估计量和同空间的重复 run 才可联合静态 β。

- inventory 读 NIfTI 头、列 run 和缺项。
- preprocess 登记已成功的 fMRIPrep 结果与输入哈希。接受已登记格式的 completion JSON 和启动参数。STC 从实际 BOLD JSON 核对。开启 STC 时用 StartTime + k×TR。
- audit 检查事件与 TR、设计秩、VIF、自由度、运动尖峰与混杂的重叠。
- glm 生成 N0、静态 GLM＋时间导数、DoubleScanner 动态拟合和 FIR 敏感性结果。0/2/3 mm 的非零核各从同一源平滑一次。AR(1) 从未额外平滑的静态模型残差估计，在该 run 和各尺度固定。
- compare 从已保存的 β 按同一网格作描述性左减右。DoubleScanner_AB 比较回退后的 β，DoubleScanner 比较回退前的 β。更改比较不重新拟合。
- report 生成自动 HTML。正式 A/B 统计输出和报告放在指定个体访视下。4D 查看器需要显式指定同网格、已登记的 T1 及掩膜。

调用示例为 python run_bold.py audit --config path/to/local.json。可选 --verify-full 对已登记文件重算哈希。Notebook 入口见 notebooks/Task_BOLD.ipynb。真实路径、事件数据和审查证据留在本地。

## 动态 A/B 的时序与方法证据

启用 A/B 要显式设定完整原空间 R² 网格、AR(1)、Bootstrap 随机种子和执行设备。CPU 是基础实现；CUDA 重抽需要另行安装兼容的 PyTorch。A/B 在声明的主核宽执行，其他核宽仍保留 N0、GLM＋导数、原动态及 FIR。

真实 run 启用 A/B 还需登记 method_review，记录 case_id、run_id、bold_sha256、observed_stc、observed_sdc、reviewer、registration_review、status=RUN_METHOD_REVIEWED，以及已登记的衍生目录相对 report 路径和 report_sha256。registration_review 由操作者实看 BOLD→T1 及标准空间对应后填写。该回执绑定既有报告和衍生输入，现有结构配准审核工具可独立使用。本入口消费已经执行的 SDC/STC，注册阶段不启动校正或空间配准。未登记方法审核时静态路线可继续，A/B 会注明缺项。

条件授权时序需确认原始刺激日志的核验状态，并绑定逐 run 的事件与 TR 证据哈希及预处理指纹。若研究使用相对扫描启动提前 3 TR 的事件，这一约定属于该次经确认的条件参考。示例配置不会移动事件，也不会根据任务标签自动使用固定时间。

## 统计与显示

1. Lag −12 到 +10 秒，间隔 0.5 秒；时长比例 0.20 到 1.50，间隔 0.05。固定逐体素 AR(1) 做 GLS，按原空间 R² 选 1,215 个候选。先判触边，四角归时长上下边，再判原空间两项增益。相对 N0 至少减少 5% SSE，且相对 (0 秒, 1 倍) 原空间 R² 至少增加 0.01。回退后仍估计带符号 β、t 与名义双侧 p。
2. 一个 run 的整个有效掩膜是一个 BY 家族，包含扫描与回退。A 用 min(1,1214p)，B 用名义 p，另存乘 1215 的审计结果。B 的 q 和稳定标签用于探索。A 对应本项目批准的乘数约定，不能视为已经在独立样本验证的通用选后推断。
3. 对 B 通过且扫描保留的体素，共享 1,000 轮中心化白化残差重抽；A 复用其中的相同体素和轮次。Manhattan 三格内且重选仍在网格内部的比例达到 80% 才标稳定。报告保留精确回选和邻近回选频率。稳定率以当前噪声和任务参考为条件。
4. HTML 保存原始 argmax、触边和效应量回退的热图。同步查看器可切换 HRF 前后、β、t 和模型来源，时间按实际采样时钟显示。共同起点与终点独立生成描述性掩膜，橙黄色用于 T1 三视图，绿色用于 3D 点云。

## 接入验收范围

冻结样例的两条 run 分别验收。保存的 A/B 体素数为 MI 5,845／37,595，PM 3,701／25,824；对应的 A/B 稳定集合在这两条 run 内分别完全相同，为 1,581 和 1,196。将 A 从 1214 改为 1215，在 PM 少 2 个通过体素，稳定集合未变。上述结果来自只读复算。可移植入口另做合成端到端测试。单例条件参考下的激活和时间峰群不外推到新的受试者。
