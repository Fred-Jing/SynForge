# SynForge 中文说明

项目用于把蛋白设计转化为可合成的基因、目标质粒及引物，再用测序证据筛选可用克隆。目标是高效、稳定的合成和 DNA 层面尽可能好的表达潜力；预测分数不能代替实验。

## 公开核心

- 结构序列提取、非标准残基处理、同义密码子搜索与表达相关指标。
- 当前专用 pET22B 骨架下的 overlap PCR / Gibson 图谱和引物设计。
- 交替 PCA 的布局搜索、热力学筛查及独立 DNA 输入入口。
- 配置驱动的 Sanger 质量、双向序列和峰图分析；导出全部达标样本，不限每构建株数。
- 全质粒环状序列比较和 BAM 检查的基础函数。
- Azenta 已适配布局的无损订单模板写入器。

PCA 的完整结构到载体流程、点突变自动设计和全质粒通用批处理入口尚未统一；请勿将批次中做过的任务当成已经产品化的功能。

## 安装与试用

在 Python 3.12、uv 环境中，从仓库根目录执行：

```sh
uv sync --locked
uv run scripts/make_demo.py --output demo
uv run scripts/analyze_sanger_batch.py --config demo/batch.json
uv run scripts/design_pca.py --target demo/pca_target.fasta --output demo/pca_layout.json
uv run scripts/design_project.py --project demo/design --chain A --map-start DEMO-1 --primer-start DEMO-1-1-1-FW --candidates 2 --workers 1
```

`demo/` 中是人工生成的序列、理想峰图及序列型 mmCIF，无真实实验结果。Sanger 预期为三个通过、一个亲本等位基因不符、一个未收到；详见 `results/usable_samples.md` 和离线 `report.html`。使用新目录运行示例，不覆盖旧结果。

默认骨架无 pelB，尾部为 `GGSGGSG-C-SLE-His6-TGA`；特殊残基默认 TAG。该骨架为特定实验体系配置，更换载体不能只改一个文件名而忽略插入位置和注释要求。

PCA CLI 在固定 DNA 上生成首轮拼接寡核苷酸，当前支持 ≤60 nt；不自动生成外侧扩增引物和载体拼接设计，不自动修改密码子。Sanger 仅确认配置指定的覆盖区域，不能称为全质粒验证或表达验证。

## 导航

|文档|内容|
|---|---|
|[英文首页](../README.md)|功能状态、安装、输入输出|
|[Sanger 工作流](SANGER_WORKFLOW.md)|样本身份、参考、质量策略及交付|
|[订单写入](ORDER_UPLOAD.md)|显式输入、模板保留与网页导入|
|[方法与限制](METHODS.md)|设计默认值、代理指标和可解释边界|
|[贡献约定](../CONTRIBUTING.md)|验证与变更规范|

运行检查：`uv run python -m unittest discover -s tests -v`。实际批次、测序原始数据、联系方式、历史订单和供应商原模板保留在实验工作区；本仓库维护可复用核心及人工示例。

自动化配置见 `ci/github-actions.yml`。当前 GitHub 连接拒绝写入 Actions 工作流，因此尚未启用远端 CI；本地 49 项测试已通过。具备相应权限后，将该文件放入 `.github/workflows/tests.yml` 即可启用。
