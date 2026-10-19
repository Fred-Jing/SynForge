# Sanger 结果分析与克隆挑选

本入口按环状质粒参考定位，用于局部 Sanger 验证，包括点突变后的目标片段确认；与全质粒共识/BAM 分析分别管理（本公开版提供基础函数，尚无通用全质粒批处理 CLI）。原始测序目录只读，结果保存在项目批次内。无需安装新依赖，使用已有 uv 环境。

## 一次完整运行

1. 盘点原始 AB1 和同名 SEQ，确认样本—构建—菌落—测序引物对应关系。每个方向是一次独立测序，SEQ 与同名 AB1 不算两份证据。
2. 选择**实际实验用的**亲本图谱，明确预期突变和必须验证的区间。不得用后来未采用的设计替代实际参考；编辑写入配置，不能覆盖原图谱或原始测序。
3. 为该批次准备 `batch.json`，运行分析：

   ```sh
   uv run scripts/analyze_sanger_batch.py \
     --config demo/batch.json
   ```

   默认输出到配置同目录的 `results/`。目录非空时拒绝覆盖；同一配置生成的结果可显式加 `--overwrite` 重跑，也可用 `--output` 写新目录。输出不能位于原始输入目录中。

4. 先看 `usable_samples.md` / `usable_samples.csv` 的全部达标名单（每个构建不限数量），再看 `construct_summary.csv` / `report.md`，再用 `report.html` 复核推荐株、异常株、突变位点双向峰图及逐碱基质量。推荐来自筛选通过的克隆，不从待复核样本中硬选。
5. 若未覆盖、双向冲突、可疑混峰或可疑 indel 无法由现有证据排除，保留复核状态；补测或重新挑克隆后使用新结果重新分析。

## 配置约定

运行 `scripts/make_demo.py --output demo` 后，以生成的 `demo/batch.json` 为可执行示例，再替换为实际参考与样本。

- `schema_version: 1`；`batch_name` 为报告名称。
- `input_directory`：原始目录。相对路径均相对于配置文件所在目录。
- `sample_regex`：完整匹配 AB1 文件名，必须包含命名组 `sample`、`construct`、`colony`、`primer`。可映射任意批次命名，不限 WT/ARS/T7。
- `primers`：引物 ID → 实际 5′→3′ 序列。每条必须在所对应目标参考中唯一匹配，方向由匹配结果决定。重复测序可配置不同引物 ID，但重复同一样本/引物 ID 会拒绝，避免默默丢掉重复文件。
- `constructs`：构建 ID → `reference`、`edits`、`required_regions`。参考支持 SnapGene、GenBank、单条 FASTA。当前**预期编辑**只支持带明确 `ref`/`alt` 的单碱基替换；**观察到的异常**支持 SNV、插入和缺失。多点替换可列多个编辑，预期插入/缺失尚无配置入口。
- `edits`：`position1` 使用参考正链的 **1-based** 坐标。必须校验源等位基因，并处于 required_regions 内。
- `required_regions`：每段包含 `start1`、`end1`（1-based、两端包含）。跨环状起点时拆成两段，避免混淆。通常应包含完整插入片段、标签、stop 和两端连接区，不能只验证目标突变一个碱基。
- `expected_samples`：可选，用于报告未收到的样本；这类缺项不计作实验失败。必须注明预期数量的依据。
- `parameters`：覆盖默认参数中的指定项；未知参数拒绝，防止拼写错误悄悄失效。

## 默认筛选与边界

本地默认策略：每个必验位点至少一条 Q30 预期匹配，至少 90% 必验位点由两个方向 Q20 一致支持；预期编辑位点必须双向 Q30 且没有触发该位点副峰线索。参数不是实验平台通用验收标准。

- **PASS_BIDIRECTIONAL**：必验区域所有位点都有正反向 Q20 一致支持。
- **PASS_COMPLEMENTARY**：部分位点由单方向 Q30 补足，同时满足上述全覆盖、双向比例和编辑位点要求。报告须保留这种证据差别。
- **REJECT_DIFFERENCE**：有双方向高质量支持的非预期碱基或 indel；并不据此推断异常产生在哪个实验阶段。
- **INCOMPLETE_COVERAGE / SINGLE_STRAND_OR_LOW_Q**：不能确认，不等于克隆一定错。
- **REVIEW_TRACE_OR_DISCORDANCE / REVIEW_INPUT_OR_TRACE**：双向冲突、双向同位点副峰、输入/峰图等疑点，需复核。
- **REVIEW_OUTSIDE_REGION**：目标区满足筛选，但读到的外侧区域有 Q30 差异线索，暂不优先使用；单向外侧差异不当作已证实突变。

峰高比是启发式线索，不是纯度或等位基因频率测量。单向副峰与同向系统性背景不能替代双向证据。选择排序依次考虑通过类别、外侧疑点、副峰线索数、双向覆盖比例和质量；排序差异不等于表达差异。

引物后窗口做保留方向的半全局比对；质量修剪和坐标转换始终保留原始 trace 索引。两端未测到的参考片段不当作缺失，内部 gap 单独检查。重复区 indel 可能有等价坐标，未统一左归一化；位置不一致不能强行合并为双向确认。

## 输出与复核

- `report.html`：离线可交互报告，无外部脚本或网络请求。点样本行后可切换预期突变/异常位置，查看原始方向峰图及已统一到参考方向的碱基表。
- `usable_samples.csv` / `usable_samples.md`：所有通过筛选的可用样本，数量与全部 PASS 状态严格一致；排名不影响纳入。
- `samples.csv`、`construct_summary.csv`、`reads.csv`、`variants.csv`：逐克隆、逐构建、逐读段和差异记录。
- `evidence.json.gz`：完整逐碱基质量、方向、原始读段索引和峰形筛查值。
- `expected_references.fasta`：**预期**完整参考序列，不是实测共识。
- `observed_validated_regions.fasta`：实际读段支持的必验位点集合；未知/冲突为 N，删除位点为 N，插入需查 variants。多个必验区间按坐标连接，不能拿它冒充全长质粒。
- `manifest.json`：配置、原始 AB1/SEQ、参考和代码 SHA-256；软件、参数、样本计数和缺失项。运行后复查原始输入及参考哈希。

最后人工复核推荐株的突变位点和高风险区，并按实际验证范围表述结论。单对 T7/T7T 只确认其高质量覆盖范围；整质粒 PCR 可能在其他位置产生错误，若需要全质粒确认应使用全质粒测序。

本次实现的针对性测试覆盖方向/环状坐标、参考编辑、质量不足、不用参考填空、双向冲突、混峰、indel 和外侧差异；统一运行 `uv run python -m unittest discover -s tests -v`。
