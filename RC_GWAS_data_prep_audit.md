# RC（残粒胆固醇）双 GWAS 数据准备审计报告

审计日期：2026-10-18　对象：`RC.R` 第 30–103 行（data source and format + meta 读取）
结论一句话：**两个输入 GWAS 的身份、N 口径、坐标系统（GRCh37）都是对的，已发表结果是自洽的；但当前 `RC.R` 跑不起来（3 处硬错误），而且现存 txt / METAANALYSIS2.TBL / 当前代码三者互不自洽，用现在的代码无法复现已发表结果。**
补充（第二轮）：等位基因方向已抽查完毕，**没有发现问题**；同时新查出一个独立隐患（Karjalainen 9,049 个重复 rsid，影响 19 个显著 SNP、0 个 lead）并搞清了 METAL 的等位基因规范化规则（§12）。

---

## 1. 输入文件身份核实

| 项目 | GCST90302074 | GCST90499341 |
|---|---|---|
| 研究 | Karjalainen MK et al., **Nature 2024-03-06**, PMID 38448586，"Genome-wide characterization of circulating metabolic biomarkers" | Zoodsma M et al., **Nat Genet 2025-10-03**, PMID 41044249，"A genetic map of human metabolism across the allele frequency spectrum" |
| 性状 | remnant cholesterol（残粒胆固醇） | remnant cholesterol（同一性状 ✓） |
| 样本 | 4,435 东亚(中国) + 11,340 南亚(英国) + 120,241 欧洲 = **136,016** | **434,646** 英国白人（UK Biobank） |
| 队列 | 33 个队列（24 项研究；LOLIPOP 8 个子队列、Rotterdam 3 个子队列）；**不含 UK Biobank**（论文原文："absence of UK Biobank from our GWAS meta-analysis"，UKB 只用于重复验证/MR） | UK Biobank 单队列 |
| 文件行数 | 13,354,183 | 11,375,404 |
| 坐标 | `variant_id` = **GRCh37**；`base_pair_location` = **GRCh38**（99.13% 行不一致） | `base_pair_location` = **GRCh37**（与 rs_id 内嵌坐标 96.6% 一致） |
| 样本量列 | 无 → 代码硬编码 136,016（正确） | 逐行 `n` ∈ {434,646: 11,021,840 行；434,219: 353,564 行} |
| rsid | 完整（0 缺失） | 完整（0 缺失；5.25% 为 `chr:pos_alleles` 形式 ID） |

**两个队列无样本重叠**（Karjalainen 发现集明确排除 UKB），这一点可以直接写进方法学。

---

## 2. 阻断性缺陷：3 处 `[列名 != 23, ]` 在 R 里永远报错

| 位置 | 代码 | 实际报错 |
|---|---|---|
| L38（Karjalainen） | `Karjalainen_RC[chr != 23, ]`（在 `separate()` 之前，列还不存在） | `Error: object 'chr' not found` |
| L64（Zoodsma） | `Zoodsma_RC[chromosome != 23, ]`（列已被 L59 改名为 `chr`） | `Error: object 'chromosome' not found` |
| L102（meta 段） | `RC[chr != 23, ]` | 同类错误（会在写 RC.txt 前中断） |

**根因不是"顺序错了"，而是 R 的 `df[i, ]` 根本不在数据框环境里找列名**（只有 data.table 的 `DT[i]` 才会）。最小复现：

```r
mtcars[mpg > 20, ]        # Error: object 'mpg' not found   ← R 的真实行为
d2[chr != "23", ]         # Error: object 'chr' not found（d2 确有 chr 列）
d2[d2$chr != "23", ]      # 正常
```

因此把 filter 挪到 `separate()` 之后**也救不了**（本轮 runB v1 实测仍然报 `object 'chr' not found`）。正确写法只有三种：`df$chr`、`dplyr::filter(df, chr != 23)`、或 `setDT(df)[chr != 23]`。

**实证记录**
- `Temp/runA.log`：原样代码跑 → `Error: object 'chr' not found`，`Execution halted`（fread 之后第 6 行就死）
- `Temp/runZ.log`：Zoodsma 段 100 行最小复现 → `Error: object 'chromosome' not found`
- `Temp/test_subset.R`：`mtcars[mpg > 20, ]` 报错，证明是 R 语义而非列名问题
- `Temp/runB.log`：只把 filter 挪到 separate 之后（不改索引方式）→ 依然 `Error: object 'chr' not found`

---

## 3. 现存产物与当前代码脱节（复现性核心问题）

### 3.1 两个 txt 的行数完全可解释

| 文件 | 行数（数据行） | 来源 |
|---|---|---|
| `Karjalainen_RC.txt` | **9,094,470** | 13,354,183 行中 EAF>0.01 = 9,094,470（awk 实测，完全相等） |
| `Zoodsma_RC.txt` | **9,712,126** | 11,375,404 行中 EAF>0.01 **且 n=434,646** = 9,712,126（EAF>0.01 里 n=434,219 另有 302,773 行，被排除） |

### 3.2 但 meta 的输入不是现在这个 Zoodsma txt

`METAANALYSIS2.TBL` 的 Weight（有效样本量）分布：

| Weight | SNP 数 | 含义 |
|---|---|---|
| 570,662 | 7,647,315 | 136,016 + 434,646（两队列都有） |
| 434,646 | 2,064,652 | 仅 Zoodsma，n=434,646 |
| 136,016 | 1,438,095 | 仅 Karjalainen |
| **434,219** | **302,773** | **仅 Zoodsma，n=434,219** |

→ 已发表的 meta **使用了** n=434,219 的 302,773 个 SNP；而现存 `Zoodsma_RC.txt` 里 n 全部是 434,646（这 302,773 行被排除了）。**现存 txt 不是 meta 的输入文件。**

### 3.3 时间线

| 时间 | 事件 |
|---|---|
| 2024-06-22 | 下载 GCST90302074（GRCh38 重映射版） |
| 2025-10-11 | 下载 GCST90499341（GRCh37 版） |
| 2025-10-22 | 生成 `Karjalainen_RC.txt` |
| 2025-12-03 | 生成 `METAANALYSIS2.TBL`（输入含 n=434,219 行） |
| 2025-12-16 20:11 | 生成 `RC.txt`（11,071,938 行 / 19,414 个 P<5e-8，与手稿一致 ✓） |
| 2025-12-16 20:21 | **重新生成** `Zoodsma_RC.txt`（此时已变成只含 n=434,646） |

→ 20:21 那次重跑用的代码和 12-03 做 meta 时不同；此后 `RC.R` 又被改过（加入 L59-60 改名、L38 前移），形成今天的死代码。

### 3.4 Zoodsma 的 N 从未被赋值

当前代码 Zoodsma 段**没有** `$N` 赋值，txt 里只有小写 `n`。也就是说：**用现在的代码重跑，METAL 拿不到 Zoodsma 的样本量**（要么报错，要么按缺失处理）。旧版用的是逐行 `n`。
修复：`Zoodsma_RC$N <- as.numeric(Zoodsma_RC$n)`。

### 3.5 Karjalainen 的 N 是对的

`4,435 + 11,340 + 120,241 = 136,016`，与论文 "up to 136,016 participants"、与 txt 里 N 列（全部 136,016）、与 meta weight（136,016）**三处一致** ✓
（说明：本轮审计早期我误算成 135,981，已更正——15,775 + 120,241 = 136,016。）

---

## 4. 坐标系统：全流程 GRCh37（自洽，但很脆弱）

- **Karjalainen 文件两套坐标并存**：`variant_id` 是 GRCh37（研究原始提交），`base_pair_location` 是 GRCh38（GWAS Catalog 重映射）。
  证据：rs1988726 → dbSNP GRCh37 = 1:564,862，GRCh38 = 1:629,482；文件里 variant_id=`1_564862_C_T`、base_pair_location=629,482。
  rs141330630 → Ensembl GRCh38 = 22:50,795,620 = 文件 base_pair_location；variant_id 22_51234048 为 GRCh37。
  13,354,183 行中 **13,238,498 行（99.13%）两套坐标不一致**，偏移随染色体/区段变化（liftover 特征：chr1 前段 +64,620；chr22 后段 −438,428）。
- **Zoodsma 文件是 GRCh37**：rs370652263 文件 22:51,237,712，Ensembl GRCh38 = 22:50,799,284（差 438,428）；rs202228854 同样差 438,428。
- **交付结果用的是 GRCh37**：19,414 个显著 SNP 溯源——
  - 16,515/16,515 个 Karjalainen 来源 SNP：`RC.txt` 的 pos == variant_id pos（GRCh37）✓
  - 19,253/19,253 个 Zoodsma 来源 SNP：`RC.txt` 的 pos == 其 base_pair_location（GRCh37）✓
  → **RC.txt 坐标全为 GRCh37，与 1kgv3/EUR、FUMA、ANNOVAR 一致** ✓

**风险提示（必须写进代码注释）**
1. 谁要是"顺手"把 `separate(variant_id)` 改成用 `base_pair_location`，Karjalainen 立刻变 GRCh38，与 Zoodsma（GRCh37）混用 → 位点整体错位（chr1 约 65 kb，chr22 约 438 kb）。
2. 谁要是重新下载 GCST90499341（Catalog 若已重映射到 GRCh38），同样会混 build。
3. 建议在脚本里显式断言：`median(base_pair_location - variant_id_pos) != 0` 且 `chr` 取 1–22（修正版脚本已内置该检查）。

---

## 5. 可复现性缺口

1. **METAL 脚本在整个项目里不存在**（`find` 无任何 `*metal*` 文件）。`RC.R` 只读 `METAANALYSIS2.TBL`，meta 本身是黑盒 → 仓库上传后别人无法复现 meta。
2. `plink_bin`、`bfile=1kgv3/EUR` 为个人机器绝对路径（已在上传包中参数化）。
3. 注释错误：Zoodsma 段注释为 "Nightingale Health 2020 met-d"，实际文件是 Zoodsma et al. Nat Genet 2025（UKB 434,646，remnant cholesterol）。

---

## 6. 口径问题（不阻断，但建议明确声明）

| # | 问题 | 数据 | 建议 |
|---|---|---|---|
| 1 | `effect_allele_frequency > 0.01` 是单侧过滤，保留了 EAF>0.99（MAF<1%）的 SNP | Karjalainen 22,063 个、Zoodsma 19,964 个；meta 层 `eaf > 0.01` 同样单侧 | 若本意是 MAF>1%，应写 `>0.01 & <0.99`；**改动会改变结果，需明确声明后再改** |
| 2 | `pop="EUR"` + 1kgv3 EUR 参考，用于含 4,435 东亚 + 11,340 南亚的 meta 做 clumping/LD | — | 在方法里说明"LD 参考为欧洲人群"，或对东亚/南亚 lead 做敏感性分析 |
| 3 | 两队列等位基因方向：16,384 个 meta 显著且双队列共有的 SNP 中，同向 16,253、**互换 107**、**不一致 24**（例：rs2074550 Karjalainen T/A vs Zoodsma C/A）；队列内显著交集 6,385 中 18 互换、2 不一致 | 见 `allele_orientation_flags.csv` | **已查清，无错误**：24 个"不一致"是标记脚本的假阳性，107 个互换 METAL 全部正确处理，131/131 的 meta `Direction` 与实际 beta 方向一致。详见 §12 |
| 4 | `pval == 0` 兜底为 2.2e-308 | Karjalainen 1 行、Zoodsma 13 行 | 已处理 ✓ |
| 5 | **se/beta 公式多乘了一个 f**：`1/sqrt(2*f*(1-f)*f*(N+Z²))`，标准式 `1/sqrt(2*f*(1-f)*(N+Z²))`（Var(基因型)=2f(1-f)，只应出现一次） | 全表 | z = beta/se 不变 → **p 值、clumping、FUMA/MAGMA/PoPS/FLAMES/fine-mapping 都不受影响**；但 beta、se 的绝对量纲同时放大 1/√f（f=0.25 时 2 倍）→ **SMR 输入、ANNOVAR INFO 里的 BETA/SE 会偏**。请核对 Zhu 2016（PMID 27019110）原文；确认后改一行，SMR 结果需重跑 |
| 6 | 两队列无样本重叠 | 论文明确排除 UKB ✓ | 可正面写入方法 |

---

## 7. 手稿层面（与本轮数据核对相关）

1. **"251,094（UKB + Karjalainen）"是错的**：251,094 = 136,016（Karjalainen 发现集）+ 115,078（Karjalainen 论文里 UKB **重复验证**样本量）。UKB 又正是 Zoodsma 队列 → 既重复计数，也与"发现集不含 UKB"的事实矛盾。meta 真实 N = 136,016 + 434,646 = **570,662**（= TBL 最大 Weight）。
2. **"European ancestry"** 与 Karjalainen 多祖源（东亚 4,435 + 南亚 11,340）不符 → 建议写 "predominantly European-ancestry (93.0% European, 8.3% South Asian, 3.3% East Asian)" 或列出祖源构成。
3. 已核对一致、无需改动的数字：**19,414** 个 P<5e-8、125 个本地 lead、FUMA 112 个位点。
4. 【2026-10-18 更新】`RC.txt` SNP 总数 **要改**：改双侧 EAF 后由 11,071,938 变成
   **11,058,177**（−13,761）。其余数字（19,414 / 125 / 112）不变，因为这 13,761 个
   MAF<1% 的 SNP 里没有一个是显著的。手稿改动见 §11.4 的"最终决策"。

---

## 8. 修复清单

**P0（不改就跑不起来 / 结果不可复现）—— 全部已修**

| # | 位置 | 修法 | 状态 |
|---|---|---|---|
| 1 | L38 | `Karjalainen_RC[Karjalainen_RC$chr %in% as.character(1:22), ]`，移到 `separate()` **之后** | ✅ 已改 RC.R |
| 2 | L64 | `Zoodsma_RC[Zoodsma_RC$chr %in% as.character(1:22), ]`（改名之后用 `chr`） | ✅ 已改 RC.R |
| 3 | L102 | `RC[RC$chr %in% as.character(1:22), ]` | ✅ 已改 RC.R |
| 4 | Zoodsma 段 | 加 `Zoodsma_RC$N <- as.numeric(Zoodsma_RC$n)` | ✅ 已改 RC.R |
| 5 | n=434,219 | 保留（`KEEP_N434219 <- TRUE`），Zoodsma_RC.txt 目标行数 **10,014,899** | ✅ 已定为默认 |

**P1（复现性）—— 已补**

6. ✅ METAL 脚本已反推并写入 `metal_rc.txt`（原脚本确认全项目不存在）。详见 §11。
7. ✅ plink/bfile 已参数化（上传包）；Zoodsma 注释已修正。
8. ✅ **genomic control 步骤已补**（2026-10-18，方法由用户确认 = `median(χ²)/0.454`）。
   实测 λ_Kar = 1.1488、λ_Zoo = 1.2277，与 TBL 反推值差 < 0.7%；meta 后 λ = 1.0245。
   已写入 `RC.R`（`gc_adjust()`）与 `RC_data_prep_fixed.R`（`apply_gc()`，开关 `GC_APPLY`）。

**P2（口径）—— 已处理**

9. ✅ se 公式已改：`1/sqrt(2*Freq1*(1-Freq1)*(Weight+Zscore^2))`（删掉多余的 `Freq1*`）。
10. ✅ EAF 已改双侧（`EAF_TWO_SIDED <- TRUE`），决策与实测见 §11.4。
11. ✅ 手稿数据集描述、祖源构成、meta N 已改（输出 `RC_GWAS_manuscript_updated_20261018.docx`）。
12. ✅ **107 个 swapped + 24 个 mismatch 位点已抽查**（2026-10-18）。结论：24 个 mismatch 全为假阳性，107 个互换 METAL 全部正确处理，2 个标记位点是 lead SNP 且方向无误。同时新查出一个真实的独立隐患（Karjalainen 9,049 个重复 rsid）并搞清了 METAL 的等位基因规范化规则。详见 **§12**。

---

## 9. 5 个问题的决策与执行结果

| # | 问题 | 你的决定 | 执行情况 |
|---|---|---|---|
| 1 | Zoodsma 302,773 个 n=434,219 SNP | **保留** | ✅ `KEEP_N434219 <- TRUE`；`Zoodsma_RC.txt` 应为 10,014,899 行 |
| 2 | METAL 脚本 | **不在** | ✅ 已反推出 `metal_rc.txt`（§11） |
| 3 | se 公式 | **改** | ✅ RC.R L76 已改；⚠️ SMR 输入与 ANNOVAR 注释里的 BETA/SE 需重跑 |
| 4 | EAF 双侧 vs 单侧 | **双侧** | ✅ `EAF_TWO_SIDED <- TRUE`；`RC.txt` 行数 11,071,938 → 11,058,177 |
| 5 | 手稿 N 与祖源 | **按报告改** | ✅ 已改并另存新文件 |
| 6 | GC 的 λ 算法（追加提问） | **`median(χ²)/0.454`** | ✅ 已验证 λ 1.1488 / 1.2277，已写入代码 |

---

## 10. 证据文件

| 文件 | 内容 |
|---|---|
| `Temp/runA.log` | 原样代码运行 → `object 'chr' not found` |
| `Temp/runZ.log` | Zoodsma 段最小复现（100 行） |
| `Temp/test_subset.R` | `mtcars[mpg > 20, ]` 亦报错，证明 R 语义 |
| `Temp/runB.log` / `runB2.log` | 修正版复跑（v1 证明换顺序不够；v2 用 `$` 索引） |
| `Temp/awk_checks.log` | txt 的 N 列（全 136,016）、chr 列（1–22，无 23/X）、Zoodsma n 列（全 434,646） |
| `Temp/awk_zcrosstab.log` | EAF×n 交叉表（9,712,126 = EAF>0.01 ∩ n=434,646） |
| `Temp/awk_offset.log` | variant_id 与 base_pair_location 偏移直方图（liftover 特征） |
| `Temp/awk_zconsist.log` | Zoodsma rs_id 内嵌坐标与 base_pair_location 一致性（96.6%） |
| `Temp/sig_rc.txt` / `sig_kar.txt` / `sig_zoo.txt` | 19,414 个显著 SNP 的坐标溯源原始数据 |
| `allele_orientation_flags.csv` | 19,414 个显著 SNP 的双队列等位基因方向标记（same/swapped/mismatch/单队列）——**其中 24 个 mismatch 已被 §12 证伪，勿据此改数据** |
| `Temp/audit_alleles.py` | 第一轮方向审计：131 → TBL 交叉核对、lead 交集、meta Direction 一致性 |
| `Temp/audit_alleles2.py` | 第二轮行感知重标：按"是否存在相同/互换等位对"重新判定，并统计两文件重复 rsid |
| `Temp/data_checks.log` | 两源文件逐项统计（EAF 极值、pval=0、rsid 缺失、显著 SNP 数等） |
| `Temp/gc_impact.py` | 逐位点重建 GC 后的 meta Z，量化"修复后已发表结果会变多少"（§13） |

---

## 11. METAL 反推（2026-10-18 补充）

原始 METAL 脚本确认在项目内不存在（`find` 无任何 metal 文件），因此从 `METAANALYSIS2.TBL`
的**输出指纹**反推配置，结论写入 `metal_rc.txt`。

### 11.1 从 Weight 列反推 scheme

TBL 的 Weight 列**只有 4 个取值**，且全部是 per-marker N 的整数和：

| Weight | 行数 | 解释 |
|---|---|---|
| 570,662.00 | 7,647,315 | = 136,016（Karjalainen）+ 434,646（Zoodsma） |
| 434,646.00 | 2,064,652 | 仅 Zoodsma，n=434,646 |
| 136,016.00 | 1,438,095 | 仅 Karjalainen |
| 434,219.00 | 302,773 | 仅 Zoodsma，n=434,219 |
| **合计** | **11,452,835** | |

两条硬结论：

1. **不可能是 `SCHEME STDERR`**。STDERR 方案下 Weight 是 Σ(1/SE²)，量级 1e4–1e5、非整数、且每个
   marker 都不同。只有 p 值方案（`SCHEME SAMPLESIZE`）+ `WEIGHT` 指向样本量列才会得到 ΣN。
   这与 RC.R 里 `se = 1/sqrt(2*Freq1*(1-Freq1)*(Weight+Zscore^2))` 的假设完全自洽——该公式把
   Weight 当 N 用（PMID 27019110）。
2. **302,773 个 n=434,219 的 marker 确实进入了 meta，且只来自 Zoodsma**：若它们同时存在于
   Karjalainen，Weight 会出现 570,235（= 136,016+434,219），而全表**没有这个值**。
   这正是"① 保留"这个决定的直接证据。

### 11.2 Direction 列证实输入文件只有两个、顺序为 Karjalainen→Zoodsma

- 仅 Karjalainen 有的 marker：`Direction = "-?"`
- 仅 Zoodsma 有的 marker：`Direction = "?+"`

`?` 是 METAL 对"该队列缺此 marker"的标记，不是缺失方向。两个字符 ⇒ 两个输入文件，
且第一列对应 Karjalainen。

### 11.3 新发现：genomic control 被做过，但代码里没有（已确认）

对 10 个"单队列 marker"（meta Z 必须等于该队列自身 beta/SE）做核对：

| 队列 | \|Z_tbl\| / \|beta/se\| | 样本 |
|---|---|---|
| Karjalainen | **0.9338**（6/6 一致） | rs146253013, rs79621693, rs189111605, rs144382691, rs147850602, rs144875014 |
| Zoodsma | **0.9055**（4/4 一致） | rs7929618, rs111161815, rs7648150, rs78656624 |

两个**不同**但**极其恒定**的缩放比 ⇒ 每个队列各自被乘了一个常数，而不是随机误差。
λ = 1/缩放² ⇒ λ_Kar = **1.147**，λ_Zoo = **1.220**。

双队列行独立验证：rs2326918，z_Kar = −0.2387、z_Zoo = −0.3058，按 METAL 的 √N 权重
（Σ√N·z / √ΣN）代入 GC 后的 z 得 |Z| = **0.3505**，TBL 实测 **0.351**。吻合到小数点后三位。

**λ 算法已由用户确认：`median(χ²)/0.454`**。按此公式在现有两个队列 txt 上实测：

| 文件 | n | median(χ²) | λ = /0.454 | 从 TBL 反推 | 差 |
|---|---|---|---|---|---|
| Karjalainen_RC.txt | 9,094,470 | 0.52154 | **1.1488** | 1.1468 | 0.17% |
| Zoodsma_RC.txt | 9,712,126 | 0.55737 | **1.2277** | 1.2197 | 0.66% |
| RC.txt（meta 后） | 11,071,938 | 0.46512 | **1.0245** | — | — |

两个队列 λ 都对得上（Zoodsma 的 0.66% 差异正好可由"meta 当时用的是含 302,773 行
n=434,219 的完整版文件、而现在盘上只有 9,712,126 行"解释）。链路闭合。

**实现要点**：METAL 走 `SCHEME SAMPLESIZE`，读的是 **PVALUE** 而不是 beta/se，
所以 GC 必须落在 p 值列上才起作用。代码里保持自洽的写法是：

* `beta` 不变（GC 只动检验统计量）
* `se ← se × √λ`
* `pval ← 2·Pnorm(−|z/√λ|)`

已写入 `RC.R`（`gc_adjust()`）与 `RC_data_prep_fixed.R`（`apply_gc()`，开关 `GC_APPLY`）。

**这同时修正了 §3.2 的旧结论**：不是"txt 被换过"，而是 meta 的输入是 GC 后的版本。

### 11.4 EAF 单侧 vs 双侧：实测差异与最终决策

`data source/RC.txt`（meta 后、单侧 `EAF > 0.01` 过滤后）实测：

| 指标 | 数值 |
|---|---|
| 总 SNP | 11,071,938 |
| EAF > 0.99（即 MAF < 1% 但被单侧过滤放进来了） | **13,761**（为什么它们以 eaf>0.99 的形式出现：见 §12.5，METAL 的 `ATCG` 等位规范化） |
| EAF < 0.01 | 0（上游已滤） |
| 其中 P < 5e-8 | **0** |
| 全基因组显著 SNP 总数 | 19,414 |
| 125 个 lead 中 EAF > 0.99 的 | **0** |

源文件层面实测（awk 直扫两个 .gz）：

| 文件 | 单侧行数 | 其中 EAF ≥ 0.99 | 双侧行数 |
|---|---|---|---|
| Karjalainen | 9,094,470 | 22,227 | **9,072,243** |
| Zoodsma | 10,014,899 | 19,965 | **9,994,934** |
| （其中 n=434,646 / n=434,219） | | 19,341 / 624 | 9,692,785 / 302,149 |

meta 层面对应的变化：`METAANALYSIS2.TBL` 应从 11,452,835 行变成 **11,439,074** 行
（−13,761）；四个 Weight 值本身不变，但每个 Weight 的行数会变，且**无法用手算出来**
（一个被剔除的 marker 可能只在一个队列存在，也可能两个都有），必须跑完再看。

**最终决策（2026-10-18）：改双侧**，`EAF_TWO_SIDED <- TRUE`。理由不是它改变结果
——它一个数都不改变——而是手稿写的是"MAF of less than 1% were excluded"，而单侧过滤
实际放进了 13,761 个 MAF<1% 的 SNP，方法与代码对不上。代价只是手稿里 11,071,938
这个数字要跟着改。
（若要与已发表数字逐位对齐，把开关设回 `FALSE` 即可，两者都没有统计后果。）

### 11.5 因此：需要有个清楚的重跑计划

GC 与双侧 EAF 两个开关都已进代码，因此下一步是**从头重跑**，而不是打补丁：

1. 用 `RC_data_prep_fixed.R`（`KEEP_N434219=TRUE`、`EAF_TWO_SIDED=TRUE`、`GC_APPLY=TRUE`）
   重建两个队列 txt —— 目标行数 Karjalainen **9,072,243**、Zoodsma **9,994,934**
   （含 n=434,219 的 302,149 行）；若数不对，先查开关。
2. 用 `metal_rc.txt` 重跑 meta —— 四个 Weight 值必须仍是
   {570,662 / 434,646 / 136,016 / 434,219}，总行数 11,439,074。
3. 在 meta 结果上重算 `se`/`beta`（修正后的公式）；
4. 下游全部重跑：clump → novelty → ANNOVAR → SuSiE → FUMA → FLAMES → SMR → coloc。

已发表的旧数字不会逐个位点对齐（GC + 双侧 EAF 都会动检验统计量），这是有意为之。
`metal_rc.txt` 末尾附了验收检查。

---
---

## 12. 等位基因方向抽查（2026-10-18 执行）

对象：19,414 个 P<5e-8 的 SNP 中被标为 `swapped`（107 个）和 `mismatch`（24 个）的
131 个位点。工具：`Temp/audit_alleles.py`（对 `METAANALYSIS2.TBL` 交叉核对）+
`Temp/audit_alleles2.py`（行感知重标，读**送进 METAL 的 txt** 而不是原始下载文件）。

### 12.1 结论总表

| 指标 | 数值 | 判读 |
|---|---|---|
| flagged markers | 131 | swapped 107 + mismatch 24 |
| 在 TBL 中找到 | 131/131 | 全部进入 meta 输出 |
| `swap_verified`（kar 的等位对 == zoo 的等位对互换） | **107/107** | 互换真实存在 |
| 行感知重标后仍是 mismatch | **0** | **24 个 mismatch 全是假阳性** |
| meta `Direction` 与实际方向一致 | **131/131** | METAL 的链翻转/方向处理**没有出错** |
| 其中是 lead SNP 的 | **2**（rs10231941、rs635634） | 逐一核对，均正确 |

两个 lead 的原始记录（meta `Direction` = 两队列符号拼接，`--` / `++` 均正确）：

| rsid | 位置 | Karjalainen | Zoodsma | meta Direction | 判读 |
|---|---|---|---|---|---|
| rs10231941 | chr7:100,532,540 | T/C, beta −0.00989 | C/T, beta +0.01652 | `--` | ✓ 互换已正确处理 |
| rs635634 | chr9:136,155,000 | C/T, beta −0.0618 | T/C, beta +0.0328 | `++` | ✓ 互换已正确处理 |

### 12.2 24 个 mismatch 的根因：取错行，不是数据错

`Temp/allele_flags.py` 用 `out[rsid] = (...)` 扫**原始下载文件**，同一 rsid 有多行时
**最后一行覆盖前面所有行**。GWAS Catalog 的 GCST90302074 里，一个 rsid 可以对应多个
变体行（同一位置的备选等位表示，EAF 相差很大）。

**24/24 个"mismatch"在原始文件里都有 ≥2 行**，且第一行（高频行）的等位对与 Zoodsma
**完全一致**；被误采的第二行都是低频备选表示：

| rsid | 原始文件全部行 | 被误采 | 送进 METAL 的 txt 里保留的行 |
|---|---|---|---|
| rs10738609 | G/A EAF=0.4677 ｜ T/A EAF=0.0072 | T/A | G/A（= Zoodsma ✓） |
| rs12131553 | C/A EAF=0.1637 ｜ T/A EAF=0.0304 | T/A | C/A（= Zoodsma ✓） |
| rs174585 | A/G EAF=0.1644 ｜ T/G EAF=0.0081 | T/G | A/G |
| rs2074550 | C/A EAF=0.2965 ｜ T/A EAF=0.0139 | T/A | C/A |
| rs6722374 | A/T EAF=0.0551 ｜ **G/T EAF=0.5573** | G/T | A/T |
| 其余 19 个同型 | 均有 2 行，被误采的都不是与 Zoodsma 匹配的那一行 | | 均与 Zoodsma 一致 |

注意 rs6722374 这类：被误采的**反而是更高频的那一行**，进一步说明"最后一行覆盖"与
频率无关，纯粹是行序 bug。EAF 最低的备选行（0.0072、0.0081、0.0139、0.0023…）在
数据准备时被 `EAF > 0.01` 滤掉，所以送进 METAL 的文件天然只剩正确那一行 ——
**meta 阶段从来没错过**。

### 12.3 新发现的真实隐患：Karjalainen 有 9,049 个重复 rsid

行感知重标顺带统计了两个 txt 全表的重复 rsid：

| 文件 | 不同 rsid 数 | 多于 1 行的 rsid | 其中等位对完全相同 | 位置相同 |
|---|---|---|---|---|
| Karjalainen_RC.txt | 9,085,410 | **9,049（0.10%）** | 11 | 9,036 |
| Zoodsma_RC.txt | 9,712,126 | **0** | 0 | 0 |

9,036 个是"同位置、不同备选等位"（rs10738609 的 `9_22114495_A_G` 与 `9_22114495_A_T`
即此类），13 个是不同位置共享同一 rsid。

**影响面（实测）**

| 范围 | 有重复行的数量 |
|---|---|
| 19,414 个显著 SNP | **19**（全在 Karjalainen 侧，Zoodsma 为 0） |
| 125 个 lead SNP | **0** |

**下游行为一致，且与本项目自洽**

- METAL 对重复 MarkerName **保留第一次出现的行**。实证：rs112148640 在 Karjalainen 有
  `G/C EAF=0.0585`（第一行）和 `T/C EAF=0.0233`（第二行），TBL 只有一个 rs112148640 且
  `Freq1=0.9456` = freq(C) 的两队列 N 加权（0.9415 / 0.9469），即用了第一行。
  全表 11,452,835 行 = 11,452,835 个不同 MarkerName —— TBL 里**没有任何重复 marker**。
- `RC.R` L129 `df1[!duplicated(df1$rsid),]` 同样是**首行优先**，与 METAL 一致 ✓
  （`duplicated()` 只标记第一次之后的出现）。

**建议**：不动它，结果与已发表一致；如果想彻底干净，按 `rsid + chr + pos + A1/A2`
唯一化后重跑，可预期变化的只有这 9,049 个 rsid 的取值（其中 19 个显著、0 个 lead）。

### 12.4 新发现（对复现至关重要）：METAL 会把等位基因规范化成 `ATCG` 顺序

这是本轮抽查最大的副产品。已发表 `METAANALYSIS2.TBL` 的 `Allele1` **不是**输入文件的
`effect_allele`，而是 METAL 内部按固定碱基序 **A > T > C > G** 挑出的那个，并**同步翻转
beta / Z 的符号**。

**证据（各 20,000 个 marker 的样本，非挑选）**

| 检验 | 结果 |
|---|---|
| 仅 Karjalainen 有的 marker（Weight=136,016）：A1 是否为 `ATCG` 序中较前者 | **20,000/20,000 = 100%** |
| 仅 Zoodsma 有的 marker（Weight=434,646）：同上 | 双等位 SNP **16,589/16,589 = 100%**（3,411 个例外全是 indel，保持输入顺序） |
| 输入顺序被翻转时，TBL 的 Z 符号是否随之翻转 | **9,613/9,613**（另 2 个是 `Z=0.000` 四舍五入的假象） |
| Freq1 是否 == A1 在该队列的频率 | **20,000/20,000 = 100%** |

实例 rs489022：两队列都是 `effect=C, other=T`（Kar beta=+0.0666，Zoo beta=+0.0415，
符号一致、无需任何链校正），但 TBL 输出 `Allele1=T, Allele2=C, Freq1=0.2465`、
`Z=−19.892` —— METAL 把它翻到 T 上，符号也翻了。按 Karjalainen freq(C)=0.7509、
Zoodsma freq(C)=0.754293 反推 freq(T)=0.2491/0.2457，N 加权 = 0.2465 ✓。

**为什么必须是 METAL 内部做的**：`Karjalainen_RC.txt` 自 2025-10-22 起未再改写
（mtime 未变），而 meta 跑于 2025-12-03；该文件里 rs489022 明确是 `effect=C`。
所以规范化不是上游文件预处理的结果。
（**推论**：`metal_rc.txt` 现有的 `ALLELELABELS effect_allele other_allele`
**不需要改**，METAL 自己会做。）

### 12.5 由此解释手稿里费解的 13,761 个 `eaf > 0.99`

`RC.txt` 的 `eaf` 列直接来自 TBL 的 `Freq1`（= A1 的频率）。当 A1 恰好是**高频**那个
等位时，单侧 `EAF > 0.01` 放进来的 MAF<1% SNP 就以 eaf>0.99 的形式出现在 meta 输出
里。实测（`data source/RC.txt`，11,071,938 行）：

| 指标 | 数值 |
|---|---|
| `eaf > 0.99` | **13,761** |
| 其中 indel | 1,910 |
| 其中双等位 SNP（例 rs10000222 T/C eaf=0.9981，freq(C)=0.0019） | 11,851 |
| `eaf < 0.01` | 0 |
| 这 13,761 个里 P<5e-8 的 | 0 |
| 125 个 lead 里 `eaf>0.99` 的 | 0 |

rs10000222 在 Karjalainen 只有一行 `effect=T, other=C, eaf=0.9981`，单侧过滤放行，
METAL 规范化后 A1=T（T 在 `ATCG` 序中先于 C），于是 Freq1=0.9981。
→ **双侧过滤（`eaf > 0.01 & eaf < 0.99`）正好把这 13,761 行全部剔除**，
与 §11.4 的 11,071,938 − 11,058,177 = **−13,761** 逐一对上。

### 12.6 对重跑计划的影响

1. **等位基因方向没有问题**，不需要任何"校正互换"的额外步骤；`metal_rc.txt` 维持
   `ALLELELABELS effect_allele other_allele` 即可复现已发表的 A1/A2/Freq1/Zscore。
2. 重跑后 `RC.txt` 的 `effect_allele`/`other_allele`/`eaf`/beta 符号应与已发表一致
   （规范化 + 符号翻转由 METAL 完成）；变化的只有 p 值列（GC）和被剔除的 13,761 行。
3. RC.R L129 的首行去重与 METAL 行为一致，**保持不动**即可复现；若改为按 rsid+pos+allele
   唯一化，会得到一组"更干净但不同"的数字（0 个 lead 受影响，可安全尝试）。
4. ANNOVAR / FUMA 的 A1/A2 全部来自 `RC.txt`，因此与已发表的等位标注一致，无需重做方向。

---
---

## 13. "把这些都修了，已发表的结果会变多少"——定量影响评估

先回答最要紧的一条：**磁盘上没有任何已发表产物被改动。** 本轮所有审计脚本都是只读的；
`RC.R` 的 5 处修改没有重跑；`METAANALYSIS2.TBL`、`RC.txt`、下游的 FUMA / ANNOVAR / SuSiE /
FLAMES / SMR / coloc 产物全部原样；手稿是另存新文件，原 `RC_GWAS_manuscript.docx` 未动。
所以**已经跑出来的东西照旧有效**。下面量化的是"若按修复版重跑"的影响。

### 13.1 零影响（纯解释性发现，或被下游公式免疫）

| 项目 | 影响 |
|---|---|
| 等位基因方向审计（107 swapped / 24 mismatch） | **0**。131/131 经核对 METAL 处理正确，无需任何校正步骤 |
| METAL 的 `ATCG` 等位规范化（§12.4） | **0**。它本来就是 METAL 行为，重跑会得到同样的 A1/A2/beta 符号/eaf；若手动预规范化反而会把符号翻错 |
| Karjalainen 9,049 个重复 rsid（§12.3） | **0**。METAL 与 `RC.R` L129 都是首行优先，保持一致即可复现 |
| se 公式多乘一个 f（§6 第 5 条） | **0 对检验统计量**：z = β/se 不变 → p 值、clumping 阈值、FUMA/MAGMA 输入、PoPS、FLAMES、fine-mapping 全部不受影响 |
| 手稿 6 处数字改动 | 只改了数据集描述/祖源/meta N 的**表述**，未动任何位点级结果 |

### 13.2 双侧 EAF：删 13,761 行，但 0 个结果变化

被删的 13,761 个 MAF<1% 变体里 **P<5e-8 的有 0 个**，125 个 lead 中 `eaf>0.99` 的 **0 个**
（§11.4、§12.5 实测）。唯一要动的是手稿里的 **11,071,938 → 11,058,177**。

### 13.3 Genomic control：|Z| 全表压缩，133 个 SNP 掉出显著集合、0 个 lead

用两个队列 txt 的 β/se 逐位点重建 meta Z（METAL `SCHEME SAMPLESIZE` 的 √N 加权等价于
`(√N₁z₁ + √N₂z₂)/√(N₁+N₂)`），先把两队列符号对齐（bi-cohort marker 中有 **107 个**需要翻转
——与 §12 查出的 107 个真实互换**完全吻合**，这是方向审计正确性的又一次独立验证），
再按 √λ 压缩（Kar ×1.0718、Zoo ×1.1080），与阈值 |Z|=5.451 比较：

| 结果 | 数量 |
|---|---|
| 19,412 个显著 SNP 中仍显著 | 19,279 |
| **掉出 P<5e-8** | **133（0.7%）** |
| 125 个 lead SNP 中仍显著 | **125 / 125** |
| 掉出显著集合的 SNP 的原 \|Z\| 区间 | 5.45 – 6.04（全部是贴着阈值的那一批） |

复算脚本：`Temp/gc_impact.py`。

**解读**：因为 GC 只把 Z 压掉 7–11%，只有原本落在 `[5.451, 6.04]` 这条窄带里的 SNP 会翻船。
**locus 层面没有变化（125 个 lead 全部保住）**，掉的是每个信号峰边缘最弱的 133 个 SNP。
因此 FUMA 位点定义、fine-mapping credible set、FLAMES/PoPS 的基因排序结论都稳。
要跟着改的是：19,414 → **19,281**（近似值，含双侧 EAF 的 −0 与 GC 的 −133，实际要以重跑为准）、
novelty 统计里与显著 SNP 计数相关的数字。

### 13.4 真正需要"返工"的是三件事（执行记录见 §14、§15）

se 公式 bug 让 β、se 的绝对量纲被放大 `1/√f`（f=0.25 时 2 倍），而 z=β/se 不变。
凡是**只用 z 或 p** 的下游全部免疫；凡是**把 β、se 分开用**的下游都要重做：

| # | 下游 | 为什么受影响 | 实测影响 |
|---|---|---|---|
| 1 | **SMR / HEIDI** | 直接用 b_GWAS、se_GWAS | p_SMR 一位不差；但 **p_HEIDI 被系统性打坏**，Liver 的 SMR+HEIDI 显著基因 **7 → 23**（§14.2.1） |
| 2 | **coloc**（`varbeta = se²`） | Wakefield ABF 里 `r = W/(W+V)` 依赖 V 本身 | PP.H4 变化 −0.0002 ～ −0.024，但 **SNP 级 PP 变化明显**（CSK 头号 SNP 0.329 → 0.224，§15） |
| 3 | **ANNOVAR INFO 的 BETA/SE** | 表里直接印 β、se | 重做后发现磁盘上那份 VCF 还是更早一版 meta 的产物（§14.6） |

### 13.5 汇总表

| 环节 | 是否需要返工 | 说明 |
|---|---|---|
| 数据准备（数据不变） | 否 | 只是加 GC 与双侧过滤 |
| meta | 是（有意为之） | p 值全变；−13,761 行 / −133 个显著 SNP |
| clumping / lead SNP | 否 | 125/125 保住，位点定义不变 |
| novelty 统计 | 可能要 | 与显著 SNP 计数挂钩的数字 |
| FUMA / MAGMA / ANNOVAR 位点注释 | 否 | A1/A2 与坐标不变 |
| SuSiE / fine-mapping | 否（除非用 p 值过滤） | lead 未变 |
| PoPS / FLAMES | 否 | 依赖 lead SNP 与基因重叠 |
| **SMR / HEIDI** | **是** | p_SMR 不变，但 HEIDI 坏了 → 显著基因名单大幅变化（§14） |
| **coloc** | **是** | PP.H4 小幅变，SNP 级 PP 明显变（§15） |
| **ANNOVAR INFO** | **是** | 且旧文件本身就是旧 meta（§14.6） |
| 手稿 | 6 处已改；计数类数字待重跑 | 11,071,938→11,058,177；19,414→19,281 |

---

## 14. 返工执行记录：SMR + ANNOVAR（2026-10-18）

按 §13.4 的结论，se 公式 bug 的唯一实际后果是 SMR 输入与 ANNOVAR INFO 的 BETA/SE。
本节记录这两项的实际执行与验证。

### 14.1 修正后的 SMR GWAS 输入

| 文件 | 行数 | 说明 |
|---|---|---|
| `smr/RC_smr_fixed.txt` | 11,071,938 | 修正后全量输入，列 `SNP A1 A2 freq b se p n` |
| `smr/RC_smr_fixed_besd.txt` | 2,868,162 | 按 49 个 besd 的 SNP 并集子集化（见 14.3） |

变换只有两行：`k = sqrt(freq)`，`se_ok = se*k`、`b_ok = b*k`。
与闭式正确式 `1/sqrt(2f(1-f)(N+Z²))` 逐行对照，**最大相对误差 5.300e-15**（rs79265268），
行数与 SNP 集合与原 `RC_smr.txt` 完全一致。例：

```
1:100004463_TA_T   f=0.2360   旧 b=1.871e-04  se=5.198e-03
                             新 b=9.092e-05  se=2.526e-03   (比值 0.486 = √f)
```

### 14.2 理论上 SMR 的 p 值不变，只有效应量缩放

SMR（Zhu et al. 2016）的检验量：

```
b_SMR  = b_GWAS / b_eQTL
se_SMR = sqrt( se_GWAS² + b_SMR² · se_eQTL² ) / |b_eQTL|      (delta 法)
```

若 `b_GWAS` 与 `se_GWAS` 同乘 `k`，则 `b_SMR → k·b_SMR`，
而 `se_SMR²` 括号里的两项都乘 `k²`，故 `se_SMR → k·se_SMR`，
于是 `T = b_SMR²/se_SMR²` **完全不变 → p_SMR 逐位不变**。

因此：

| 量 | 修正后的变化 |
|---|---|
| `p_GWAS` | 不变（z=β/se 不变） |
| `p_SMR` | **不变** |
| `b_SMR`、`se_SMR` | **乘 √f**（f = topSNP 在 GWAS 里的频率）→ 报告的效应量必须换 |
| `p_HEIDI` | **会变**：HEIDI 跨 SNP 比较 b_SMR，每个 SNP 的 `k_i=√f_i` 不同，不是统一缩放 |

也就是说：**SMR 的 p 值不该变，但表里的 b_SMR/se_SMR 数值全是错的**，
而且 `p_HEIDI` 会变——下面 14.2.1 的实测说明，这一项的变化**远超预期**。

### 14.2.1 A/B 实测：p_SMR 一位不差，但 HEIDI 被 bug 打死了

同一组织（GTEx v8 Liver）、同一 besd、同一 LD 参考、同一 SNP 集合（`nsnp_HEIDI` 与
`topSNP` 2493/2493 完全相同），只把 GWAS 输入换成修正前的 `RC_smr_old_besd.txt`：

| 检查项 | 结果 |
|---|---|
| `p_SMR` | **2493/2493 逐字节相同**（最大相对差 0）✓ 理论完全成立 |
| `b_SMR` 新/旧 vs `sqrt(f_GWAS)` | 最大相对偏差 **7.25e-06**（输出文件的 6 位有效数字造成）✓ |
| `p_HEIDI` | **只有 156/2493 相同** |
| SMR 显著 probe（p_SMR < 2.01e-05） | 旧 46 / 新 46，**集合完全相同** ✓ |
| 其中通过 HEIDI（p_HEIDI>0.05） | 旧 **7** → 新 **23** |
| 反向（旧通过 → 新不通过） | **0** |

最典型的一行：

```
SORT1  rs12740374   p_SMR = 6.19e-41（两次完全相同）
                    p_HEIDI  旧 1.92e-18  →  新 0.195
```

**机理**：HEIDI 检验的是"同一 cis 区间内各 SNP 的 b_SMR 是否一致"。bug 把每个 SNP 的
b_SMR 乘了 `1/√f_i`，而一个位点里各 SNP 的 f 可以差 10 倍以上 → 人为制造出巨大的异质性
→ p_HEIDI 被压到 1e-18 量级 → **HEIDI 判为"异质"，该基因被丢掉**。
换句话说：**se bug 不只是让效应量偏大，它把 HEIDI 过滤器整个打坏了，
使已发表的 SMR 基因列表系统性漏报**。Liver 一个组织就漏了 16 个（23 → 7）。

（注：在全部 2493 个 probe 上，"p_HEIDI>0.05"的比例旧 2083 / 新 2123，差别不大；
破坏集中在**信号最强的那批位点**——因为 b_SMR 越大，`1/√f` 造成的绝对异质性越大。
所以受害的恰好是论文最关心的那些基因。）

**已发表的 GTEx v10 结果有同样的指纹**（磁盘上的 `smr/gtex10/`）：

| 数据集 | probe | SMR 显著 | 过 HEIDI | 通过率 | 显著集里 p_HEIDI 中位数 |
|---|---|---|---|---|---|
| GTEx v10 Liver | 3,649 | 60 | **15** | 25% | 0.0024 |
| GTEx v10 Heart_Left_Ventricle | 5,753 | 42 | **3** | 7% | 0.00034 |
| GTEx v8 Liver（修正后，本次） | 2,493 | 46 | **23** | 50% | — |
| GTEx v8 Liver（修正前，A/B） | 2,493 | 46 | **7** | 15% | — |

注意：**通过率本身不是判据**——20 个组织合并做多重校正后，通过率也是 25% 左右，
说明 HEIDI 在强信号上本来就会拒掉不少真异质性位点。判据是**同一份数据换输入前后的对比**
（15% → 50%）。v10 那两行只能算旁证：它的 p_HEIDI 分布极度偏向 0（中位数 1e-3 量级），
形状和 bug 的指纹一致，但 v10 与 v8、分染色体跑与全基因组跑的 Bonferroni 口径都不同，
不能当作定量证据。

### 14.3 子集化提速（等价性已证）

SMR 只检验 besd 与 GWAS 共有的 SNP，且只在这些 SNP 内算 LD，其余 GWAS 行是纯浪费。

| 步骤 | 数字 |
|---|---|
| 49 个 GTEx v8 besd 的 SNP 并集 | **3,313,991** |
| 单组织 SNP 数 | 74,728（Kidney_Cortex）～ 1,248,735（Thyroid） |
| GWAS 子集化 | 11,071,938 → **2,868,162**（25.9%） |
| LD 参考（g1000_eur 子集） | 22,665,064 → **2,930,100** variant |

**等价性验证**：Liver 全量输入跑 24 min，子集输入跑 12 min，
两份 `.smr` 逐行 diff —— `identical 2493 / full-only 0 / sub-only 0`，**bit-identical**。

### 14.4 GTEx v8 全 49 组织结果（修正后）

**49/49 个组织跑完**（总耗时 114.5 min），共 **224,632** 个 probe–组织检验。多组织 Bonferroni 阈值
`0.05/224632 = 2.23e-07`：**932 个显著，其中 205 个通过 HEIDI，落在 69 个基因上**。

（14.2.1 的 A/B 用的是"单组织 Bonferroni"`0.05/2493 = 2.01e-05`，所以那里 Liver 是 46/23；
换成多组织阈值后 Liver 是 21/10。两个口径都对，只是严格程度不同。）

**Liver（多组织阈值下 10 个基因）**

| Gene | topSNP | b_SMR | p_SMR | p_HEIDI |
|---|---|---|---|---|
| SORT1 | rs12740374 | −0.0471 | 6.19e-41 | 0.195 |
| CELSR2 | rs12740374 | −0.0577 | 3.00e-28 | 0.164 |
| DNAH11 | rs12670798 | 0.0254 | 7.50e-12 | 0.114 |
| RHCE | rs9689 | 0.0289 | 3.12e-10 | 0.093 |
| LIPC | rs2070895 | −0.1648 | 7.14e-10 | 0.199 |
| PKN3 | rs9697210 | −0.0092 | 4.29e-09 | 0.722 |
| ALDH1A2 | rs12908474 | −0.0785 | 9.68e-09 | 0.084 |
| GAS6 | rs6602909 | −0.0346 | 3.60e-08 | 0.159 |
| MYRF | rs198465 | −0.0127 | 8.82e-08 | 0.076 |
| FADS3 | rs174560 | 0.0814 | 1.16e-07 | 0.096 |

**其他心血管/代谢组织**

| 组织 | 基因（p_SMR） |
|---|---|
| Whole_Blood | **PCSK9** 2.8e-11、DOCK7 3.0e-19、PKN3、MYRF、GAS6-AS1、NTN5（*FUT2* 位点 rs601338）、MRPL45P2、RPL36P4 |
| Heart_Atrial_Appendage | **CETP** 7.5e-08、DOCK7 2.9e-08、SDCBP、RPL36P4 |
| Heart_Left_Ventricle | RHCE 1.3e-14、POLK 2.1e-07 |
| Muscle_Skeletal | **HMGCR** 2.7e-09（他汀靶点）、CELSR2 3.6e-52、CNIH4 |
| Artery_Aorta | TMEM57 9.3e-10、POLK、MKRN2 |
| Artery_Tibial | TMEM57 7.6e-10、ANKDD1B |
| Artery_Coronary | GEMIN7 1.1e-07、CEACAM19 1.5e-07 |
| Adipose_Visceral | CEACAM19 1.4e-09、RHCE、ULK3、ANKDD1B |
| Pancreas | RHCE 1.8e-14、CELSR2 1.8e-12、SPTY2D1、ULK3、DNAJC13 |
| Small_Intestine | DOCK7 7.2e-08 |

**生物学上说得通**：PCSK9（Whole_Blood）、HMGCR（骨骼肌，rs3846662 正是那个经典剪接位点）、
CETP（心房，rs1800775 启动子位点）、LIPC/SORT1/CELSR2（肝脏）、FADS3——全是降脂药靶或已知脂质基因。
1p13.3 的 SORT1/CELSR2/PSRC1 三个基因，在 Liver 里 **SORT1 与 CELSR2 过 HEIDI（0.195 / 0.164），
PSRC1 不过（0.002）**；PSRC1 在绝大多数组织的 p_HEIDI 都 <0.01。这与"PSRC1 只是邻近基因、
真正效应基因是 SORT1"的既有认识一致。

**一条要注意的假阳性**：`RPL36P4`（核糖体蛋白 L36 假基因 4）在 **40/49 个组织**里显著且过 HEIDI，
topSNP 在多个组织里都是 rs7261820 / rs12481365。核糖体蛋白假基因是典型的重复序列/多重比对
artifact，**不应作为 effector gene 报告**。同理 `SDHDP6`（琥珀酸脱氢酶假基因）、
`MRPL45P2`（线粒体核糖体蛋白假基因）也要打问号。

完整结果：`smr/gtex8_smr_all.tsv`（224,632 行，含所有列）、`smr/gtex8_smr_summary.tsv`
（每组织计数）、`smr/gtex8_smr_top.tsv`（每组织 top5）。汇总脚本 `Temp/summarize_smr_gtex8.py`。

### 14.5 ANNOVAR：为什么只改 INFO 就等价于重跑

`table_annovar.pl` 把 VCF 的 INFO 字符串**逐字**拷进 `Otherinfo1` 列，
而 `Func / Gene / ExonicFunc / AAChange` 只依赖 `CHROM/POS/REF/ALT`。
所以**只替换 INFO 列，与真正重跑 ANNOVAR 的结果完全相同**——已在 109 行上逐列 diff 验证：
除 INFO 外所有注释列逐字节相同。

执行时的两个坑：

1. multianno 输出的表头把所有 VCF 列都叫 `OtherinfoN`，**不能按列名定位**。
   必须按内容定位：含 `AF=...BETA=` 的列是 INFO（本例第 21 列 / `Otherinfo11`），
   形如 `rs\d+` 的列是 ID（第 16 列 / `Otherinfo6`）。
2. `REF/ALT` 与 `effect/other` 的对应关系必须先核对：**109/109 一致**（ALT = effect_allele），
   否则 BETA 会整体反号。

新 INFO 格式串：`AF=%s;BETA=%.6g;SE=%.6g;P=%s;Z=%.5g;N=%d`。
例（rs10162642）：

```
旧  AF=0.2024;BETA=-0.076139;SE=0.007805;P=1.76e-22;Z=-9.755;N=251094
新  AF=0.2073;BETA=-0.023895;SE=0.002309;P=4.208e-25;Z=-10.349;N=570662
```

### 14.6 重做 ANNOVAR 时发现的更大问题：磁盘上这份 VCF 是**旧 meta** 的产物

这一条超出了 se 公式 bug 的范围，但必须记下来。

| | 磁盘上的 `RC_clump_annovar_final.vcf` | 当前 `data source/RC.txt` |
|---|---|---|
| 文件时间 | 2025-10-25 | 2025-12-16 |
| N | 251,094（108 行）/ 136,016（1 行） | 570,662（107 行）/ 136,016（2 行） |

- 109/109 行的 INFO **全部**与 `RC.txt` 不符（AF 变 107 行、BETA/SE 109 行、P/Z/N 各 108 行）。
- `Z_new / Z_old` 的比值在 **0.632 – 1.567** 之间（中位 1.045），离散度远超 GC 的 7–11%
  → 这不是"统一缩放"，而是**两次不同的 meta**。`251,094` 正是手稿里那个错误 N（§7.1）。
- 用最终 meta 重建 INFO 后：109 个位点中 **95 个仍 P<5e-8，14 个掉出，0 个新增**。

| rsid | 旧 P | 新 P |
|---|---|---|
| rs632057 | 1.492e-09 | 1.33e-04 |
| rs799157 | 5.147e-10 | 1.93e-07 |
| rs1590381 | 4.659e-10 | 4.44e-07 |
| rs13300853 | 3.002e-09 | 1.47e-07 |
| rs2153669 | 3.264e-09 | 5.90e-07 |
| rs11057841 | 6.672e-09 | 2.28e-05 |
| rs13389219 | 1.353e-08 | 2.66e-07 |
| rs140999088 | 1.268e-08 | 3.40e-06 |
| rs11057604 | 2.381e-08 | 7.76e-06 |
| rs13287546 | 2.533e-08 | 3.25e-06 |
| rs807060 | 2.802e-08 | 3.51e-05 |
| rs1134030 | 3.098e-08 | 7.62e-08 |
| rs28929474 | 2.239e-08 | 1.80e-05 |
| rs11118618 | 5.660e-09 | 5.18e-08 |

**结论**：ANNOVAR 表不是"改一下 se"的问题，而是"整表要按最终 meta 重建"的问题。
如果论文里的变异注释表是从这份 2025-10-25 的 VCF 生成的，需要核对上面这 14 个位点
（其中 rs28929474 是 SERPINA1 的 PiZ 位点，rs632057 在 6q23.3，都是有生物学含义的位点）。
反之，如果论文表格本来就是从 `RC.txt` 重新导出的，那这份 VCF 只是中间文件，可弃用。

### 14.7 交付物

| 文件 | 内容 |
|---|---|
| `smr/RC_smr_fixed.txt` | 修正后 SMR 输入（全量，783 MB） |
| `smr/RC_smr_fixed_besd.txt` | 子集版（201 MB），实际跑批用 |
| `smr/ld_ref/gtex8.{bed,bim,fam}` | LD 参考子集（2,930,100 variant） |
| `smr/gtex8_out_sub/*.smr` | 49 组织 SMR+HEIDI 结果（224,632 个 probe–组织检验） |
| `smr/gtex8_smr_all.tsv` / `_summary.tsv` / `_top.tsv` | 合并表 / 每组织计数 / 每组织 top5 |
| `smr/run_smr_gtex8.sh` | Ubuntu 版一键脚本（`bash run_smr_gtex8.sh [tissue...]`） |
| `smr/gtex8_ab/Liver_old.smr` | A/B 对照：同一组织、修正前输入 |
| `smr_tool/LV_coloc_results_fixed.csv` | 修正后的三个 LV 位点 coloc 结果 |
| `RC_clump_annovar_final_fixed.vcf` | 按最终 meta 重建 INFO 的 VCF（109 行） |
| `my_first_anno.hg19_multianno_fixed.txt` | 对应 multianno（除 INFO 外逐字节不变） |
| `Temp/build_smr_input_fixed.py`、`Temp/subset_smr_gwas.py`、`Temp/fix_annovar_betase.py` | 三个重做脚本 |
| `Temp/ab_check_smr.py` | SMR 修正前/后的 A/B 对照脚本 |
| `smr_tool/run_coloc_LV_fixed.R`、`Temp/run_coloc_ab.R` | coloc 修正版脚本 + A/B 对照脚本 |

两个 ANNOVAR 文件已拷到 Ubuntu 的 `~/annovar/`（本项目里是 WSL 的 `\\wsl.localhost\Ubuntu\...`）。
Windows 侧没有 Perl（Rtools43 只剩环境变量残留），所以没有在本地真跑 `table_annovar.pl`；
如需二次确认，在 Ubuntu 上执行：

```bash
cd ~/annovar
table_annovar.pl RC_clump_annovar_final_fixed.vcf humandb/ \
    -buildver hg19 -out my_first_anno_fixed \
    -protocol refGeneWithVer -operation g -nastring .
```

---

## 15. 顺带查出来的第三处返工：coloc

§13 原先把 se bug 的影响限定在 "SMR + ANNOVAR"。做 §14 时回头核对下游，发现 **coloc 也中招**，
而且它同样是"β 和 se 分开用"的类型。

### 15.1 哪些 coloc 用了 buggy 的 se

| 脚本 | 代码 | 用到的 se |
|---|---|---|
| `smr_tool/run_coloc_LV.R` L30-33 | `d2 <- list(beta=d$beta, varbeta=d$se^2, ...)` | `RC.txt` 的 `se` ✗ |
| `RC.R` L1527-1530（INHBB、HGFAC 等 pQTL 位点，共 8 段） | `data$VAR.x = data$se^2` → `data1` 的 `varbeta` | `RC.txt` 的 `se` ✗ |

即：**所有把 `RC.txt$se` 当 varbeta 的 coloc 都要重跑**。eQTL/pQTL 一侧的 se 是原始数据，不受影响。

### 15.2 为什么 z 不变但 ABF 会变

coloc 的 Wakefield 近似：

```
lABF = 0.5 * ( log(1-r) + z² · r ),   r = W / (W + V),   V = varbeta = se²
```

`z = β/se` 确实不变，但 `r` 直接依赖 `V`。bug 让 `V → V/f`，
低频位点的 `V` 被放大得最厉害（f=0.01 时放大 100 倍），`r` 显著变小 → ABF 被压低。
所以**低频 SNP 在旧结果里被系统性地低估了"因果"资格**。

### 15.3 实测（三个 LV 位点，用真 coloc 5.2.3 跑两遍）

`coloc` 5.2.3 装在 R 4.3 的库里，R 4.5.1 直接 `.libPaths()` 指过去可以正常加载。
先验证用旧 se 能逐位复现已发表的 `LV_coloc_results.csv`（PP.H4 0.971502 / 0.510643 / 0.940059，
完全一致），再换成修正后的 se：

| 基因 | nSNP | PP.H4（旧 → 新） | PP.H3（旧 → 新） | ΔPP.H4 |
|---|---|---|---|---|
| CWF19L1 | 5,014 | 0.9715 → 0.9713 | 0.0278 → 0.0279 | −0.0002 |
| **CSK** | 4,069 | **0.5106 → 0.4865** | 0.4894 → 0.5135 | **−0.0241** |
| THOP1 | 7,616 | 0.9401 → 0.9343 | 0.0594 → 0.0651 | −0.0058 |

**结论层面基本稳**：三个位点里两个仍是很强的 H4（>0.93），只有 **CSK 从 0.511 掉到 0.487**，
跨过了 0.5 这条线——从"略偏共享因果"变成"略偏两个独立信号"。如果论文里对 CSK 的表述是
"coloc 支持共享因果（PP.H4=0.51）"，需要改成"证据两可"。

**但 SNP 级的结果变化明显**：

| 基因 | 头号 SNP（旧 → 新） | SNP.PP.H4 | \|ΔlogABF\| 中位 / 最大 | 位移 >1 的比例 |
|---|---|---|---|---|
| CWF19L1 | rs17668357 → rs11597086 | 0.101 → 0.100 | 0.268 / 1.893 | 1,155 / 5,014（23%） |
| CSK | rs34933034（不变） | **0.329 → 0.224** | 0.293 / 1.897 | 1,034 / 4,069（25%） |
| THOP1 | rs2260414（不变） | 0.277 → 0.253 | 0.353 / 1.916 | 2,026 / 7,616（27%） |

即：**"哪个变异是共享因果"这个问题的答案没变，但"有多确信"普遍被旧结果高估了**
（CSK 的 0.329 是高估，真实约 0.22）。四分之一左右的 SNP 的 logABF 位移超过 1 个单位，
这些几乎全是 MAF<5% 的低频位点。

### 15.4 交付

- `smr_tool/run_coloc_LV_fixed.R`：`run_coloc_LV.R` 的修正版（只多了 `k=sqrt(eaf)` 两行缩放），
  跑完写 `smr_tool/LV_coloc_results_fixed.csv`。
- `Temp/run_coloc_ab.R`：A/B 对照脚本（同一份数据跑两遍，输出上面的表）。
- `RC.R` 里那 8 段 pQTL coloc 没有一并重跑——它们各自要读 `3617_80_*.txt.gz` 等外部 pQTL 文件，
  且位点不同。**重跑方法完全一样**：在 `data$VAR.x = data$se^2` 之前插一行
  `RC$se <- RC$se * sqrt(RC$eaf); RC$beta <- RC$beta * sqrt(RC$eaf)` 即可。
