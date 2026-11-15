# 第二轮内部审稿意见的修改说明

分支：`review-round2-fixes`，基于 `f7bcce0`。
对应的审稿意见：`review_lrd-width-luminosity.md`（第一轮）、`review_round2_lrd-width-luminosity.md`（第二轮）。下文用 M1–M5、N1–N7 等编号指代意见条目。

本文逐条列出改了什么，并给出改前与改后的对照。引用的正文已把宏展开成数值，便于阅读。行号指本分支的 `paper/ms.tex`。

---

## 0. 概要

**结论层面的改动只有一处：头条上限的表述。** 基线上限 0.075 dex 保留，改称“基线上限”（baseline limit）。另立一个“保守上限”0.10 dex，覆盖窄线分解、LSF 和样本定义的各种合理选择。0.10 是这些变体里最大的单侧 95% 上限，来自 [O III] 绑定、并按 §2.3 规则重新做质量筛选的那个变体。第 5 节的所有物理约束仍以 0.075 为准，同时注明在 0.10 下放宽多少。重复光谱给出的 0.06 dex 不再出现在摘要中，改为“取决于误差模型，0.049–0.062 dex”。

**需要两位作者决定的事项**
1. **标题**：已改为 *Nearly Uniform Broad Hα Widths in Little Red Dots at Fixed Luminosity and Their Consequence for Virial Black Hole Masses*，理由见 M4。如不同意，改回原标题只需改 `ms.tex` 第 25、29 行。
2. **期刊与篇幅**：用同一方法统计，正文（不含摘要、图注、表格、附录、参考文献）从 **3880 词增至 4617 词**，摘要从 205 词增至 241 词，PDF 从 13 页增至 15 页。投 ApJ 没有问题。如果投 ApJL，需要再删约 1100 词，建议把第 4 节的变体段落和 A2744-QSO1 段落移入附录（N4）。
3. **保守上限的取法**：现在取的是“有逐源输入、能算单侧 95% 上限的变体中的最大值”。±8000 km/s 拟合窗口只存了分位数，单侧上限算不出来；它的 95% 区间上端是 0.119 dex，正文已单独写明（第 227、232 行）。如果希望保守上限把这个变体也覆盖进去，可以改用 0.12，第 5 节的放宽比例相应变为 60%。

**验证**
- `make analysis`（14 个脚本）在仓库副本上完整重跑。除 `scattering_kernel.json`（及由它派生的 `review_checks.json` 中两个数）有 10⁻⁸ 量级的浮点差外，所有结果文件与提交的版本逐字节一致。
- `make tables numbers figures` 和 `scripts/check_reproduction.py` 全部通过。`paper_reference_hashes.json` 已更新为新的 `numbers.tex` 与 `numbers.json` 的哈希。
- `analysis/sample/make_selection_audit.py` 生成的 `selection_audit.json` 与仓库文件逐字节一致。
- 论文编译无未定义引用，也没有卡住的浮动体。双盲版本中已不含 Zenodo DOI，构建脚本里加了断言检查。
- 所有新增数值都由脚本从仓库内存储的逐源测量值算出，没有手写数值。

**无法完成的事项**：凡是需要原始光谱重新拟合的，本次都没做，见第 9 节。

---

## 1. 标题与摘要

### 标题（M4）
- 改前：A Tight Broad-line Width–Luminosity Relation in Little Red Dots and Its Consequence for Virial Black Hole Masses
- 改后：Nearly Uniform Broad Hα Widths in Little Red Dots at Fixed Luminosity and Their Consequence for Virial Black Hole Masses
- 理由：斜率只由 3 个最亮源决定，而且依赖窄线分解（b = 0.09–0.17）。作者自己也写了 “None of our conclusions relies on its value”。所以标题不宜突出“关系”和斜率。

### 摘要第 3 句：头条上限（M1、M2、N1、N7）——第 43 行
- 改前：At fixed line luminosity, LRD broad Hα widths have an intrinsic scatter below 0.075 dex (95%), and below 0.06 dex and consistent with zero once the variance between repeat spectra is included.
- 改后：At fixed line luminosity, LRD broad Hα widths have an intrinsic scatter of 0.054 dex, below 0.075 dex (95%) in our baseline fits and below 0.10 dex under alternative narrow-line decompositions, line-spread functions, and sample definitions; repeat spectra attribute part or all of it to measurement variance.
- 理由：
  - 0.06 dex 是手写的取整值。它依赖误差模型：加性地板下上限为 0.049，乘性模型下为 0.062。
  - “consistent with zero” 只在加性模型下成立。乘性模型给出 0.041（0.025–0.067）。

### 摘要新增：常数线宽（M4）——第 44 行
- 改后：Below the three most luminous LRDs, the widths are consistent with a constant ≃1760 km s⁻¹.

### 摘要第 4–5 句：对照组与误差放大（M5、次要 9）——第 45 行
- 改前：Other JWST broad-line sources scatter by 0.15 dex and identically fitted SDSS quasars by 0.12 dex. Five uniform width definitions and independent fits to the same public spectra give the same tightness.
- 改后：Other JWST broad-line sources scatter by 0.15 dex and identically fitted SDSS quasars by 0.12 dex, a contrast that persists when all measurement errors are inflated by up to a factor of 2.8.
- 理由：
  - 误差放大测试是对照组结论最有力的支撑，第一轮建议写入。
  - “give the same tightness” 删去：非参数 FWHM 为 0.076，比基线大约 40%。正文与 Summary 中改为 “similar”。

### 摘要维里质量句（N6）——第 46 行
- 改前：…scatter intrinsically by less than 0.16 dex at fixed Hα luminosity…
- 改后：…scatter intrinsically by less than 0.16 dex (0.21 dex conservatively) at fixed Hα luminosity…

### 摘要光球句（M3）——第 47 行
- 改前：…yet at fixed Hα luminosity the widths do not follow the photospheric scaling.
- 改后：…while at fixed Hα luminosity we detect no dependence on photospheric luminosity in 13–15 sources.
- 理由：这个检验约为 2σ，而且去掉 RUBIES-UDS-29813 后不再排除光球标度（见第 3 节）。

### 摘要 λ 范围（次要 3）——第 48 行
- 改前：λ ≤ 0.16–0.30
- 改后：λ ≤ 0.10–0.30
- 理由：与 §5.3、Summary 和表 3 统一（`\lamVirLo`–`\lamIsoHi`）。

---

## 2. 第 1–3 节

| 位置 | 意见 | 改前 | 改后 |
|---|---|---|---|
| 第 60 行 | 文献：Juodžbalis | the dynamical mass of A2744-QSO1 **matches** its virial mass | …**is consistent with** its virial mass（与原文摘要措辞一致） |
| 第 93 行 | 第一轮代码审计 4（§2.2 与表 B1 矛盾）；次要 9 | does not compress parent scatter on winged profiles…, and gives **the same scatter as** every other… | does not compress parent scatter on winged profiles **at the observed widths**…, and gives **a scatter similar to** every other… |
| 第 88 行 | 代码问题 (b) | …from **80** pixel-noise bootstrap realizations… | …from **200**…（实际使用的协方差来自 200 次，`covariance_scope` 字段 44/44 一致） |
| 第 106 行（新增） | 次要 12 | — | CEERS-746 and Skyfire-299 lie 35″ apart at the same redshift (z = 5.623) in separate slitlets of one MSA mask; we treat them as distinct objects. |
| 第 115 行 | 文献：BX26 去红化 | …or by **class-median** reddening otherwise | …or by **their sample-median** reddening otherwise |
| 第 159 行（新增） | M4 | — | The remaining 26 LRDs, spanning log L_Hα = 41.8–42.9, are consistent with a constant width: fixing b = 0 gives σ_int = 0.035 dex (0.015–0.067) around a median FWHM of 1760 km s⁻¹. |
| 第 162 行 | N2 | A2744-QSO1 (open diamond…) **lies outside our parent catalogs** and is discussed in Section 4. | A2744-QSO1 (diamonds…) **is not in our sample** and is discussed in Section 4. |
| 第 182–183 行（新增） | M5 | — | Inflating every formal error by factors of 1.5–2.8, or adding 0.044–0.08 dex in quadrature, keeps Δσ_int at −0.12 to −0.08 dex with 95% intervals excluding zero in both samples… The comparison class is heterogeneous, however, and may include outflow-dominated or low-S/N broad components. |
| 第 190–191 行（新增） | M5（SDSS） | — | Removing the χ² rule adds 31 quasars and gives 0.126 dex. The ±5000 km s⁻¹ window spans only about ±1.5 FWHM for these quasars (those rejected for low S/N have median FWHM 7090 km s⁻¹), so truncated wings, if anything, lower the SDSS scatter and width offset. |

---

## 3. 第 4 节：误差标定与筛选检验

### 重复光谱（M2、N7）——第 200–205 行
- 改前：
  > A joint random-effects fit to all 43 passing LRD spectra … gives σ_int < 0.049 dex (< 0.057 dex with RUBIES-EGS-42232 reduced to its primary spectrum). The exact floor applicable to the 29 primary spectra is uncertain: … indicating smaller unmodeled variance in high-S/N primary spectra than in secondary repeat exposures. We therefore treat σ_int < 0.075 dex (formal errors) as a conservative ceiling and σ_int < 0.06 dex (consistent with zero) as the repeat-calibrated bound.
- 改后：
  > A joint random-effects fit …, with an additive per-spectrum floor (posterior median 0.053 dex), likewise gives σ_int < 0.049 dex (< 0.057 dex …). Writing the unmodeled variance instead as a multiple of C_yy (median factor 2.0), so that high-S/N spectra receive less of it, gives σ_int < 0.062 dex. The repeats also reveal a grating offset: in all seven LRDs with passing G395H and G395M spectra the G395H width is narrower (median −0.06 dex; a free offset in the random-effects model gives 0.054 ± 0.022 dex without changing the limit), and the highest-S/N spectrum we retain is on average −0.056 dex narrower than its repeats, as expected because narrower lines reach higher S/N at fixed flux. The 22 LRDs with G395M primary spectra give σ_int = 0.048 dex (< 0.072 dex). The floor applicable to the 29 primary spectra is uncertain: … The repeat spectra therefore limit σ_int to 0.049–0.062 dex depending on the error model, consistent with zero only for an additive floor, and we use formal errors below.
- 说明：
  - 乘性模型、自由 G395H 偏移模型都已加入 `analysis/13_parent_and_repeat_checks.py`。
  - 第二轮记作“7 个源中 6 个”的 G395H 更窄，按每个源取中位后实际是 **7 个全部更窄**（RUBIES-EGS-29489 只差 −0.008 dex）。
  - 删除了 “conservative ceiling” 的自相矛盾表述（N6）。

### 注入测试的位置（第二轮“做得好”第 1 条）——第 210 行
- 改前：49 模板的结果放在 “Including the rejected spectra changes little.” 段落之后。
- 改后：移到模拟段落中，紧接 29 模板的结果，作为主要论据：“The quality cuts are also insensitive to width: on all 49 in-band parent templates, injections shifted +0.10 dex above the locus pass 0.99–1.04 times as often…, and σ_parent = 0.10 dex yields σ_int ≤ 0.054 dex in none of 100 Gaussian and two of 100 exponential trials after the full quality cuts.”

### “in passband” 的定义（N5）——第 212 行
- 改前：49 **lie inside the nominal filter passband**
- 改后：49 **have Hα longward of the order-blocking filter cut-on**
- 说明：判据只检查长通滤光片的截止波长，不检查光栅的名义范围。

### 9 个 S/N ≥ 3 失败源的局限（第二轮“做得好”第 2 条）——第 213 行
- 改后加入：although seven of the nine have ΔBIC < 6 against the narrow-only model

### 窄线分解（N1）——第 216–219 行
- 改前：
  > In two of the three luminous LRDs … the narrow Gaussian reaches its 600 km s⁻¹ FWHM bound… Excluding these two objects (N = 27), fixing every narrow-line FWHM at 150 km s⁻¹, or anchoring it to [O III] λ5007 lowers the slope to 0.09–0.11 and gives σ_int = 0.040–0.072 dex, with 95% ceilings of 0.062–0.097 dex, comparable to halving the LSF…
- 改后：
  > With its FWHM free, the narrow Gaussian sits at the 30 km s⁻¹ lower bound in 17 of the 29 LRDs, and in two of the three luminous LRDs … it reaches the 600 km s⁻¹ upper bound…; [O III] λ5007 is unresolved (at the 30 km s⁻¹ bound) in 21 of the 23 LRDs where it is measured. Anchoring the narrow Hα FWHM to [O III], the standard choice, changes individual widths by 0.06 dex rms (up to 0.21 dex) and gives σ_int = 0.072 dex (< 0.097 dex), or 0.075 dex (< 0.102 dex; N = 27) after re-applying the absorption-bound rule of Section 2.3. Fixing every narrow FWHM at 30 or 150 km s⁻¹ gives 0.070 and 0.046 dex, and excluding the two bound objects gives 0.040 dex. These variants lower the slope to 0.09–0.11, and Δσ_int excludes zero in every case; they hold membership fixed, although whether a broad component is required also depends on the split.
- 说明：
  - 不再把 150 km/s 这一任意选择与物理上最有依据的 [O III] 绑定并列成 “0.040–0.072”。
  - 质量筛选的重新应用：[O III] 重拟合中 RUBIES-EGS-29489 和 RUBIES-EGS-55604 的 τ_abs 撞到边界，按 §2.3 应剔除。
    - 仓库没有存变体的协方差，所以剔除带来的变化用基线协方差计算，再加到存储的精确结果上。
    - 这种近似对全部 29 个源的复现误差小于 0.001 dex（0.0712 对存储值 0.0719），表 2 注释 e 已说明。

### A2744-QSO1（N2）——第 221–224 行
- 改前：
  > A2744-QSO1 … **lies outside our parent catalogs: it is lensed and in the Abell 2744 field.** Fitted with the same pipeline, its free narrow Gaussian also reaches the 600 km s⁻¹ bound … for magnification μ = 7.5 ± 1.5 (Furtak2025). … its single-Gaussian FWHM is 1223 km s⁻¹ (S/N = 20.9, −0.08 ± 0.03 dex from Equation (3)), and adding it to the [O III]-anchored fit leaves σ_int unchanged…
- 改后：
  > A2744-QSO1 … **is not in our sample: it lies in the Abell 2744 field, which none of our three classification catalogs covers, and with a free narrow line its broad component has S/N = 3.4, below our cut.** In that fit the narrow Gaussian reaches the 600 km s⁻¹ bound … for the magnification μ = 7.5 ± 1.5 **of image B** … its single-Gaussian FWHM is 1223 km s⁻¹ (S/N = 20.9, −0.08 ± 0.03 dex from Equation (3)), **but this fit is not preferred over a single 600 km s⁻¹ narrow component without a broad line (ΔBIC = −2.6), so its single-Gaussian width is poorly defined.** Adding it to the [O III]-anchored fit leaves σ_int unchanged…
- 说明：
  - de Graaff+2026 是母星表之一，它的样本包含 QSO1，原来的说法因此不对。
  - 真实的落选原因有两个，两者都可核实：
    - 三个分类星表（Taylor、Kocevski、Hviding）只覆盖 EGS 和 UDS，§2.3 规定缺少分类的源要排除。
    - 窄线自由时 S/N = 3.4。
  - Furtak 的放大率注明了是 image B。

### 其他变体（M1、次要 11）——第 226–230 行
- 改前：
  > Halving the nominal LSF width gives σ_int = 0.069 dex (95% limit 0.094 dex), widening the fit window to ±6800 km s⁻¹ gives 0.062 dex, and admitting every source labeled an LRD in any of the four parent catalogs gives 0.059 dex across 35 objects.
- 改后：
  > Halving the nominal LSF width gives σ_int = 0.069 dex (< 0.094 dex), and widening the fit window to ±6800 km s⁻¹ or to ±8000 km s⁻¹ with [S II] masked gives 0.062 and 0.083 dex (0.059–0.119). Retaining the lowest-S/N passing spectrum of each object gives 0.059 dex (< 0.083 dex), and admitting every source labeled an LRD in any of the four parent catalogs gives 0.059 dex across 35 objects (0.069 dex across 42 with all rules relaxed, mostly from one source with reduced χ² = 66). The primary spectrum of RUBIES-EGS-42232 (G235M/F170LP) places Hα at 3.91 μm, beyond the nominal G235M range; replacing it by its best G395M spectrum gives 0.052 dex (< 0.072 dex).
- 说明：第一轮 M1 列出的三组未报告结果（±8000 窗口、最低 S/N 光谱、N = 42）现已全部报告。42232 的替换结果是新算的。

### 第 4 节结论段（新增，M1、N1、N6）——第 232–233 行
> Across these fit, narrow-line, and sample variants (Figure 2), one-sided 95% limits reach 0.10 dex, set by the [O III]-anchored fit with the absorption rule re-applied, and only the widest fit window extends the 95% interval further, to 0.119 dex. We therefore quote σ_int < 0.075 dex as the baseline limit and 0.10 dex as a conservative ceiling that covers these analysis choices.

---

## 4. 第 5 节：物理约束

| 位置 | 意见 | 改前 | 改后 |
|---|---|---|---|
| 第 238 行 | N1、N6 | …evaluated at the **conservative** formal-error 95% limit σ_int = 0.075 dex; … the [O III]-anchored ceiling (0.097 dex) loosens them by 30%. | …evaluated at the **baseline** 95% limit σ_int = 0.075 dex; bounds linear in σ_int loosen by **35% at the conservative 0.10 dex ceiling** … and tighten by a quarter at the same-grating repeat-spectrum limit (0.057 dex). |
| 第 247 行 | N1 | …scatters by at most 0.15 dex. | …scatters by at most 0.15 dex (**0.20 dex at the conservative ceiling**). |
| 第 248 行 | 次要 7（Bentz） | …scatters by 0.13–0.19 dex at fixed 5100 Å luminosity, **which alone would contribute** 0.065–0.095 dex to σ_int, comparable to or above our 95% ceiling… | …at fixed 5100 Å luminosity, **including measurement error**; **if a scatter this large held at fixed L_Hα with M independent of R, it would contribute up to** 0.065–0.095 dex to σ_int, comparable to our 95% limit… |
| 第 258 行 | N1 | …2.8 times the entire 0.075 dex ceiling… | …2.8 times the entire 0.075 dex limit (**2.1 times the conservative ceiling**)… |
| 第 259 行 | N1 | (leaving σ_M/R ≤ 0.13 dex at H/R = 1) | (…; **H/R = 0.5 suffices at the conservative ceiling**) |
| 第 290 行 | M3 | the photospheric parameters **do not govern** source-to-source width variations | …**show no detected influence on**… |
| 第 293 行（新增） | M3 | — | This null test rests largely on RUBIES-UDS-29813, which lies 0.5 dex below the L_ph–L_Hα trend; without it the partial slopes become 0.02 (−0.16 to 0.22) and −0.04 (−0.33 to 0.24), consistent with both zero and the simulated photospheric medians (0.14 and 0.21). |
| 第 294 行 | M3 | …to 0.031–0.057 dex, **twice as tightly** as v_K(R_ph)… | …to 0.031–0.057 dex, **1.3–1.4 times more tightly** than v_K(R_ph)…（按 dex 比值；原文的 “twice” 只在按方差算时近似成立） |
| 第 303 行 | 次要 4 | Gas bound … moves no faster than **√2 v_K at R ≥ R_ph** | …moves no faster than **the escape speed, √2 v_K λ^{−1/2} at R_ph and less beyond it** |
| 第 305 行 | 次要 1 | …gives **the identical range** 0.16–0.17 | …gives 0.16–0.17 |
| 第 310 行 | 文献：Curtis | …at most 0.29 (λ < 1.41) at **the H⁻ opacity-valley floor** κ = 10⁻⁴ cm² g⁻¹ | …at most **0.30 (λ < 1.42)** for κ ≈ 10⁻⁴ cm² g⁻¹, **the surface opacity outside the H⁻ band in the atmosphere models of Curtis+2026 (Appendix A)** |
| 第 318 行 | 次要 4 | broadened by at least **(2.5–3.0)/k λ^{1/2}** | broadened by at least **λ^{1/2}(2.5–3.0)/k**: …**at λ = 1** |
| 第 319 行 | 文献：Inayoshi | Electron scattering at τ_e ≃ 2–3 … **produces this factor** | …**can maintain lines this broad** |
| 第 325–326 行 | N3、次要 5 | …scattered widths scale as L^{0.09–0.39} across τ_e = 0.5–3, **overlapping our slopes** (b = 0.09–0.17…) and bracketing that of BX26 (0.25). | …scale as **L^{0.14–0.39} at the τ_e ≃ 2–3 required above**, and only as L^{0.09} at τ_e = 0.5, which supplies too little broadening. This matches our baseline slope (0.17) and that of BX26 (0.25), **but the slopes after correcting the narrow-line split (0.09–0.11, with 95% upper ends 0.13–0.17) lie at or below the lower edge of this range.** |
| 第 332 行 | 次要 8、N6 | …less than 0.16 dex (0.11 dex at the formal-error median, and < 0.12 dex with the same-grating repeat floor). | …less than 0.16 dex (**2.06 σ_int rounded up**; 0.11 dex…, < 0.12 dex…, **and < 0.21 dex at the conservative ceiling**). |
| 第 333 行 | 次要 6 | The virial masses tabulated by BX26 **confirm this directly**: … 0.17 dex for LRDs versus 0.42 dex for comparison sources. | …**are consistent with this**: … 0.42 dex for comparison sources, **although the latter value also reflects heterogeneous mass estimators.** |

δ 的精确值（N6）：`make_numbers.py` 改用含 (1 − R_ph/r_out) 因子的精确式（`*_finite_rout100` 字段）。κ = 10⁻⁴ 时 δ 由 0.29 变为 0.30，λ 由 1.41 变为 1.42；κ = 10⁻³ 时 δ 仍为 0.003–0.030。

---

## 5. Summary（第 339–342 行）

- 第 1 条：
  - 改前：
    > …(one-sided 95% upper limit 0.075 dex); calibrating unmodeled variance from duplicate LRD spectra lowers the 95% limit below 0.06 dex, consistent with zero scatter. Five uniform width definitions … yield the same tightness…
  - 改后：
    > …(one-sided 95% upper limit 0.075 dex, and 0.10 dex across alternative narrow-line decompositions, line-spread functions, and samples); calibrating unmodeled variance from duplicate LRD spectra gives 95% limits of 0.049–0.062 dex depending on the error model. Below the three most luminous LRDs the widths are consistent with a constant 1760 km s⁻¹. Five uniform width definitions … yield similar tightness…
- 第 2 条（M3）：
  - 改前：
    > …yet at fixed L_Hα the widths show no dependence on L_ph (partial slope −0.01, −0.16 to 0.15).
  - 改后：
    > …yet at fixed L_Hα we detect no dependence on L_ph (partial slope −0.01, −0.16 to 0.15), a null result that becomes consistent with photospheric scaling when its most leveraged source is removed.
- 第 4 条（N6）：括号中加入 “< 0.21 dex at the conservative ceiling”。

---

## 6. 附录

- **附录 B（第一轮代码审计 4，第二轮“做得好”第 1 条），第 394 行新增：**
  > Two simplifications limit these suites: injected lines carry no absorption, although five LRDs (including the three most luminous) are fitted with absorption, and the injected noise equals the formal pixel errors, so the reduced-χ² rule never rejects a trial and the unmodeled variance seen in repeat spectra (Section 4) is not represented.
- **附录 C（N3），第 404 行新增：**
  > For τ_e rising as L_Hα^{1/4} to L_Hα^{1/2} (Section 5.3), the predicted width–luminosity slope is D/4 to D/2: 0.17–0.39 for pure wings and 0.14–0.37 with core dilution across τ_e = 2–3, and steeper with a transmitted core.

---

## 7. 图

| 图 | 意见 | 改前 | 改后 |
|---|---|---|---|
| 图 1a | 第一轮图表意见；N2 | QSO1 只画窄线自由的值（2107 km/s，空心菱形）；没有标出决定斜率的亮源 | QSO1 同时画两种分解：空心为窄线自由，实心为窄线固定 150 km/s，两者用虚线相连。3 个最亮 LRD 标出名字（EGS-42046、EGS-55604、EGS-49140）。图注相应修改 |
| 图 1b | 第一轮图表意见 | y 轴 540–7600 km/s，主体被压扁；未说明为何没有误差棒 | y 轴改为 640–5200 km/s（数据范围 1020–4090，下方留出图例位置）；图注注明 “error bars are omitted for clarity” |
| 图 2 | 第一轮图表意见；N1 | 单栏图，只有 6 种宽度定义加 LSF×0.5；红色带 0.05–0.08 容易被误读为置信区间 | 改为双栏两联图。(a) 6 种宽度定义加 BX26。(b) LRD 的 19 个分析变体，分 4 组（拟合、窄线、样本、重复光谱），空心三角标出单侧 95% 上限。删去红色带，改为保守上限 0.10 的点线。数据来自新文件 `data/results/review_checks.json` |
| 图 3 | 第一轮图表意见；N1 | 缺 45°<i<90° 的线（0.042），只有中位数和 95% 上限两条曲线 | 加上 45°<i<90° 的线和保守上限 0.10 的点线；图注同步 |
| 图 4c | M3；第一轮图表意见 | 参考线只画理论斜率 0.25，孤立点未点名 | 加两条点划线，分别是样本 A、B 在实测误差下光球模拟的中位偏斜率 0.14 和 0.21；标出 RUBIES-UDS-29813。图注同步 |

---

## 8. 表

- **表 2**（`scripts/make_tables.py`）：
  - 新增三行：“[O III]-anchored, absorption rule re-applied”（N = 27）、“Fit window ±8000 km s⁻¹, [S II] masked”、“Lowest-S/N passing spectrum”。
  - 注释 a 补充了 0.102（[O III] 重新筛选）、0.083（最低 S/N）以及随机效应 0.049 和 0.062 四个上限，并把 “[O III] hits the lower bound in 21 of the 29 LRDs” 改正为 “21 of the 23 LRDs where it is measured”。
  - 新增注释 e，说明近似方法。
  - SDSS 行原来手写 `$0.123$ ($0.112$--$0.136$)`，仍在表中；正文中的 `\Ssg` 宏改为从 `data/sdss/regression_results.json` 读取。
- **表 3**：
  - 第 2 行：
    - 改前：σ_int < 0.075 dex (< 0.06 dex with repeat floor)
    - 改后：σ_int < 0.075 dex (baseline; < 0.10 dex across analysis variants, < 0.049–0.062 dex with repeat-spectrum variance)
  - 第 3 行散射壳一栏（N3）：
    - 改前：predicts FWHM ∝ L^{0.09–0.39} across τ_e = 0.5–3 …, overlapping b
    - 改后：predicts FWHM ∝ L^{0.14–0.39} at τ_e ≃ 2–3 depending on R(L_Hα): matches b = 0.17, but the narrow-line-corrected b lies at or below its lower edge
  - 第 5 行：加入 “< 0.21 dex at the conservative ceiling”。
  - 表注：
    - 改前：“…the [O III]-anchored ceiling (0.097 dex) loosens them by 30%”
    - 改后：“…the conservative ceiling (0.10 dex) loosens them by 35%”

---

## 9. 统计与代码（我负责的部分）

### 新脚本 `analysis/14_review_checks.py`
从存储的逐源测量值计算审稿意见要求的全部新数值，写入 `data/results/review_checks.json`，并已加入 `Makefile` 的 `analysis` 目标。内容如下：

1. 19 个分析变体的中位数、95% 区间，以及在有逐源输入时的单侧 95% 上限。供图 2(b) 和第 4 节使用。
2. [O III] 绑定并重新应用吸收边界规则（近似方法见第 3 节）。
3. 常数线宽（b = 0）模型。
4. 去掉 RUBIES-UDS-29813 后的光球偏斜率，以及它在固定 L_Hα 下的 L_ph 偏离：样本 A 为 −0.48，样本 B 为 −0.50，在两个样本中都是最负的点。
5. 误差放大 1.5、2.0、2.8 倍，以及加性地板 0.044、0.06、0.08 dex，对本文样本和 BX26 的 Δσ_int 的影响。
6. 重复光谱诊断：G395H − G395M 宽度差、主光谱选择偏差、只用 G395M 主光谱的子样本、42232 换成 G395M 光谱。
7. SDSS 去掉 χ² 规则的结果，以及因 S/N 落选源的 FWHM 中位数。
8. A2744-QSO1 相对式 (3) 的残差复核，见第 11 节。
9. τ_e = 2–3 下散射壳的斜率预测。
10. CEERS-746 与 Skyfire-299 的角距离：35.2″。

### `analysis/13_parent_and_repeat_checks.py`
- **Σ₀₀ 的 bug（N7）**：
  - 改前：多光谱对象的协方差矩阵在两次减去 b·C_xy 之后，又执行了 `Sigma[..., 0, 0] += (B_sub * cxy0)`，使主光谱对角元只剩 −b·C_xy。
  - 改后：删除这一行，对角元为正确的 −2b·C_xy，与单光谱分支一致，并加了注释。
  - 影响：加性模型的上限从 0.04915 变为 0.04923，正文显示的 0.049 不变。
- **新增乘性误差模型**：`error_model="multiplicative"`，未建模方差写成 f²·C_yy，f 在 1–6 上取均匀先验。结果为 f 中位数 2.04，σ_int 中位数 0.041，95% 上限 0.062；42232 只保留主光谱时为 0.067。
- **新增自由光栅偏移模型**：函数 `ranef_free_offset`。G395H 偏移为 0.054 ± 0.022 dex，σ_int 上限 0.049，不变。
- **`exterior_potential_check.json` 的手写文献注释**：移进脚本，作为常量 `LITERATURE_OPACITY_CITATIONS`，并注明是手工录入。现在脚本输出与仓库文件一致。
- **`repeat_floor_slicings.tsv`**：用当前脚本重新生成。数值不变，行标签改为与 JSON 一致。

### `scripts/make_numbers.py`
- **删除手写的 `FsigCeilRep = '0.06'`。**
- **新增约 60 个宏**，全部从 `review_checks.json` 或 `repeat_random_effects.json` 读取。
- **`Ssg`、`Ssglo`、`Ssghi`、`NSsg`**：原来是手写的 `'0.123'` 等，改为从 `data/sdss/regression_results.json` 读取，数值相同。
- **`nBootNp`**：80 改为 200。
- **δ 与 λ 改用精确式**（N6）。
- **`narCeilLoosenPct`**：改为相对保守上限计算，30 变为 35，断言同步更新。
- **恢复每个宏的出处注释**：原脚本在写出前把注释清空了（`M = [(…, '') for …]`）。现在每行形如 `\newcommand{\FsigCons}{0.10}  % [C] conservative ceiling: …`，文件头说明了 R/C/A/D/L/V/X 各标记的含义。投稿用的去注释版本由 `make_submission_source.py` 自动去掉这些注释。

### 其他
- `scripts/make_submission_source.py`（N4）：双盲版本中，用正则把 `archived at\dataset[doi:…]{…}` 替换为 “archived in a public repository (DOI withheld for review)”，并断言输出中不含 “zenodo”。
- bib 键名：`DEugenio2025`、`Juodzbalis2025`、`Inayoshi2025` 改为 `…2026`，与 year 字段一致。修改范围是 `refs.bib`、`refs_ads_20261113.bib` 和 `ms.tex`。
- `scripts/check_reproduction.py`：通过时打印的 “767 macro values and 497 definitions” 是写死的计数，改为不含具体数字的提示。`data/measurements/paper_reference_hashes.json` 已更新。

---

## 10. 代码与可复现性（第二轮第 5 节的问题与上一轮的代码问题 (a)–(d)）

### 三个无法运行的新脚本（第二轮代码审计第 1 条）
- **改前**：三个脚本都引用了未定义的变量 `T2_DIR`、`DATA_DIR`、`CODE_DIR`、`HERE`、`M_DIR`（`parent49_recovery.py` 在导入时就失败），并依赖仓库中不存在的 `t1_output.json`、`S4_QSO1_fit.json`、`working/`。
- **改后**：所有路径改为相对仓库，每个缺失的输入都对应到仓库中已有的文件：

| 原来的输入 | 改为读取 | 一致性核对 |
|---|---|---|
| `t1_output.json["failed_21_audit"]` | `data/widths/parent_and_repeat_covariances.json["failed_clean_lrd_21"]` | 与 `parent_covariance_audit.json` 完全相同 |
| `S4_QSO1_fit.json` | `data/qso1/qso1_fit.json` | QSO1 的 `orig_*` 字段一致 |
| `working/<fid>.json` | `data/linefits/uniform_linefits.tsv` | null 模型的 BIC 取 `main_BIC + main_delta_BIC_null` |
| `M_DIR/M1_frozen_sources.json` | `analysis/injection/templates.py` 与 `frozen_sources.tsv` | 带光谱和误差数组的哈希校验 |
| `targets[fid]["path"]`（原为 KeyError） | `["input_path"]` | — |

  每个脚本在开始拟合前都会检查所需光谱是否存在、SHA-256 是否一致，缺失时提示运行 `fetch_jwst.py --download --qso1`。

- **核实**：
  - pyflakes 不再报告任何未定义名。三个脚本都能导入并读入全部仓库输入，然后在缺少光谱时正常停止。
  - 用 `qso1_fit.json` 中存储的 QSO1 像素做端到端重拟合：FWHM 相差在 4×10⁻⁵ km/s 以内；W150 与 W30 的重拟合和 `narrow_degeneracy_refit.json` 的差别在 10⁻⁵ 以内。
- **注释与代码不一致**（第二轮审计第 7 条），`narrow_core_refit.py:509–515`：注释已改为描述实际代码，即窄线固定的无宽线拟合中 `F_narrow` 作为通量上限，下限为 L = 10⁴⁰ erg/s。
  - 这一规则与 `data/jwst/parent49_criteria.json` 中的描述（“…ceiling / empirical LRD draw”）**不一致**：代码中并没有任何经验抽样。预注册文件带时间戳，我们没有改动它，请作者决定如何说明。
- `parent49_recovery.py` 并没有真正导入 `stage2_sloped_injection.py`，只在注释里提到过它。`sloped_injection_recovery.json` 在仓库中仍然没有生成脚本。

### (a) 质量判据与样本复现（第一轮代码问题 1）
- **改前**（`lrdwidth/fit_batch.py:46`）：`all_parameters_away_from_bounds=not f['bound_parameters']`，即任何参数撞边界都判为不合格。按这一判据，26/29 个 LRD 不合格。
- **改后**（`fit_batch.py:25–26、51–59`）：新增 `quality_tests()`，只在宽线宽度撞边界，或者采用了吸收屏且 `tau_abs`、`W_abs` 撞边界时才剔除，与 §2.3、SDSS 重拟合和注入代码一致。原来的标记保留为诊断字段，不参与筛选。
- **新增 `analysis/sample/make_selection_audit.py`**：从 `uniform_linefits.tsv` 重建 `selection_audit.json`。
  - **我已核实与仓库文件逐字节一致**，选出的光谱与 `frozen_sources.tsv` 相同，顺序也相同：70 → 50 → 44，即 29 个 LRD 加 15 个其他源。
  - 加 `--reduced-chi2-max 3` 时得到 79/54/46，与 `rule_variants.json` 一致。
  - 唯一的例外：文件中两个 SHA-256 字段指向开发环境中的文件，只能照录，无法重算。

### (b) 自助抽样次数（第一轮代码问题 2）
- `fit_batch.py`：80 改为 `NBOOT_NP=200`。正文的 `\nBootNp` 也改为 200（第 9 节）。
- 核实：`source_fits` 中 44 个源 × 2 种非参数宽度的 88 条记录，`N_requested` 都是 200，可用次数 179–200。80 次只对旧规则下通过的 9 个源跑过。
- 仍需说明：仓库中没有对观测光谱运行 `measure.np_bootstrap` 的那次运行记录和随机种子。

### (c) 变体缺少生成脚本（第一轮代码问题 3；第二轮审计第 3 条）
- `lrdwidth/linefit.load_line` 新增 `window`、`mask_rest_um`、`mask_halfwidth_kms` 三个参数，默认值下行为不变。
- `analysis/jwst/refit_spectra.py` 新增 `--lsf-scales`、`--window-kms`、`--mask-sii`、`--mask-halfwidth-kms`、`--nboot-np`、`--frozen-only` 等选项，每个输出都会记录所用设置。
- 新增 `analysis/jwst/make_lsf_half_inputs.py`，从 `source_fits_lsf_half/` 逐字节重建 `regression_inputs_lsf_half.json`。
- **仍未解决**：
  - ±8000 km/s 运行时 [S II] 掩膜的宽度没有记录；把窗口运行结果汇总成 `window_and_wing_checks.json` 的脚本也不在仓库中。
  - `parent_and_repeat_covariances.json` **无法从存储数据重建**。缺少的是 3/43 条 QC 光谱和 11/21 条失败光谱的通量–宽度协方差、`W_narrow`，以及需要光谱才能算的 `err_Ha_over_inband_median`。其余标量字段都与存储数据完全一致。

### (d) 定标表版本（第一轮代码问题 5）
- **改前**：`scripts/fetch_calibration.sh` 从 msaexp 的 `main` 分支下载，没有校验。
- **改后**：
  - 固定到 commit `9d342cd85f2dae2435f14c5201733a5e98d8d9a5`，它是最后一次修改这三张表的提交。
  - 每个文件下载后都对照 `third_party/external_data.tsv` 中记录的 SHA-256 校验，不一致就停止。这个 commit 下的三个文件与记录的哈希完全一致，因此版本是确定的，不是猜测。
  - `third_party/external_data.tsv` 中的 URL 也改为指向这个 commit。

### 字段名（第一轮次要 13）
- `B15_nonparam_fesc_pct_median` 改名为 `B15_nonparam_flux_below_vesc_pct_median`，因为 74% 是束缚部分的比例，100 − 74 = 26，落在正文的 23–28% 之内。
- 改动范围是 `data/widths/window_and_wing_checks.json`。`ancillary_checks.json` 由 `analysis/12` 原样复制这一块，重跑后同步更新。

### 文档
- 更新了 `README.md`、`analysis/REPRODUCING_MEASUREMENTS.md`、`data/README.md` 和 `analysis/reproduction_provenance.json`。
  - 原来写的是“所有数字自动生成”，§5 声称可以从外部 DJA 光谱复现，这两处都已改写。
  - 现在的表述是：所有统计分析都能从存储的逐源测量值复现；谱线重拟合和注入只需要 DJA 光谱和定标表。仍缺的生成脚本，如 `sloped_injection_recovery.json`、`membership_full_covariance.tsv` 等，单独列在新增的 §6 “Not scripted” 中。
- `make_numbers.py` 中仍有少数手写常数（如 `fKepNp` 38–44%、`fEscNp` 23–28%），README 已如实说明。

---

## 11. 对第二轮意见的更正

- **N2 最后一条有误。** 第二轮说正文 “−0.08 ± 0.03 dex from Equation (3)” 的参照关系写错了，实际是相对窄线固定的关系（a = 3.280，b = 0.109）。这次复核发现，宏 `\QSOresFix` 取的是 `residual_Eq3_delensed_dex = −0.0830`，它确实是相对式 (3)（a = 3.346，b = 0.170）计算的，我独立重算得到 −0.0830。−0.0796 是仓库中另一个字段，正文并没有用它。所以**这一条不需要修改**，正文保持原样。
- **M2 中的 “6/7”** 应为 7/7，见第 3 节。

---

## 12. 未处理的事项

以下几项需要原始光谱或开发环境中的中间文件，本次无法完成：
- **“LSF×0.5 + [O III] 绑定” 组合变体**（N1）。
- **窄线变体的精确质量筛选。** 目前用基线协方差近似，而且只重新应用了吸收边界规则；ΔBIC_null 的参照模型问题只在正文中加了一句说明。
- **在注入中加入吸收和超出形式误差的噪声。** 现在只在附录 B 中写明这两点局限。
- **49 模板注入通过率比值的不确定度**（第二轮代码审计第 6 条）。
- **第一轮第 6 条（BX26 维里质量）。** 用 GH05 链复算 “仅由宽度预期的弥散”（0.15 对 0.31）需要额外的推导，本次只在措辞上注明对照组含有异质的质量估计器。
- **期刊和篇幅。** 见第 0 节，由作者决定。
