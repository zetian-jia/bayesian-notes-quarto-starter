# DSS-single 论文的 Q1–Q5 框架总结

| | 论文内容 | 框架对应 |
|---|---|---|
| **Q1 潜变量** | $p_{ij}$（真实甲基化比例）、$\mu_{ij}$（beta 均值，进一步等于 $f_j(l_i)$）、$\phi_{ij}$（beta 离散度） | 潜变量本身是**分层的**：$p_{ij}$ 是最底层的"个体"潜变量，而 $\mu_{ij}$、$\phi_{ij}$ 是控制 $p_{ij}$ 分布的"超参数"，其中 $\mu_{ij}$ 又被进一步结构化为空间坐标 $l_i$ 的函数 $f_j(l_i)$ |
| **Q2 观测模型** | $X_{ij}\mid N_{ij},p_{ij}\sim\text{Binomial}(N_{ij},p_{ij})$ | 标准二项抽样，覆盖度 $N_{ij}$ 从这里进入模型，是测序深度不均一的技术噪声来源 |
| **Q3 先验结构** | 三种先验同时使用：(a) $p_{ij}\sim\text{Beta}(\mu_{ij},\phi_{ij})$——**收缩到 $\mu_{ij}$**，捕捉生物学 overdispersion；(b) $\mu_{ij}=f_j(l_i)$——**光滑**（沿基因组坐标）；(c) $\phi_{ij}\sim\text{log-normal}(m_{j0},r_{j0}^2)$——**收缩到组内共同分布**（跨位点信息共享，本质是经验贝叶斯的先验） | 这篇论文的核心结构特点是：**光滑先验（Q3-光滑）被用来间接解决 Q1 层面本该由"重复"提供的信息缺口**——这是它区别于普通 DSS 的关键 |
| **Q4 后验计算** | $\hat\mu_{ij}$ 用移动平均（非贝叶斯闭式解，是一种确定性/频率派的平滑估计）；$\phi_{ij}$ 用经验贝叶斯（矩估计 + 跨位点信息借用，而非完整后验采样）；最终推断用 **Wald 检验**而非完整后验分布 | 计算策略明显是"轻量级"路线：不做 MCMC/VI，全部用近似矩估计和渐近正态（Wald）拼接而成，牺牲了严格贝叶斯的完整后验，换取全基因组规模下的可扩展性 |
| **Q5 推断验证** | 模拟研究（已知 ground truth 的 spiked-in DMR）+ 真实数据的独立外部基准（DHS、CGI、DEG、ChromHMM）交叉验证 | 用**外部生物学证据**代替传统的置换检验/数据划分，这是基因组学里常见但统计上"软"的验证方式——它验证的是"生物学合理性"而非严格的统计校准（calibration） |

---

## 创新点

1. **用空间平滑先验替代生物学重复**（Q3 光滑先验去补 Q1 的信息缺口）
   这是全文最核心的思想转折。经典 beta-binomial 模型（DSS 原版）依赖多个生物学重复来估计 $\phi_{ij}$；但 DSS-single 意识到，如果假设 $\mu_{ij}=f_j(l_i)$ 是空间光滑的，那么相邻 CpG 位点在均值已知的前提下就可以彼此充当"伪重复"（pseudo-replicates）。这本质上是把"生物学重复"（跨样本的信息）替换成了"空间重复"（跨位点的信息）——一种维度替换的技巧。

2. **经验贝叶斯离散度估计不要求重复**
   关键洞察在于：EB 方法需要的是 $\mu_{ij}$ 已知，而不是"多次独立观测同一个 $p_{ij}$"。只要 $\hat\mu_{ij}$ 能被平滑估计出来，单个观测值 $X_{ij}$ 也能提供关于方差的信息（因为 beta-binomial 的方差结构中 $\phi_{ij}$ 和 $\mu_{ij}$ 是可分离的）。

3. **工程上的效率优先**：用简单滑动窗口平均代替 BSmooth 的样条平滑，证明"简单方法几乎一样好，但快得多"——这是方法论上一个重要但容易被忽视的贡献：**不是每个环节都需要用最复杂的模型**，Q4 层面的近似选择要服务于全基因组规模的可扩展性。

## 难点/关键点

1. **信息借用的合法性依赖于"光滑"假设本身是否成立**
   整个方法的地基是"甲基化水平沿基因组连续变化"，如果某处真实存在陡峭的边界（比如一个很短、很极端的 DMR），平滑窗口会把边界"抹平"，导致 $\hat\mu_{ij}$ 估计有偏——这是 Q3 光滑先验固有的偏差-方差权衡（bias-variance tradeoff），论文中"最小区域长度"这类后处理规则某种程度上是在弥补这个问题，而不是从模型层面解决。

2. **Wald 检验的方差修正是全文技术难度最高的部分**
   文中提到"modify the variance calculation to account for smoothing effects"，但细节推到了补充材料。这提示了一个关键点：**一旦 Q4 阶段用了平滑（一种非独立的估计过程），Q5 阶段的检验统计量方差就不能再用"朴素"公式计算**，因为平滑引入了位点间的估计相关性，若忽略这一点会导致 $p$ 值虚假地过小（这正是框架里 Q5 独立轴反复强调的"用同一份数据先平滑（相当于'看过'数据）再检验"的经典陷阱——某种程度上这里存在轻微的"双重使用数据"的风险，只是论文声称已经做了修正）。

3. **无重复本质上是信息量的硬约束，光滑先验只能部分补偿**
   论文自己也承认"still preferable to have biological replicates"——光滑先验能提供的信息终究是有限的替代品，尤其在真实生物学变异（而非测序噪声）本身就很大的位点，用相邻位点做伪重复无法反映该位点独有的生物学随机性。

4. **Q5 的验证方式是间接的**
   模拟部分（已知 ground truth）能验证统计性质，但真实数据部分依赖 DHS/CGI/DEG 等外部标记的"重叠富集"作为代理验证——这只能说明检测结果"生物学上合理"，并不能严格验证 $p$ 值的统计校准是否正确（即 Q5 意义上"检验是否诚实"这件事，本质上仍是悬而未决的，只是被生物学证据侧面支持了）。



# Modeling whole genome BS-seq data

For WGBS data from two groups and one replicate in each group, we use the following notation: At the *i*th CpG site and *j*th treatment group (*j* = 1, 2), let $X_{ij}$ be the count of methylated reads, and let $N_{ij}$ be the total read count. Denote the underlying 'true' methylation proportion by $p_{ij}$. We have previously shown that it is reasonable to model $X_{ij}$ as a beta-binomial distribution, which captures both the biological and technical variation in the counts. The beta distribution is parameterized by mean ($\mu_{ij}$) and dispersion ($\phi_{ij}$), where $\phi_{ij}$ represents the biological variance among replicates in the same treatment group. Further, a log-normal prior is imposed on $\phi_{ij}$ in order to borrow information from all CpG sites in estimating the site-specific dispersions. To incorporate the spatial correlation in methylation levels, we assume that the underlying mean of the beta distribution, $\mu_{ij}$, varies smoothly across the genome. That is, we assume $\mu_{ij} = f_j(l_i)$, where $l_i$ denotes the genomic coordinate of the *i*th CpG site, and $f_j$ is a smooth function. Putting all of these pieces together, the data generated from WGBS experiments can be described by the following hierarchical model:

$$
X_{ij} \mid N_{ij}, p_{ij} \sim \text{Binomial}(N_{ij}, p_{ij})
$$

$$
p_{ij} \mid \mu_{ij}, \phi_{ij} \sim \text{Beta}(\mu_{ij}, \phi_{ij})
$$

$$
\phi_{ij} \sim \text{log-normal}(m_{j0}, r_{j0}^2)
$$

$$
\mu_{ij} = f_j(l_i)
$$


**公式与内容逐条解释**

- **$X_{ij}$**：第 $i$ 个 CpG 位点、第 $j$ 个处理组中，甲基化读段的计数（methylated read count）

- **$N_{ij}$**：第 $i$ 个 CpG 位点、第 $j$ 个处理组中，总读段计数（total read count）

- **$p_{ij}$**：潜在的"真实"甲基化比例（true methylation proportion），是无法直接观测的潜变量

- **$X_{ij} \mid N_{ij}, p_{ij} \sim \text{Binomial}(N_{ij}, p_{ij})$**
  给定总读段数 $N_{ij}$ 和真实甲基化比例 $p_{ij}$，观测到的甲基化读段数服从二项分布。这是测序过程本身的抽样模型（技术层面的随机性）

- **$\mu_{ij}$**：Beta 分布的均值参数，代表位点 $i$、组 $j$ 的期望甲基化水平

- **$\phi_{ij}$**：Beta 分布的离散度参数，代表同一处理组内、不同生物学重复之间甲基化水平的变异程度（即生物学变异）

- **$p_{ij} \mid \mu_{ij}, \phi_{ij} \sim \text{Beta}(\mu_{ij}, \phi_{ij})$**
  真实甲基化比例 $p_{ij}$ 本身也是随机的，服从以 $\mu_{ij}, \phi_{ij}$ 为参数的 Beta 分布。这一层把生物学重复间的变异（overdispersion）纳入模型，使得 $X_{ij}$ 边际上服从 beta-binomial 分布，而非简单二项分布

- **$m_{j0}, r_{j0}^2$**：log-normal 先验的位置参数与尺度参数（组特异，$j=1,2$ 各自不同）

- **$\phi_{ij} \sim \text{log-normal}(m_{j0}, r_{j0}^2)$**
  对离散度参数 $\phi_{ij}$ 施加 log-normal 先验，目的是"借用"全基因组所有 CpG 位点的信息（across all sites in the same group），使得单个位点的离散度估计更加稳健、避免因单点数据稀疏而估计不准

- **$l_i$**：第 $i$ 个 CpG 位点在基因组上的坐标（genomic coordinate）

- **$f_j$**：一个平滑函数（smooth function），依赖于组别 $j$

- **$\mu_{ij} = f_j(l_i)$**
  假设均值参数 $\mu_{ij}$ 是基因组坐标 $l_i$ 的平滑函数值，而不是每个位点独立估计。这一设定编码了甲基化水平沿基因组的**空间相关性**（相邻 CpG 位点甲基化水平相近）

This is a comprehensive model that captures all three of the important characteristics in the WGBS-seq data discussed above. The binomial distribution captures the random sampling process of the BS-seq experiment, the beta distribution models the biological variation among replicates, and the smooth function accounts for the spatial correlation among nearby CpG sites. The log-normal prior of the dispersion combines information from CpG sites genome-wide, which provides a basis for information sharing and helps the estimation of dispersion.

## Smoothing procedure

We use a simple moving average procedure of the collapsed counts to estimate $f_j$. Specifically, at the *i*th CpG site, we estimate the mean by $\hat{\mu}_{ij} = \sum_{l \in S_i} X_{lj} / \sum_{l \in S_i} N_{lj}$, where $S_i$ is a set of CpG sites within a user-defined window of size $w$ (500 base pairs by default), e.g. $S_i = \{m : |l_m - l_i| < w\}$. We will show below that the simple procedure performs almost as well as more complicated, spline-based smoothing from BSmooth, yet it is much more computationally efficient.

- **$\hat{\mu}_{ij}$**：第 $i$ 个 CpG 位点、第 $j$ 个组的均值估计值（估计得到的甲基化比例）

- **$S_i$**：以第 $i$ 个位点为中心、用户自定义窗口宽度 $w$（默认 500 bp）内的 CpG 位点集合，即 $S_i = \{m : |l_m - l_i| < w\}$

- **$\hat{\mu}_{ij} = \sum_{l \in S_i} X_{lj} / \sum_{l \in S_i} N_{lj}$**
  用窗口内所有位点的甲基化读段计数之和，除以总读段计数之和，得到该窗口的**移动平均**（moving average）估计，作为位点 $i$ 均值的估计

- **方法特点**：这是一种简单的滑动窗口平均法，用来估计平滑函数 $f_j$；相比 BSmooth 中更复杂的样条平滑（spline-based smoothing），效果几乎相当，但计算效率更高


## Dispersion estimation

With $\hat{\mu}_{ij}$, the dispersion parameters $\phi_{ij}$ are estimated through an empirical Bayes (EB) procedure as proposed in (11). The procedure borrows information from all CpG sites, and provides more accurate estimates of the dispersion. The EB procedure does not require replicated data as long as $\mu_{ij}$ is available. This makes sense because with the mean methylation known, even one observed data point can provide some information for the variance. Taking advantage of the spatial correlation in methylation levels, the means can be estimated through smoothing. So intuitively, when there is no replicate, data from nearby CpG sites can serve as 'pseudo-replicates'. We will show below that the procedure works well in both simulation and real data analysis. Although it is still preferable to have biological replicates, DSS-single works better than methods that completely ignore the variance, e.g. methods that use the differences in means or Fisher's exact test to call DMRs.

- **经验贝叶斯（EB）估计**：在已知 $\hat{\mu}_{ij}$ 的前提下，用经验贝叶斯方法估计离散度参数 $\phi_{ij}$，该方法借用所有 CpG 位点的信息，使估计更准确

- **无需重复的原因**：EB 方法只要求 $\mu_{ij}$ 已知，即使只有一个观测数据点，在均值已知的情况下也能提供关于方差的信息

- **"伪重复"（pseudo-replicates）思想**：利用甲基化水平的空间相关性，通过平滑得到均值估计后，相邻 CpG 位点的数据可以充当"伪重复"，弥补没有真实生物学重复的不足

- **优越性**：尽管有真实生物学重复仍然更好，但 DSS-single 比那些完全忽略方差的方法（如仅比较均值差异，或用 Fisher 精确检验判定 DMR）表现更优


## DMR calling algorithm

We use the following algorithm for DMR detection from WGBS data. The inputs for the algorithm are $X_{ij}$, $N_{ij}$ and $l_j$. We first perform local smoothing on estimated methylation proportions to obtain estimates for $\mu_{ij}$. Next, we estimate the dispersions through the EB procedure described above. We then identify DML (differentially methylated loci) by performing a hypothesis test: $H_0: \mu_{i1} = \mu_{i2}$ for the equality of mean methylation levels in two groups at each CpG site. To do this, we adopt the Wald test procedure proposed in (11), and modify the variance calculation to account for smoothing effects (details provided in Supplementary Materials). The Wald test is performed at each CpG site, and *P*-values are obtained from the test statistics. Finally, a user-defined *P*-value threshold and additional criteria such as minimum region length are applied to define DMRs.

- **输入数据**：算法输入为 $X_{ij}$（甲基化读段计数）、$N_{ij}$（总读段计数）和 $l_j$（基因组坐标）

- **步骤一：局部平滑**：对估计的甲基化比例进行局部平滑，得到 $\mu_{ij}$ 的估计值

- **步骤二：离散度估计**：通过前述经验贝叶斯（EB）方法估计离散度参数

- **步骤三：假设检验（识别 DML）**：在每个 CpG 位点检验原假设 $H_0: \mu_{i1} = \mu_{i2}$（即两组均值甲基化水平相等），采用 **Wald 检验**，并针对平滑带来的方差变化进行了修正（细节见补充材料）

- **步骤四：得到 P 值**：每个 CpG 位点进行 Wald 检验后，从检验统计量得到 P 值

- **步骤五：定义 DMR**：结合用户设定的 P 值阈值，以及最小区域长度等附加条件，最终划定差异甲基化区域（DMR）

## WGBS data sets

We analyze several WGBS data sets generated by Roadmap Epigenomics projects (22), including H1 (human embryonic stem cells) and IMR90 (human fibroblasts) cell lines, as well as human liver and hippocampus. The H1 and IMR90 data were obtained from the Gene Expression Omnibus (GEO) with accession number GSE16256, and the hippocampus and liver data are also from GEO with accession number GSE64577. There are two replicates available for each sample, but we limit our analyses to single-replicate comparisons. We perform analyses between H1 versus IMR90, and liver versus hippocampus to evaluate the DMR calling results. The H1 data are also used as template to simulate realistic WGBS data.

For H1-IMR90 comparison, the benchmarks are created as follows. After obtaining the DNase-seq data for H1 and IMR90 from ENCODE data (23), we apply MAnorm (24) to compare them and generate differential DNase I hypersensitive sites (DHSs). CpG island (CGI) data were downloaded from UCSC genome browser (25), with CGI shores defined as ± 1000 base pairs of each side of a CGI. Gene expression for H1 and IMR90 was obtained from (6), and the RPKM values are downloaded from the Human DNA Methylome website at the Salk Institute. We define differentially expressed genes (DEGs) as regions with absolute log2 fold changes of RPKMs greater than 1. The promoter regions for DEGs are defined as the regions ± 5000 base pairs of the transcriptional start sites. The genome segmentation by ChromHMM for H1 cells are obtained from ENCODE.

For liver-hippocampus comparisons, we also utilize available gene expression data and define DEGs using the same approach. However, since the DNase-seq data for these samples are unavailable, differential DHSs cannot be defined. Instead, we use the list of DNase I Hypersensitivity Clusters (DHCs) obtained from ENCODE, which is based on the union of DNase-seq peaks from 125 cell types. We use this list as a benchmark to assess the DMR calling accuracies under the assumption that this list contains active genomic regions for all biological conditions, and thus the DMRs are more likely to overlap with these regions.

- **Roadmap Epigenomics projects (22)**：数据来源项目，提供多种细胞/组织类型的表观基因组数据

- **H1 / IMR90**：两种细胞系，H1 为人胚胎干细胞（human embryonic stem cells），IMR90 为人成纤维细胞（human fibroblasts），二者用于 DMR 检测方法的评估比较

- **liver / hippocampus**：人肝脏与海马组织样本，作为另一组比较对象

- **GSE16256 / GSE64577**：GEO 数据库中对应的数据集编号（accession number），分别对应 H1/IMR90 数据和肝脏/海马数据

- **单重复比较（single-replicate comparisons）**：尽管每个样本实际有两个重复，但本文分析中刻意只使用单个重复，以模拟无重复场景（呼应 DSS-single 的设计目的）

- **DNase-seq / DHS（DNase I hypersensitive sites）**：DNase 高敏感位点数据，用于标记开放染色质区域，作为评估 DMR 检测结果生物学意义的独立基准（benchmark）

- **MAnorm (24)**：用于比较两组 DNase-seq 数据、识别差异 DHS 的统计方法/工具

- **CGI（CpG island）与 CGI shores**：CpG 岛及其"岸"区域（定义为 CGI 两侧各 ±1000 bp），是甲基化研究中常用的功能性基因组区域标记

- **DEG（differentially expressed genes，差异表达基因）**：定义为两组间 RPKM 值的 log2 倍数变化绝对值 > 1 的基因

- **RPKM**：基因表达定量指标（每千碱基每百万读段的读段数），来自 Salk 研究所人类甲基化组网站

- **promoter region（启动子区域）**：定义为转录起始位点（TSS）两侧各 ±5000 bp 的区域，用于关联 DEG 与甲基化变化

- **ChromHMM**：基因组分割（染色质状态注释）工具，为 H1 细胞提供功能区域标签，来自 ENCODE

- **DHC（DNase I Hypersensitivity Clusters）**：当某些样本（如肝脏、海马）缺乏自身 DNase-seq 数据时，改用 ENCODE 提供的、汇总 125 种细胞类型 DNase-seq 峰的"聚类"列表作为替代基准，假设该列表能代表所有生物学条件下的活跃基因组区域


## Simulation settings

All simulations are based on WGBS data from the H1 cell line. In comparisons of smoothing procedures and dispersion estimates, we select 20,000 contiguous CpG sites on chromosome 1, smooth the counts using BSmooth with different spans, and treat the smoothed values as the true $\mu$ in a hypothetical genome region. In each simulation, we generate counts based on the beta-binomial model described above, with $\varphi$ independently generated from the log-normal(−2.5, 1) distribution and an average read depth of 10x.

For comparisons of DML and DMR calling in simulations, data are generated for 100,000 CpG sites. We first obtain the true $\mu$ parameters for the first treatment group by smoothing the data from H1 ESC using BSmooth with smoothing span of 500 bps. We then generate the true $\mu$ parameters for the second treatment group with 100 DMRs 'spiked in' as follows. Using the original H1 data, we randomly generate 100 DMRs with lengths uniformly distributed between 5 and 50 CpGs. To generate these DMRs, we first randomly select 100 regions as 'target regions', and then select another set of 100 random regions as 'source regions'. We obtain the counts (X and N) from the source regions and replace values in the target regions with the counts from the source regions. We then smooth the data as above using BSmooth, and take the results as true $\mu$ parameters for the second group. We use this approach to guarantee that the true $\mu$ in the second group is smooth even after spiking in DMRs. Under this simulation setting, about 5% of the CpG sites lie in the DMRs. As above, we then generate counts from the true $\mu$'s in each group by independently simulating the dispersion parameters $\varphi$ from a log-normal(−2.5, 1) distribution, and then generating the counts based on the beta-binomial model.

- **20,000 个连续 CpG 位点**：用于比较不同平滑方法（smoothing）和离散度估计方法的模拟数据规模，取自 1 号染色体

- **BSmooth**：一种基于样条的平滑方法，用于生成模拟中的"真实" $\mu$（即把 BSmooth 平滑后的结果当作 ground truth）

- **span（平滑跨度/窗口宽度）**：BSmooth 平滑时使用的窗口大小参数，模拟中比较了不同 span 设置下的表现

- **$\varphi \sim \text{log-normal}(-2.5, 1)$**：模拟中离散度参数的生成分布，位置参数为 −2.5，尺度参数为 1（与前述 DSS 模型中 $\phi_{ij}\sim\text{log-normal}(m_{j0},r_{j0}^2)$ 的先验设定相呼应）

- **平均测序深度 10x**：模拟读段总数 $N_{ij}$ 的期望覆盖度设置

- **100,000 个 CpG 位点**：用于比较 DML/DMR 检测性能的更大规模模拟数据集

- **第一组真实 $\mu$**：直接对真实 H1 ESC 数据用 BSmooth（平滑跨度 500 bp）平滑得到

- **100 个"人为植入"的 DMR（spiked-in DMRs）**：模拟中人为制造的差异甲基化区域，长度在 5–50 个 CpG 位点之间均匀分布

- **target regions / source regions（目标区域/源区域）**：分别随机选取 100 个区域作为"目标区域"和另外 100 个区域作为"源区域"；用源区域的读段计数 $(X, N)$ 替换目标区域原有数值，从而人为制造差异

- **对第二组结果再次平滑**：将替换后的数据再次用 BSmooth 平滑，以保证"植入 DMR 后"第二组的真实 $\mu$ 仍然是空间平滑的，符合模型假设

- **约 5% 的位点落在 DMR 内**：该模拟设置下 DMR 占全部 CpG 位点的比例

- **计数生成流程**：由每组真实 $\mu$ 出发 → 独立生成离散度 $\varphi \sim \text{log-normal}(-2.5,1)$ → 依据 beta-binomial 模型生成最终的读段计数，与正文层级模型完全对应