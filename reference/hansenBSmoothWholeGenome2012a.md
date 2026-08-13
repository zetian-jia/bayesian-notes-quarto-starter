# BSmooth × Q1–Q5

| | BSmooth 内容 | 说明 |
|---|---|---|
| **Q1 潜变量** | $\pi_j$（单点均值）；DMR阶段 $\pi_{i,j}=f_i(l_j)=\alpha(l_j)+\beta(l_j)X_i+\varepsilon_{i,j}$ | 潜变量做了**函数型分解**：基线 $\alpha$ + 组间差异 $\beta$ + 生物学噪声 $\varepsilon$，$\beta(l_j)$ 才是检验目标 |
| **Q2 观测模型** | $M_j\sim\text{Binomial}(N_j,\pi_j)$ | 标准二项抽样，与 DSS 系列一致 |
| **Q3 先验结构** | (a) 均值光滑：局部似然（二次多项式+精度加权+tricube核）；(b) 方差光滑：经验SD→75分位截断→滑动平均 | **双重光滑**（均值+方差都光滑），比 DSS-single 更进一步；实现方式是局部加权GLM，非简单移动平均 |
| **Q4 后验计算** | 全程频率派：加权GLM → 均值差 $\hat\beta$ → 截断+滑动平均得 $\hat\sigma$ → $t=\hat\beta/\hat\sigma$ | **无贝叶斯后验**，比 DSS-single（还有EB收缩）更彻底地走近似/矩估计路线 |
| **Q5 推断验证** | 阈值 $c$ 按 $t$ 的经验分布选取；DMR判定叠加启发式规则（连续性、间距、最小长度、最小效应量） | **无正式 $p$ 值**，用"弱推断+强规则"代替，回避了双重使用数据的问题，但牺牲统计严谨性 |

## 关键创新点

1. **函数型分解**：$\beta(l_j)$ 独立建模为"组间差异曲线"，DMR = $\{\beta(l_j)\neq0\}$，比 DSS-single 的组均值曲线更结构化
2. **双重光滑**：均值和方差都做空间平滑，不止 DSS-single 那样只平滑均值
3. **精度+距离双重加权**的局部似然平滑，比简单滑动窗口平均统计效率更高
4. **核心论断**：覆盖度降低技术噪声，重复数降低生物学噪声 —— 二者不可互相替代
5. **启发式阈值+区域规则**代替严格假设检验，工程实用主义路线


## Materials and methods

### Datasets

The Lister data are from a WGBS experiment on the IMR90 fibroblast cell line. Six different library preparations were sequenced individually on an Illumina sequencer using up to 87 bp single-end reads and subsequently pooled to yield 25× coverage of CpGs.

The Hansen data are from a WGBS experiment on three paired tumor-normal colon samples, sequenced on ABI SOLiD using 50 bp single-end reads with a CpG coverage of 4×. These data were prepared and sequenced in the laboratory of AP Feinberg.

The Hansen-capture data comprise the same six samples as the Hansen data sequenced on an Illumina sequencer with up to 80 bp single reads, using a bisulfite padlock probe (BSPP) capture protocol, yielding a CpG coverage of 11× to 57× of 40,000 capture regions (one sample had substantially lower coverage than the rest, and the capture regions varied in efficiency). These data were prepared and sequenced in the laboratory of K Zhang.

The Tung data are from a WGBS experiment on peripheral blood mononuclear cells from six rhesus macaque individuals, three of high social rank and three of low social rank. The data were sequenced using an Illumina sequencer with 75 bp single end reads, yielding a CpG coverage of 11× to 14×.

The Lister data were created in the following way: we obtained the raw reads from the IMR90 cell line and aligned against the hg19 genome using Merman with iterative trimming. Prior to alignment, two bases were trimmed from the start of the read and one base from the end of the read. Based on our M-bias plots, we furthermore filtered the last ten bases of every read (based on its trimmed length), when we summarized the methylation evidence. Based on the quality control plots, the flowcells marked ECKER_1062 were discarded. These data form the basis for all analysis of the Lister data in the manuscript as well as Figures S1 to S4 in Additional file [1](https://link.springer.com/article/10.1186/gb-2012-13-10-r83#MOESM1).

## Smoothing

We denote the number of reads associated with the $j$ th CpG being methylated and unmethylated with $M_j$ and $U_j$, respectively. The CpG-level summary is simply the proportion $M_j/N_j$, with $N_j = M_j + U_j$ the coverage for the $j$ th CpG. We assume each $M_j$ follows a binomial distribution with success probability $\pi_j$. The success probability represents the true proportion of cells for which the $j$ th CpG is methylated in the sample being assayed. The proportion $M_j/N_j$ is an unbiased estimate of $\pi_j$ with standard error

$$
\text{SE}(\hat{\pi}_j) = \sqrt{\frac{\pi_j(1-\pi_j)}{N_j}}
$$

and we denote $\hat{\pi}_j$ the single-CpG methylation estimate of $\pi_j$.

- **$M_j$**：第 $j$ 个 CpG 位点上，支持"甲基化"状态的读段数（methylated reads）

- **$U_j$**：第 $j$ 个 CpG 位点上，支持"未甲基化"状态的读段数（unmethylated reads）

- **$N_j = M_j + U_j$**：第 $j$ 个 CpG 位点的总覆盖度（coverage），即甲基化读段与未甲基化读段之和

- **$M_j/N_j$**：CpG 位点层面的甲基化比例，是最直接、最粗糙的观测统计量

- **$\pi_j$**：第 $j$ 个 CpG 位点的**真实**甲基化比例——即样本中该位点实际处于甲基化状态的细胞比例，是无法直接观测的**潜变量**

- **$M_j \sim \text{Binomial}(N_j, \pi_j)$**（隐含的观测模型）
  假设 $M_j$ 服从以 $N_j$ 为试验次数、$\pi_j$ 为成功概率的二项分布——这是测序读段计数的抽样过程模型

- **$M_j/N_j$ 是 $\pi_j$ 的无偏估计**：即 $E[M_j/N_j] = \pi_j$，符合二项比例估计的基本性质

- **标准误 $\text{SE}(\hat\pi_j) = \sqrt{\pi_j(1-\pi_j)/N_j}$**
  单个位点估计值的不确定性，随覆盖度 $N_j$ 增大而减小——这是后续加权平滑时"权重"的来源：覆盖度越高（标准误越小），该位点在平滑中占的权重越大

- **$\hat\pi_j$**：单位点甲基化估计值（single-CpG methylation estimate），即基于 $M_j/N_j$ 得到的 $\pi_j$ 点估计

We furthermore assume that $\pi_j$ is defined by a smoothly varying function $f$ of the genomic location, that is, for location $l_j$,

$$
\pi_j = f(l_j)
$$

- **$\pi_j = f(l_j)$**
  假设真实甲基化比例 $\pi_j$ 是基因组坐标 $l_j$ 的**光滑函数** $f$ 的取值——与 DSS/DSS-single 中 $\mu_{ij}=f_j(l_i)$ 的空间平滑先验思想完全一致，都是用"空间光滑性"来正则化/借用邻近位点信息

We estimate $f$ with a local-likelihood smoother [28]. We start by choosing a genomic window size $h(l_j)$ for each $l_j$. The window is made large enough so that 70 CpGs are included but at least 2 kb wide. Within each genomic window we assume

$$
\log\left[\frac{f(l_j)}{1-f(l_j)}\right]
$$

- **局部似然平滑（local-likelihood smoother）**：估计 $f$ 的具体方法，属于非参数平滑技术的一种（区别于 DSS-single 用的简单滑动窗口平均）

- **窗口大小 $h(l_j)$**：以位点 $l_j$ 为中心的基因组窗口宽度，动态选取，要求窗口内至少包含 70 个 CpG 位点、且窗口宽度至少 2 kb——这是一种"自适应带宽"（adaptive bandwidth）策略，在 CpG 密度低的区域自动扩大窗口以保证信息量

- **$\log\left[\dfrac{f(l_j)}{1-f(l_j)}\right]$**
  对 $f(l_j)$ 做 **logit 变换**（对数几率），将取值范围为 $(0,1)$ 的比例映射到整个实数轴，便于用多项式回归建模——这是二项 GLM（logistic 回归）的标准处理方式

is approximated by a second degree polynomial. We assume that data follow a binomial distribution and the parameters defining the polynomial are estimated by fitting a weighted generalized linear model to the data inside the genomic window. For data points inside this window, indexed by $l_k$, weights are inversely proportional to the standard errors of the CpG-level measurements, $\text{SE}(\hat{\pi}_k)$, and decrease with the distance between the loci $|l_k - l_j|$ according to a tricube kernel (Figure 3a,b). Note that the smoothness of our estimated profile $\hat{f}$ depends on genomic CpG density. We recommend users adapt the algorithm's parameters when applying it to organisms other than human.

- **二次多项式近似（second degree polynomial）**：在每个窗口内，假设 logit 变换后的 $f$ 可以用关于位置的二次多项式局部逼近——这是局部多项式回归（local polynomial regression）的核心假设，比简单的"窗口内取平均"更灵活，能捕捉局部趋势（斜率、曲率）

- **加权广义线性模型（weighted GLM）**：在窗口内，用二项似然（binomial likelihood）拟合上述二次多项式的系数，属于 GLM 框架下的最大似然估计（而非纯粹频率派的滑动平均）

- **权重的两个来源**：
  1. 与标准误 $\text{SE}(\hat\pi_k)$ **成反比**——覆盖度越高、估计越可靠的位点，权重越大（精度加权，precision weighting）
  2. 与位点间距离 $|l_k - l_j|$ 按**三次核函数（tricube kernel**递减——距离中心位点越远，权重越小，这是局部回归中常见的核平滑权重函数

- **$\hat f$ 的平滑程度依赖于 CpG 密度**：由于窗口大小是按"CpG 个数"（70 个）而非固定物理距离定义的，在 CpG 稀疏区域窗口会自动变宽，因此平滑程度并非均匀，而是随基因组局部 CpG 密度自适应变化——这也是作者提示"应用到人类以外物种时需调整参数"的原因（不同物种 CpG 密度差异很大）

## Identification of differentially methylated regions

To find regions exhibiting consistent differences between groups of samples, taking biological variation into account, we compute a signal-to-noise statistic similar to the *t*-test. Specifically, we denote individuals with $i$ and use $X_i$ to denote group; for example, $X_i = 0$ if the $i$ th sample is a control and $X_i = 1$ if a case. The number of controls is denoted $n_1$ and the number of cases $n_2$. We assume that the samples are biological replicates within a group.

Similar to the previous section, we denote the number of reads for the $i$ th sample associated with the $j$ th CpG being methylated and unmethylated with $M_{i,j}$ and $U_{i,j}$, respectively. We assume that $Y_{i,j}$ follows a binomial distribution with $M_{i,j} + U_{i,j}$ trials and success probability $\pi_{i,j}$, which we assume is a sample-specific smooth function of genomic location $l_j$:

$$
\pi_{i,j} = f_i(l_j)
$$

- **$i$**：样本（个体）索引

- **$X_i$**：分组指示变量，$X_i=0$ 表示第 $i$ 个样本为对照组（control），$X_i=1$ 表示为病例组（case）

- **$n_1, n_2$**：对照组、病例组的样本数

- **"样本是组内生物学重复"**：这是本节方法与前述 DSS-single（无重复）方法的关键区别——这里假设每组内部有多个独立的生物学重复，因此可以直接估计组间/组内的生物学变异

- **$M_{i,j}, U_{i,j}$**：第 $i$ 个样本、第 $j$ 个 CpG 位点上，甲基化读段数与未甲基化读段数

- **$\pi_{i,j}$**：第 $i$ 个样本在第 $j$ 个 CpG 位点的真实甲基化比例（潜变量）

- **$\pi_{i,j} = f_i(l_j)$**
  假设每个**样本**都有自己独立的、随基因组坐标平滑变化的甲基化曲线 $f_i$——即平滑函数是**样本特异的**（sample-specific），而不是像前面章节那样只对组求一个共同的平滑函数


Furthermore, we assume that $f_i$ has the form

$$
f_i(l_j) = \alpha(l_j) + \beta(l_j) X_i + \varepsilon_{i,j}
$$

- **$f_i(l_j) = \alpha(l_j) + \beta(l_j)X_i + \varepsilon_{i,j}$**
  这是核心的**函数型线性模型**（functional linear model），把每个样本的甲基化曲线分解为三部分：
  - **$\alpha(l_j)$**：基线甲基化谱（baseline methylation profile），代表对照组（$X_i=0$ 时）的平均甲基化水平曲线
  - **$\beta(l_j)$**：两组间的真实差异函数（true group difference），这是本方法真正关心的目标量——**$\beta(l_j)\neq 0$ 的区域就是 DMR**
  - **$\varepsilon_{i,j}$**：残差项，代表样本间的生物学变异（biological variability），不是测序噪声

Here $\alpha(l_j)$ represents the baseline methylation profile and $\beta(l_j)$ the true difference between the two groups. The latter is the function of interest, with non-zero values associated with DMRs. The $\varepsilon_{i,j}$s represent biological variability with the location-dependent variance

$$
\text{var}(\varepsilon_{i,j}) \equiv \sigma^2(j)
$$

- **$\text{var}(\varepsilon_{i,j}) \equiv \sigma^2(j)$**
  假设生物学变异的方差本身也是基因组位置的**光滑函数**——即不同区域的生物学变异程度可以不同，但变化是连续、平滑的

- **关键结论："增加测序深度不能降低 $\varepsilon$ 带来的变异，只有增加生物学重复数才能"**
  这是全文最重要的统计学直觉之一：覆盖度（coverage）只影响单样本内部对 $\pi_{i,j}$ 估计的**技术噪声**（二项抽样误差），而 $\varepsilon_{i,j}$ 代表的是**样本间真实存在的生物学差异**，二者是两个完全不同来源的不确定性，只能靠"更多独立样本"来降低后者

assumed to be a smooth function. Note that increasing coverage does not reduce the variability introduced by $\varepsilon$; for this we need to increase the number of biological replicates.

We use the smoothed methylation profiles described in the previous section as estimates for the $f_i$, denoted $\hat{f}_i$. We estimate $\alpha$ and $\beta$ as empirical averages and difference of averages:

$$
\hat{\alpha}(l_j) = \frac{1}{n_1}\sum_{i:X_i=0}\hat{f}_i(l_j)
$$

$$
\hat{\beta}(l_j) = \frac{1}{n_2}\sum_{i:X_i=1}\hat{f}_i(l_j) - \frac{1}{n_1}\sum_{i:X_i=0}\hat{f}_i(l_j)
$$

- **$\hat f_i$**：用前一节局部似然平滑方法得到的、每个样本各自的平滑甲基化曲线估计值，作为 $f_i$ 的估计

- **$\hat\alpha(l_j)$**：对照组内所有样本平滑曲线 $\hat f_i(l_j)$ 的简单算术平均——基线水平的估计

- **$\hat\beta(l_j)$**：病例组平均曲线减去对照组平均曲线——组间差异的估计，是**均值差**（difference of means），而非某种复杂的模型系数

- **标准差估计的三步流程（提高稳健性的处理）**：
  1. 先计算两组内的**经验标准差**（empirical SD）作为 $\sigma(l_j)$ 的初始估计
  2. **下限截断（flooring）**：把标准差在其第 75 百分位数处"封底"（即小于该分位数的值都提升到该分位数），这是借鉴 SAM（significance analysis of microarrays）类方法的经验做法，用于防止某些位点因偶然估计出极小的方差、从而在后续统计量中被人为放大

  3. **滑动平均平滑（running mean，窗口宽度 101）**：对截断后的标准差序列再做一次窗口平滑，得到最终的 $\hat\sigma(l_j)$，进一步降低估计噪声、保证空间连续性


To estimate the smooth location-dependent standard deviation, we first compute the empirical standard deviation across the two groups. To improve precision, we used an approach similar to [30]: we floored these standard deviations at their 75th percentile. To further improve precision, we smoothed the resulting floored values using a running mean with a window size of 101. We denote this final estimate of local variation with $\hat{\sigma}(l_j)$.

We then formed signal-to-noise statistics:

$$
t(l_j) = \frac{\hat{\beta}(l_j)}{\hat{\sigma}(l_j)}
$$

- **$t(l_j) = \hat\beta(l_j)/\hat\sigma(l_j)$**
  **信噪比统计量**（signal-to-noise statistic），形式上类似经典 *t*-检验统计量：分子是组间差异估计（信号），分母是局部标准差估计（噪声）。$|t(l_j)|$ 越大，说明该位点存在真实差异的证据越强

To find DMRs, that is, regions for which $\beta(l_j) \neq 0$, we defined groups of consecutive CpGs for which all $t(l_j) > c$ or $t(l_j) < -c$ with $c > 0$ a cutoff selected based on the marginal empirical distribution of $t$. We adapted our algorithm so that CpGs further than 300 bp apart were not permitted to be in the same DMR.

We recommend including in the procedure only CpGs that have some coverage in most or all samples. Furthermore, we recommend filtering the set of DMRs by requiring each DMR to contain at least three CpGs, have an average $\beta$ of 0.1 or greater, and have at least one CpG every 300 bp.

- **DMR 判定规则**：寻找连续 CpG 位点组成的区域，使得区域内所有位点都满足 $t(l_j) > c$ 或都满足 $t(l_j) < -c$，其中阈值 $c>0$ 根据 $t$ 统计量的边际经验分布来选取（一种数据驱动的阈值设定方式，而非固定的理论分位数）

- **300 bp 间距限制**：算法规定，若相邻两个满足条件的 CpG 位点之间物理距离超过 300 bp，则不允许划入同一个 DMR——防止把基因组上相隔很远、偶然都显著的位点错误地合并成一个区域

- **后处理过滤标准（推荐规则）**：
  - 只纳入在**大多数或全部样本中都有一定覆盖度**的 CpG 位点（保证数据质量）
  - 每个 DMR 至少包含 **3 个 CpG 位点**（避免单点噪声被误判为区域性差异）
  - 平均 $\hat\beta \geq 0.1$（要求差异幅度有实际生物学意义，而不仅是统计显著）
  - 区域内每 300 bp 至少有一个 CpG 位点（保证区域内部密度、避免"虚假拼接"的稀疏区域）