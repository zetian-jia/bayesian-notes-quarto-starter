# MethSCAn × Q1–Q5

| | 内容 | 说明 |
|---|---|---|
| **Q1 潜变量** | 单细胞层面 $x_{ij}\in\{0,1,\text{NA}\}$ → 位点均值 $\bar x_i$/$\tilde x_i$ → 区间收缩残差 $r_{Ij}$ → PCA 隐空间 $\mathbf{x}_j^P$ | 潜变量是**多层聚合链**：细胞×位点(稀疏0/1/NA) → 位点(跨细胞) → 区间(跨位点，单细胞) → 低维隐空间(跨区间)，逐层降维压缩稀疏信号 |
| **Q2 观测模型** | $x_{ij}\in\{0,1,\text{NA}\}$，无显式概率分布，NA 直接编码缺失 | 与 DSS/BSmooth 不同，**不建二项似然**，直接对二值+缺失做确定性代数运算（矩估计式），跳过参数化观测模型 |
| **Q3 先验结构** | (a) 位点均值光滑：tricube 核加权（同 BSmooth）；(b) 区间残差收缩：分母 $n{+}1$ 伪计数收缩到 0（隐式贝叶斯先验） | **光滑 + 收缩到常数** 两种先验词汇同时使用；收缩机制专门应对单细胞覆盖度极不均匀问题，是 bulk 方法没有的设计 |
| **Q4 后验计算** | 全程无 MCMC/VI：核平滑、收缩均值、Welch t、PCA 均为**确定性闭式/迭代算法**；缺失值用**迭代 PCA**（EM 式自洽迭代）填补 | 计算策略是"轻量确定性算法+迭代自洽"，PCA 缺失值填补本质是一种 **EM 思想**（用低秩重构填补缺失，再重新分解） |
| **Q5 推断验证** | VMR/DMR 用滑窗+前2%阈值+区间融合；DMR 显著性用**置换检验**构造经验零分布，BH 法的无 p 值形式控制 FDR | **不依赖理论分布**（因融合步骤破坏独立性假设），改用置换重建整个"扫描+融合"流程的空分布——是本文 Q5 层面最严谨的部分；明确声明推断范围仅限"样本内"（cell 是统计单位，非 sample） |

---

## 关键创新点

1. **稀疏矩阵表示 + 双索引结构（$C_i$, $G_j$）**：专为单细胞数据设计的存储与索引范式，CSR 压缩支持高效区间查询
2. **两级聚合**：先跨细胞平滑得 $\tilde x_i$（位点趋势），再用残差 $r_{ij}=x_{ij}-\tilde x_i$ 衡量单细胞偏离——把"总体趋势"和"个体偏差"显式分离
3. **收缩均值（+1 伪计数）**：轻量级贝叶斯收缩，专门解决单细胞区间覆盖稀疏导致的估计不稳定
4. **VMR 概念**：无监督地找"细胞间变异大"的区域，独立于分组比较，是单细胞特有的探索性分析工具（bulk 方法没有对应概念）
5. **迭代 PCA 缺失值填补**：用 PCA 自身的低秩重构迭代填补缺失，避免标准 SVD 算法无法处理 NA 的限制
6. **置换检验+经验 FDR**：正面解决"扫描+区间融合"破坏传统统计假设的问题，是全文统计严谨性的核心
7. **邻居保真度评分（mean neighbor score）**：用 $k$ 近邻同类型比例量化降维/聚类效果，作为下游分析质量的量化基准

## Methods

### Raw data

Let us write $x_{ij}$ for the methylation status of CpG $i$ in cell $j$. The index $i$ runs over all CpG positions present in the genome, and the index $j$ over all cells in the assay. We write $x_{ij} = 0$ if position $i$ was found to be unmethylated in cell $j$ by bisulfite sequencing, $x_{ij} = 1$ if it was methylated and $x_{ij} = \text{NA}$ if position $i$ was not covered by reads from cell $j$ and the methylation status is therefore not available (NA).

---

- **$x_{ij}$**：第 $j$ 个细胞中、第 $i$ 个 CpG 位点的甲基化状态（这是单细胞层面的**原始观测**，不同于 bulk WGBS 中的读段计数 $X_{ij}/N_{ij}$）

- **索引含义**：
  - $i$：遍历基因组上所有 CpG 位点
  - $j$：遍历检测中的所有细胞（single cell）——这是单细胞方法与前面 bulk 方法（DSS、BSmooth）最本质的区别：数据结构从"位点 × 处理组"变成了"位点 × 细胞"

- **$x_{ij}$ 的三种取值**：
  - $x_{ij}=0$：亚硫酸氢盐测序检测到该位点在该细胞中**未甲基化**
  - $x_{ij}=1$：检测到**已甲基化**
  - $x_{ij}=\text{NA}$：该位点未被该细胞的任何读段覆盖，状态**缺失**——单细胞数据的核心特征是**极度稀疏**（大量 NA），这是后续所有统计建模都要面对的关键挑战

---

These values can be readily obtained from single-cell bisulfite sequencing data using tools such as Bismark.

If multiple reads from the same cell cover a position, these will typically be PCR duplicates of each other and hence agree. Of course, the two alleles of a CpG may rarely differ in their methylation status. While it is, in principle, possible that one obtains discordant reads stemming from the same position on both the paternal and the maternal chromosomes of the same cell, this is so unlikely that we can ignore such cases. Hence, whenever the methylation caller reports multiple reads covering the same position in the same cell, we set $x_{ij}$ to 0 or 1 whenever all reads agree. When there is disagreement, we put $x_{ij} = \text{NA}$ by default, or optionally follow the majority of reads whenever possible.

---

- **PCR 重复与等位基因一致性假设**：
  同一细胞若有多条读段覆盖同一位点，通常视为 PCR 扩增的重复读段，理应结果一致；理论上父源/母源两条染色体的甲基化状态可能不同，但这种情况概率极低，因此模型中忽略了等位基因特异性差异，将同一细胞同一位点视为单一状态

- **多读段冲突时的处理规则**：
  - 若同一细胞同一位点的所有读段结果一致 → 直接取该值（0 或 1）
  - 若读段结果不一致（disagreement）→ 默认设为 $\text{NA}$（保守处理，宁可丢弃也不武断判断），也可选择"多数投票"（majority rule）作为备选方案——这体现了数据预处理阶段对**测序噪声 vs 真实生物学信号**的权衡

---

For later use, we define $C$ as the set of all cells in the assay (that is, $C$ is the index set for the cell indices $j$). Moreover, we define $C_i \subset C$ as the set of all those cells $j$ that have reads covering position $i$:

$$
C_i = \{j \in C : x_{ij} \neq \text{NA}\}.
$$

Conversely, we define $G_j$ as the set of all the CpG positions $i$ covered by reads from cell $j$:

$$
G_j = \{i : x_{ij} \neq \text{NA}\}.
$$

---

- **$C$**：所有细胞索引 $j$ 构成的集合，即数据集中全部细胞的总集

- **$C_i \subset C$**：在**位点 $i$** 上有读段覆盖的所有细胞组成的子集，定义为
  $$C_i = \{j \in C : x_{ij} \neq \text{NA}\}$$
  直观理解：固定一个 CpG 位点，问"哪些细胞在这个位置上有数据"——这个集合的大小（$|C_i|$）就是该位点在单细胞层面的"**跨细胞覆盖数**"，类似 bulk 数据中的覆盖度 $N_{ij}$，但这里衡量的是"多少个细胞提供了信息"而非"多少条读段"

- **$G_j$**：在**细胞 $j$** 上有读段覆盖的所有 CpG 位点组成的集合，定义为
  $$G_j = \{i : x_{ij} \neq \text{NA}\}$$
  直观理解：固定一个细胞，问"这个细胞的测序读段覆盖了基因组上哪些位点"——由于单细胞测序通常覆盖度很低，$G_j$ 往往只占全基因组 CpG 位点的一小部分

- **$C_i$ 与 $G_j$ 的对偶关系**：
  这两个集合本质上是同一个"位点 × 细胞"稀疏矩阵（覆盖矩阵）的**行索引集**和**列索引集**——$C_i$ 是按列（位点）切片时的非缺失行集合，$G_j$ 是按行（细胞）切片时的非缺失列集合。这种"双重索引"的定义方式，是为后续需要"聚合某位点在哪些细胞间的信息"或"聚合某细胞在哪些位点间的信息"（例如估计位点均值、细胞整体甲基化水平等）打基础，对应你的框架中 **Q1 潜变量结构**——单细胞数据的潜变量组织方式（每个细胞、每个位点各自可能有潜在参数）比 bulk 数据更复杂，而 $C_i, G_j$ 正是描述这种稀疏结构的记号工具

---


## Data storage

The function 'methscan prepare' reads a set of methylation files (for example, produced by Bismark) and produces one file per chromosome. These files store the matrix $x$, where each column represents a cell and each row represents a base pair, in a space-efficient format as follows: $x$ is represented as a SciPy sparse matrix [27], encoding the actual values 0, 1 and NA as $-1$, 1 and 0, respectively. Since the vast majority of values in this matrix are missing owing to the sparsity of scBS data (and because rows for base pairs not corresponding to a CpG site contain no data), we encode missing values as zero and then store the data in compressed sparse row format. This format does not explicitly store zeroes (here, missing values) and is optimized for row-wise access, which results in substantial compression and allows fast access to the methylation status of genomic intervals. In all that follows here, any mention of $x$ will, however, always mean the encoding as $x_{ij} \in \{0, 1, \text{NA}\}$.


---

- **`methscan prepare`**：读取一批甲基化文件（如 Bismark 输出）、按染色体分别生成文件的预处理函数

- **矩阵 $x$ 的物理布局**：列代表细胞（cell）、行代表碱基对（base pair）——注意这里的行是**碱基对级别**而非仅 CpG 位点，说明存储结构是对整条基因组坐标轴的通用表示，只是非 CpG 位点行不含数据

- **SciPy 稀疏矩阵（sparse matrix）**：底层存储采用 SciPy 的稀疏矩阵结构，而非普通稠密数组——这是应对单细胞数据极度稀疏（大量 NA）的核心工程手段

- **编码映射（存储层 vs 逻辑层的区分）**：
  - 逻辑上的 $x_{ij}\in\{0,1,\text{NA}\}$（未甲基化/甲基化/缺失）
  - 实际存储时映射为 $\{-1, 1, 0\}$：$0\to-1$，$1\to1$，$\text{NA}\to0$
  - 这样做的目的是让"缺失"对应数值 **0**，从而可以利用稀疏矩阵"不存储零值"的特性来压缩数据

- **为什么用压缩稀疏行格式（CSR，compressed sparse row）**：
  1. scBS 数据本身极度稀疏（大部分位点在大部分细胞中都是 NA）
  2. 非 CpG 碱基对所在的行本来就没有数据，进一步增加稀疏性
  3. CSR 格式**不显式存储零值**（此处即缺失值），大幅节省存储空间
  4. CSR 对**按行访问**做了优化，适合快速提取"某个基因组区间内所有细胞的甲基化状态"这类操作

- **符号使用约定**：文中后续所有提到的 $x$，无论具体如何存储实现，其**逻辑含义**始终固定为 $x_{ij}\in\{0,1,\text{NA}\}$——即存储编码（$-1,1,0$）只是内部实现细节，读者理解模型时只需按原始的 0/1/NA 语义解读，不用关心底层的稀疏矩阵映射

---

## Smoothing

For each CpG position $i$, we write

$$
\bar{x}_i = \langle x_{ij} \rangle_{j\in C_i} = \frac{1}{|C_i|}\sum_{j\in C_i} x_{ij}
$$

for the average methylation at position $i$, where $\langle \cdot \rangle$ denotes averaging, and the average runs over all the cells $j \in C_i$, that is, over those cells for which methylation data are available for position $i$.

- **$\bar x_i$**：第 $i$ 个 CpG 位点的**跨细胞平均甲基化水平**（cross-cell average），是所有覆盖该位点的细胞的甲基化状态的简单算术平均

- **$C_i$**：前面定义过的集合——在位点 $i$ 上有读段覆盖的所有细胞的索引集合

- **$\bar x_i = \dfrac{1}{|C_i|}\sum_{j\in C_i} x_{ij}$**
  对固定位点 $i$，把所有"有数据"的细胞 $j\in C_i$ 上的 $x_{ij}$ 取平均——这一步本质上是把**单细胞层面的稀疏 0/1/NA 数据**，先聚合成**位点层面的一个连续值**（类似 bulk 数据中的甲基化比例），是从"单细胞分辨率"退回到"位点汇总统计量"的第一步

- **这一步对应的是"跨细胞"聚合，而非"跨位点"聚合**：先在每个位点内部，把多个细胞的信息合并成一个数（$\bar x_i$ 相当于该位点在所有测到的细胞中的平均甲基化比例）


---

We then run a kernel smoother over these per-position averages to obtain the smoothed averages $\tilde{x}_i$. Specifically, we use

$$
\tilde{x}_i = \frac{\sum_{i'} \bar{x}_{i'}\, k_h(d_{ii'})}{\sum_{i'} k_h(d_{ii'})},
$$

- **$\tilde x_i$**：对 $\bar x_i$ 序列**沿基因组坐标做核平滑**后得到的最终平滑估计值——这是第二步聚合，即"跨位点"（跨基因组邻近位置）的平滑

- **$\tilde x_i = \dfrac{\sum_{i'}\bar x_{i'}\,k_h(d_{ii'})}{\sum_{i'}k_h(d_{ii'})}$**
  这是标准的**核加权平均[Nadaraya-Waton] 型核回归**公式：分子是邻近位点 $\bar x_{i'}$ 按核权重加权求和，分母是权重归一化项，保证最终结果仍是一个合法的加权平均——与前面 BSmooth 中的局部加权 GLM 思路一致，都是"离得越近权重越大"的空间平滑，但这里更简单（直接对 $\bar x_{i'}$ 加权平均，而非局部多项式回归）

- **$d_{ii'}$**：位点 $i$ 与 $i'$ 之间的物理距离，以碱基对（bp）为单位——决定了平滑权重的空间衰减程度

- **$h$（带宽 bandwidth）**：控制平滑窗口大小的参数，默认取 $h=1{,}000$ bp——带宽越大，纳入平滑计算的邻近位点范围越广，结果越平滑（但可能过度模糊局部差异）；带宽越小则越贴近原始数据（但噪声更大）

that is, $\tilde{x}_i$ is the weighted average over the per-position averages $\bar{x}_{i'}$, taken over the CpG sites $i'$ in the neighborhood of $i$, and weighted using a smoothing kernel $k_h$ with bandwidth $h$. Here, $d_{ii'}$ is the distance between CpG positions $i$ and $i'$, measured in base pairs, $h$ is the smoothing bandwidth in base pairs (by default, $h = 1{,}000$) and $k_h$ is the tricube kernel

$$
k_h(d) =
\begin{cases}
(1 - |d/h|^3)^3 & \text{for } |d| < h \\
0 & \text{otherwise.}
\end{cases}
$$


- **$k_h$：三次核（tricube kernel）**
  $$k_h(d)=\begin{cases}(1-|d/h|^3)^3 & |d|<h\\ 0 & \text{否则}\end{cases}$$
  这是局部回归中最常用的核函数之一：
  - 当 $|d|<h$（即位点在带宽范围内）时，权重随距离增大而**平滑地衰减**到 0（在 $|d|=h$ 处权重恰好降为 0，函数本身及一阶导数都连续，没有突兀的跳变）
  - 当 $|d|\geq h$（超出带宽范围）时，权重直接为 0——即"硬截断"，超出范围的位点完全不参与平滑

- **与之前 BSmooth 中 tricube 核的呼应**：这里使用的核函数形式与 BSmooth 局部似然平滑中的加权方式完全一致（都是 tricube kernel），说明这是甲基化数据空间平滑中的**一种通用惯例**，不同工具（BSmooth、MethSCAn 等）在"如何按距离降权"这一环节上收敛到了同一套数学工具，只是在"平滑的对象是什么"（原始比例 vs. logit 变换后的量）、"是否做局部多项式回归"等细节上有所不同



---

## Methylation for an interval

Next, we discuss averaging methylation over a range of CpG sites.

Given an interval $I$ on the chromosome, we wish to quantify the average methylation $m_{Ij}$ of the CpG sites within the interval for cell $j$. If we interpret $I$ as the set of CpG positions $i$ in the interval, we may write

$$
m_{Ij} = \langle x_{ij} \rangle_{i\in I\cap C_j}.
$$

Here, the average runs over all those sites $i$ that lie within the interval $I$ and are covered by reads from cell $j$.


- **$I$**：染色体上的一个区间（interval），可以理解为一组连续 CpG 位点 $i$ 的集合——这一节要解决的问题是"如何把单细胞在一个区间内的甲基化信息汇总成一个数"

- **$G_j$**：前面定义过的集合——细胞 $j$ 有读段覆盖的所有 CpG 位点集合

- **$I\cap G_j$**：区间 $I$ 内、且被细胞 $j$ 覆盖到的那些 CpG 位点——即"在这个区间里，细胞 $j$ 到底有哪些位点是有数据的"

- **$m_{Ij} = \langle x_{ij}\rangle_{i\in I\cap G_j}$**
  细胞 $j$ 在区间 $I$ 内的**平均甲基化水平**：对区间内所有该细胞有数据的位点，把 $x_{ij}$ 取平均——这是"跨位点、单细胞内部"的聚合（区别于前一节 $\bar x_i$ 是"跨细胞、单位点"的聚合），二者是**转置关系**的两种不同汇总方向



---

If we wish to compare cells, it can be helpful to center this quantity by subtracting its average using

$$
\color{red} {z_{Ij} = m_{Ij} - \langle m_{Ij'} \rangle_{j'\in C}}
$$

- **$z_{Ij} = m_{Ij} - \langle m_{Ij'}\rangle_{j'\in C}$**
  把 $m_{Ij}$ 减去**所有细胞在该区间上的平均值**（即中心化，centering）——目的是**比较不同细胞之间**在该区间的相对甲基化水平（谁比平均水平高、谁比平均水平低），而不是看每个细胞的绝对甲基化比例

- **为什么要中心化**：不同细胞的整体甲基化背景水平可能不同（类似 bulk 数据中不同样本的"批次效应"），中心化后 $z_{Ij}$ 更能反映**细胞间的相对差异**，是后续做细胞聚类、细胞间比较分析的常用预处理手段
---

As an alternative, we suggest to consider the residuals of the individual CpG methylation values $x_{ij}$ from the smoothed average $\tilde x_i$

$$
r_{ij} = x_{ij} - \tilde{x}_i,
$$

- **$r_{ij} = x_{ij} - \tilde x_i$**
  另一种思路：不是先算区间内平均、再中心化，而是先在**单个位点**层面把观测值 $x_{ij}$ 减去该位点的**平滑后跨细胞均值** $\tilde x_i$（前一节算出来的），得到**残差**（residual）——即"这个细胞在这个位点上，相对于全体细胞在这个位点上的整体趋势，是偏高还是偏低"
---

and averaging over these, obtaining

$$
r_{Ij} = \frac{1}{|I\cap C_i|+1}\sum_{i\in I\cap C_i}\left(x_{ij}-\tilde x_i\right). \tag{1}
$$

This is a **shrunken average**, with denominator $n+1$. This extra pseudocount has the effect of shrinking the value toward the 'neutral' value 0, with the shrinkage becoming stronger if the data are 'weak' because the number $|I\cap C_i|$ of positions covered by reads from cell $j$ is low. In the extreme case of none of the reads from cell $j$ covering $I$, the sum becomes 0 and the denominator 1; that is, $r_{Ij} = 0$ in this case.



- **$r_{Ij} = \dfrac{1}{|I\cap C_i|+1}\displaystyle\sum_{i\in I\cap C_i}\left(x_{ij}-\tilde x_i\right)$**（公式 1）
  把区间 $I$ 内、细胞 $j$ 覆盖到的所有位点的残差 $r_{ij}$ 求和，再除以 $(|I\cap C_i|+1)$ 得到区间层面的汇总残差——这是本节的**核心公式**

- **分母中的 "+1"：伪计数（pseudocount）与收缩（shrinkage）**
  - 正常的平均应该除以 $n=|I\cap C_i|$（实际覆盖的位点数），但这里故意除以 $n+1$
  - 效果：当 $n$ 很大（该细胞在该区间内数据很充分）时，$n+1\approx n$，几乎不影响结果
  - 当 $n$ 很小（数据稀疏，正是单细胞数据的常态）时，$n+1$ 相对 $n$ 的比例差异更大，使得最终结果**被拉向 0（中性值）**——即"数据越少，越不敢给出极端估计，结果越保守地收缩回 0"

- **极端情况**：如果细胞 $j$ 的读段完全没有覆盖到区间 $I$ 内任何位点，则分子求和为 0（空集求和），分母因为 "+1" 变成 1，因此 $r_{Ij}=0$——这保证了"完全没有数据"时不会出现除以 0 的错误，而是优雅地退化为"中性"输出

- **与你的 Q1–Q5 框架的对应**：
  - **Q1 潜变量**：区间层面的"真实"相对甲基化水平（$z_{Ij}$ 或 $r_{Ij}$ 试图估计的目标量）
  - **Q3 先验结构**：这个 "+1" 伪计数本质上是一种**隐式的贝叶斯收缩先验**——等价于对该区间的真实值施加了一个"均值为 0、带有一个虚拟观测"的先验，是经典的"加一平滑"（add-one smoothing / Laplace smoothing）思想在这里的具体应用，只是论文没有显式写出概率模型，而是直接给出这个确定性公式
  - 这种收缩机制解决的正是单细胞数据**覆盖度极不均匀**带来的问题：覆盖度低的细胞-区间组合，其估计天然更不可靠，收缩机制自动对其做出"降权、往保守值靠拢"的处理，避免小样本导致的极端虚假信号
---


## Finding VMRs

For any interval $I$, we denote by $v_I$ the variance of its residual averages $r_{Ij}$:

$$
v_I = \frac{1}{|C_I|}\sum_{j\in C_I}\left(r_{Ij} - \langle r_{Ij'}\rangle_{j'\in C_I}\right)^2, \tag{2}
$$

where the average runs only over the set $C_I = \bigcup_{i\in I} C_i$ of those cells $j$ that have reads covering interval $I$.

To find VMRs, we define intervals $I_1, I_2, \ldots$, all of the same width, and with stepwise increasing starts, then calculate $v_1, v_2, \ldots$ for these intervals. We then mark the intervals with the 2% highest variances. We take the union of all these intervals, split the union into connected components and call each component a VMR. Putting that last step in other words: We take all the intervals with variance in the top 2 percentile, fuse intervals that overlap and call the regions thus obtained the VMRs.

- **$v_I$**：区间 $I$ 内**收缩残差 $r_{Ij}$ 在细胞间的方差**（variance），衡量的是"这个区间在不同细胞之间甲基化水平差异有多大"

- **$C_I=\bigcup_{i\in I}C_i$**：覆盖区间 $I$ 内**任意**位点的所有细胞的并集——注意这与 $I\cap C_i$（某一位点的覆盖细胞）不同，这里是把区间内所有位点各自的覆盖细胞集合取并集，范围更宽

- **公式(2)：$v_I=\dfrac{1}{|C_I|}\sum_{j\in C_I}(r_{Ij}-\langle r_{Ij'}\rangle_{j'\in C_I})^2$**
  标准的样本方差计算：先算所有细胞在该区间收缩残差的均值，再算每个细胞与该均值的偏差平方，取平均——即"细胞间变异程度"的度量

- **VMR（variably methylated region，可变甲基化区域）的定义思路**：
  1. 用固定宽度、**逐步滑动**（stepwise shifted，即滑动窗口）的方式在基因组上生成一系列重叠区间 $I_1,I_2,\dots$
  2. 对每个区间计算方差 $v_1,v_2,\dots$
  3. 取方差排名**前 2%** 的区间（即"细胞间差异最大"的那些区间）
  4. 把这些区间取并集，再按**连通性**拆分成若干连通片段——每个连通片段即为一个 VMR

- **"连通分量"的直观理解**：多个互相重叠的高方差区间会被"融合"成一个更大的区域，而不是各自独立报告——这类似图论中"取并集后找连通分量"的操作，避免因滑窗重叠而重复报告同一片区域

- **VMR 与 DMR 的区别**：VMR 只衡量"细胞间变异大不大"，**不涉及分组比较**（是一种无监督/探索性分析，找"哪里在细胞之间波动剧烈"），而 DMR 是**有监督**的（比较两组细胞之间的系统性差异）

---

## Finding DMRs

When searching for DMRs, we compare two groups of cells, whose index sets we denote by $C_A$ and $C_B$. For a given interval $I$, we calculate the mean each of the mean shrunken residuals $r_{ij}$ (equation (1)) over the cells $j$ in each of the two groups:

$$
\mu_I^g = \langle r_{Ij}\rangle_{j\in g}, \qquad g = \text{A, B}.
$$


$$
v_I^g = \frac{1}{|C_g|-1}\sum_{j\in g}\left(r_{Ij}-\mu_I^g\right)^2, \qquad g = \text{A, B}.
$$

---

- **$C_A, C_B$**：待比较的两组细胞各自的索引集合（例如两种细胞类型，或两种处理条件下的细胞）

- **$\mu_I^g=\langle r_{Ij}\rangle_{j\in g}$，$g=A,B$**：区间 $I$ 内，**组 $g$ 内所有细胞**收缩残差 $r_{Ij}$ 的均值——分别代表 A 组、B 组在该区间的"平均相对甲基化水平"

- **$v_I^g=\dfrac{1}{|C_g|-1}\sum_{j\in g}(r_{Ij}-\mu_I^g)^2$**：组 $g$ 内部、区间 $I$ 上收缩残差的**组内方差**（分母用 $|C_g|-1$，是无偏样本方差的标准写法）


---

From this, we calculate Welch's *t* statistic as usual:

$$
t_I = \frac{\mu_I^A - \mu_I^B}{\sqrt{\dfrac{v_I^A}{|C_A|} + \dfrac{v_I^B}{|C_B|}}}.
$$

- **Welch's $t$ 统计量**：
  $$t_I=\frac{\mu_I^A-\mu_I^B}{\sqrt{v_I^A/|C_A|+v_I^B/|C_B|}}$$
  分子是两组均值之差（信号），分母是两组各自的标准误合并（噪声）——这是经典**异方差双样本 t 检验**（不假设两组方差相等），与前面 BSmooth 中的 $t(l_j)=\hat\beta(l_j)/\hat\sigma(l_j)$ 思路一致，但这里是标准 Welch 形式，且比较对象是"细胞"而非"样本"

---

To find candidate DMRs, we again define overlapping and stepwise-shifted intervals $I_1, I_2, \ldots$ as for the VMRs and calculate $t$ statistics $t_1, t_2, \ldots$ for these. As before, we take the top 2 percentile of these values, fuse intervals that overlap and call the regions thus obtained candidate DMRs. We repeat the procedure for the bottom 2 percentile to get the candidate DMRs for the other sign.

- **候选 DMR 的生成流程**：与 VMR 完全类似——滑动窗口生成区间 → 计算 $t_I$ → 取绝对值最大的**前 2%**（正向差异，B 大于 A）和**后 2%**（负向差异，A 大于 B）→ 合并重叠区间为候选 DMR


Next, we need to check these candidate DMRs for statistical significance. We first remind the readers here that, as this is a **within-sample** analysis, cells, not samples, are the statistical unit. Therefore, a call as significant implies that the same DMR is likely to be called again if we repeated the analysis with another set of cells taken from the same biological sample, not that it would generalize to further samples. This fact, although often overlooked, is common to all within-sample analyses in the single-cell field, for example, also to the differential expression tests performed in scRNA-seq analyses to find marker genes that differentiate clusters.

- **"细胞是统计单位，而非样本"（cell as statistical unit）**：这是全文最重要的方法论提醒——本分析是**样本内**（within-sample）分析，即所有细胞都来自**同一个**生物学样本。因此显著性结论只能说明"换一批同一样本内的细胞重新分析，结果大概率能重现"，**不能**推广到"换一个新的生物学样本/个体也会得到同样结果"——这是单细胞领域普遍容易被忽视的推断范围限制（scRNA-seq 找 marker gene 的差异表达检验也有同样问题）

It may seem that we could use the standard procedure for the Welch $t$-test here, that is, use the Welch–Satterthwaite formula to get an approximate degree of freedom and then calculate the tail probability of the corresponding $t$ distribution. However, this is unlikely to hold for two reasons: First, the Welch–Satterthwaite degrees of freedom are only based on the number of cells per group and do not account for the fact that the read coverage might vary from cell to cell. Second, the fusing of the DMRs obtained in the scanning step to obtain fused candidate DMRs would invalidate subsequent $P$ value-based adjustment for multiple testing.

- **为什么不能直接用标准 Welch $t$-检验的理论 $p$ 值**：
  1. **Welch–Satterthwaite 自由度公式**只依赖"每组细胞数"，没有考虑不同细胞之间**读段覆盖度差异巨大**这一单细胞数据特有的问题，导致理论自由度不准确
  2. 候选 DMR 是通过"**融合重叠区间**"得到的，这个融合过程本身破坏了传统 $p$ 值多重检验校正（如 Bonferroni、BH 法）所依赖的独立性假设，使得直接套用理论分布不可靠

Therefore, we have instead implemented a permutation procedure, which works as follows: We randomly shuffle the assignment of the cells in $C_A \cup C_B$ to either of the two groups and then rerun the whole procedure, that is, the scanning step, the DMR fusing and the calculation of $t$ values from the (potentially fused) candidate DMRs. This needs to be done for a sufficiently large number of permutations. Running through the whole genome for each permutation would be too computationally expensive. Instead, we go through the genome only once, but reshuffle the group labels every **2 Mb**.

- **置换检验（permutation procedure）的具体做法**：
  1. 把 $C_A\cup C_B$ 中所有细胞的分组标签**随机打乱**
  2. 用打乱后的标签，重新跑一遍完整流程（滑窗扫描 → 计算 $t$ → 融合候选 DMR）
  3. 重复足够多次置换

- **计算效率的折衷方案**：完整地对每次置换都重新扫描全基因组太耗计算资源，因此实际做法是**只扫描一次基因组**，但每隔 2 Mb 就重新打乱一次分组标签——即在空间上分段随机化，而不是对每次置换单独跑全流程，用局部随机化近似全局置换的效果

All the $t$ values obtained from this permutation procedure are taken together to obtain an empirical null distribution. Then, we can use this null distribution to control the FDR by applying the Benjamini–Hochberg procedure in its $P$ value-free form: Let us write $T$ for the set of all $t$ values obtained from the 'real' assignment of cells to group labels and $T_0$ from the set of all $t$ values obtained from the shuffled assignments, that is, the empirical null distribution. To adjust a specific $t$ value $t_i \in T$, we calculate

- **经验零分布（empirical null distribution）**：把所有置换过程中得到的 $t$ 值汇总起来，构成一个"如果两组毫无差异，$t$ 值应该长什么样"的经验分布 $T_0$，用来代替理论上的 $t$ 分布

- **无需 $p$ 值的 Benjamini–Hochberg（BH）FDR 控制**：
  - $T$：真实分组下得到的所有 $t$ 值集合
  - $T_0$：置换（打乱分组）后得到的所有 $t$ 值集合，即经验零分布

$$
p_i^{\text{adj}} = \frac{\left|\{t_j \in T_0 : |t_j| > |t_i|\}\right| / |T_0|}{\left|\{t_j \in T : |t_j| > |t_i|\}\right| / |T|}.
$$

- **调整后的"伪 $p$ 值"公式**：
  $$p_i^{\text{adj}}=\dfrac{|\{t_j\in T_0:|t_j|>|t_i|\}|/|T_0|}{|\{t_j\in T:|t_j|>|t_i|\}|/|T|}$$
  - **分子**：在零分布 $T_0$ 中，绝对值比 $t_i$ 更极端的比例——相当于"在无差异情况下，偶然出现比 $t_i$ 更极端结果的概率"
  - **分母**：在真实数据 $T$ 中，绝对值比 $t_i$ 更极端的比例——相当于"在真实数据里，有多少候选区域比 $t_i$ 更显著"
  - **比值的直观含义**：这个比值近似给出了"如果用 $|t_i|$ 作为显著性阈值，预期的假发现率（FDR）是多少"——分子是"零假设下预期能有多少假阳性"，分母是"真实数据中实际报告了多少阳性"，二者之比即为 FDR 的经验估计

---

In words, we calculate which fraction of the null $t$ values is greater than $t$ by absolute value, and which fraction of the real $t$ values is. The ratio gives us the FDR we should expect if we used the $t$ value as threshold.

- **方法论意义（对照 Q5 框架）**：这一整套流程是 **Q5 推断验证** 层面对"检验统计量本身依赖于数据驱动的预处理步骤（滑动窗口扫描+区间融合）"这一问题的正面回应——不像 BSmooth 那样用启发式规则回避严格检验，而是**用置换检验重建整个流程的空分布**，把"扫描+融合"这个复杂的、难以用解析公式描述不确定性的过程，直接纳入了重复置换的模拟中，从而得到一个经验上更可信的 FDR 控制方法，是本文在统计严谨性上比 BSmooth 更进一步的地方

---

## Calculating cell-to-cell distances

Given a set $\mathcal{V} = \{I_1^V, I_2^V, \ldots\}$ of intervals corresponding to VMRs, we get a relative methylation fraction $r_{ij}$ for each VMR $I_i^V$ and each cell $j$ from equation (1). The matrix thus obtained can then be centered and used as input for a PCA. If we calculate the top $R$ principal components, we thus obtain for each cell $j$ an $R$-dimensional principal component vector $\mathbf{x}_j^P$.

For any two cells $j$ and $j'$, we use the Euclidean distance $\|\mathbf{x}_j^P - \mathbf{x}_{j'}^P\|$ as the measure of dissimilarity of the two cells. Thus, the matrix of PC scores can be used as input to dimension reduction methods such as t-SNE or UMAP, and to clustering methods like the Louvain or Leiden algorithm, which require such a matrix as input to the approximate nearest neighbor finding algorithm that is their first step.

## PCA with iterative imputation

Whenever a region is not covered by any read in a cell, the corresponding data entry in the input data matrix for PCA will be missing. The standard approach to calculate PCA, commonly done using the implicitly restarted Lanczos bidiagonalization algorithm [28], is not suited to deal with missing data. We circumvent this issue by simply using the PCA itself to impute the missing value in an approach that we call 'iterative PCA'.

Let us write $A$ for the matrix to which the PCA is to be applied, with the features (here, regions) represented by the matrix rows. The matrix has already been centered, that is, $\sum_i a_{ij} = 0$. To establish notation, we remind the reader that performing a PCA on $A$ means finding the singular value decomposition $A = UDR^\top$, with $D$ diagonal and $U$ and $R$ orthonormal. The PCA scores are contained in $X = UD$, the loadings in $R$. To reconstruct the input data $A$ from the PCA representation, one may use $A = XR^\top$, that is, $a_{ij} = \sum_r x_{ir} r_{jr}$, where the equation is exact if $r$ runs over all principal components and approximate if it is truncated to the leading ones.

Our iterative imputation strategy is now simply the following: We first replace all missing values in the row-centered input matrix $A$ with zeroes and perform the (truncated) PCA. Then, we use the PCA predictions for the missing values, that is, the $a_{ij} = \sum_r x_{ir} r_{jr}$, as refined stand-ins for the missing values in $A$ and run PCA once more. This can now be iterated until convergence.

We note that similar approaches have also been used elsewhere [29].






## Analysis of scBS datasets for benchmarks

To analyze scBS data from Kremer et al. [11], single-cell CpG methylation reports from all conditions were first stored with 'methscan prepare' and then smoothed with 'methscan smooth' using the default bandwidth of 1,000 bp. These data were then analyzed multiple times with different combinations of analysis methods, namely four ways to divide the genome into intervals, two ways to quantify methylation in these intervals and four approaches for dimensionality reduction.

The following four sets of genomic intervals were used: VMRs, obtained with 'methscan scan' using the current default options (bandwidth of 2,000, step size of 10, variance threshold of 0.02 and minimum cell requirement of 6); adjacent tiles of 100 kb width; promoter regions as defined by the ±2 kb domain around the TSS of coding genes; and candidate cis-regulatory regions annotated by the ENCODE consortium [12]. Methylation was quantified either by averaging binary methylation values, or by calculating the shrunken mean of the residuals as described earlier.

We used four different approaches for dimensionality reduction. Three of them involve imputation of missing values followed by PCA: The first approach, iterative PCA, was described earlier. Second, 'PCA on high-quality features' imputes missing values with the mean methylation level of a given interval, while only retaining high-quality features selected as in Luo et al. [8]: We selected tiles spanning ≥20 CpG sites and with sequencing coverage in at least 70% of cells. We then imputed missing values with the mean of each tile, centered the values and performed PCA. The third approach, 'mean-imputed PCA' is identical to the second approach but without the quality-filtering step. Lastly, we used MOFA+ with default parameters instead of PCA, which does not require imputation of missing values.

In all four cases, we reduced the dimensionality of the input data to 15 PCs or MOFA factors. In some cases, MOFA+ returned a smaller number of factors, since some of the requested 15 had zero variance. For visualization, these 15 PCs or factors were subjected to UMAP with parameters min_dist = 0.2 and init = 'spca'. To flexibly adapt to datasets of different sizes, we set

$$
n\_neighbors = \frac{\sqrt{n}}{1.5}
$$

(rounded to the nearest integer), where $n$ is the total cell number.

The same analysis was repeated for three additional scBS datasets [8,10,13] and for smaller datasets generated by randomly subsampling cells separately from these datasets.

VMRs that intersect protein-coding gene bodies, CpG islands (from the University of California Santa Cruz genome browser) or cCREs were quantified by subtracting VMRs with at least 1 bp of overlap using 'bedtools subtract -A' [30] and counting the remaining VMRs.

To test our DMR detection approach, we selected oligodendrocytes and NSCs from healthy wild-type mice of Kremer et al. [11] and ran 'methscan diff' with default parameters. For GO enrichment analysis, DMRs with adjusted $P < 0.01$ were uploaded to GREAT 4.0.4 (ref. [18]) with the association rule 'basal plus extension, 0 bp upstream, 20 kb downstream, 1 Mbp max extension, curated regulatory domains included'.


## Mean neighbor score

To assess the performance of our methods, we employed a score that quantifies how well cell types (or cell states) are separated in 15-dimensional PCA space. For data from Luo et al. [8], we used cell type labels that the authors manually curated based on CH methylation. For the multi-omic dataset, we repeated the single-cell transcriptomics analysis described in Kremer et al. [11] with two adjustments: We did not filter off-target cells such as endothelial cells, and we annotated cell types using Leiden clustering with the Seurat [6] function 'FindClusters(resolution = 0.1)'.

The score, based on the Γ score [31], varies between 0 and 1, where higher scores reflect better separation of cell types. For each cell $j$, we count how many of its $k$ nearest neighbors have been assigned to the same cell type as cell $j$. We denote this count by $a_j^k$, where we have chosen

$$
k = \frac{\sqrt{n}}{1.5}
$$

rounded. The overall score used to evaluate a given combination of methods is then simply the mean of all cell-wise scores.

## Correlation of DNA methylation and gene expression

To assess the correlation between gene expression and the methylation status of VMRs or promoters, we first detected VMRs with 'methscan scan --bandwidth 1000 --var-threshold 0.05'. We then quantified DNA methylation at VMRs and promoters with 'methscan matrix'. We defined promoters as ±2,000 bp regions centered on a gene's TSS. When multiple TSSs were annotated, we chose the TSS of the 'principal' isoform [32]. Log-normalized expression values reported by Kremer et al. [11] were then correlated with methylation of the gene's promoter or with methylation of the VMR closest to the gene body. When multiple VMRs intersected the gene body, we chose the VMR with the highest methylation variance. As a measure of methylation, we used the shrunken mean of the residuals. We omitted lowly expressed genes (with scRNA-seq counts in <5% of cells) and promoters and VMRs with low scBS coverage (in <10 cells).

MethSCAn was implemented in Python 3.8 using the packages NumPy 1.20.1, SciPy 1.6.1, numba 0.53.0 and Pandas 1.2.3. Benchmarks were performed on MethSCAn version 0.6.2 using Snakemake 7.26 and analyzed/visualized with tidyverse 1.3.1 packages.