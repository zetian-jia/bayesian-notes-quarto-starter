## DSS 原始论文（Feng, Conneely, Wu 2014, NAR）套进 Q1-Q5 框架——完整版

|   | 书中概念 | 这篇论文的具体内容 |
|---|---|---|
| **Q1** | 参数 $\theta$ | $p_{ijk}$——第 $i$ 个CpG位点、第 $j$ 组、第 $k$ 个重复的真实甲基化比例 |
| **Q2** | 似然函数 $L(\theta)=p(Y\mid\theta)$ | $X_{ijk}\mid p_{ijk},N_{ijk}\sim\text{Binomial}(N_{ijk},p_{ijk})$——$X_{ijk}$是甲基化read数，$N_{ijk}$是覆盖度，覆盖度从这里显式进入模型，这是全篇反复强调的关键点 |
| **Q3** | 先验分布 $p(\theta)$、共轭先验族 | **两层嵌套先验，对应两种不同来源的不确定性**：<br>① $p_{ijk}\sim\text{Beta}(\mu_{ij},\phi_{ij})$——捕捉重复间的生物学变异，对应词汇表里的**收缩到常数**（跨重复共享$\mu_{ij}$）<br>② $\phi_{ij}\sim\text{log-normal}(m_{0j},r_{0j}^2)$——离散度本身也有先验，$m_{0j},r_{0j}^2$ **从全基因组数据估计**（不是主观设定），这是经验贝叶斯"借力"的核心一步，同样对应**收缩到常数**这一支 |
| **Q4** | 后验分布 $p(\theta\mid Y)$ 的计算 | **两步走**：<br>① 先用**矩估计（MOM）**给每个位点粗估 $\phi_{ij}$（快速摸底，闭式解）<br>② 再用**牛顿-拉夫逊法**最大化条件后验（对 $\phi_{ij}$ 求MAP，$\mu_{ij}$ 视为已知固定值——即"轮廓后验"简化）<br>③ 基于估计出的参数，用 **Wald检验**（渐近正态近似）在每个位点做假设检验，对应"MAP + 渐近正态"这条Q4路径 |
| **Q5** | 独立于Q1-Q4的一条轴 | 三层验证，一层比一层更"刁难"模型：<br>① 在自己假设成立时的**Type I error模拟**（$\mu_{i1}=\mu_{i2}$，看假阳性率是否可控）<br>② **故意违反log-normal假设**（改用Gamma分布/真实数据经验分布生成$\phi$），测试稳健性<br>③ **换基因组物种**（拟南芥数据）测试可推广性<br>④ **连Beta分布假设本身都打破**（改用截断正态生成甲基化水平），测试模型对最核心假设错误的容忍度

### 一句话总结整篇论文在Q1-Q5里做了什么

$$
\text{Q1,Q2 照搬标准框架} \;\longrightarrow\; \text{Q3是核心创新（两层先验+经验贝叶斯）} \;\longrightarrow\; \text{Q4沿用经典工具（MOM+牛顿法+Wald检验）} \;\longrightarrow\; \text{Q5用多层压力测试验证稳健性}
$$

**这篇论文的真正贡献集中在Q3**——不是发明了新的分布或新的检验方法（Beta-Binomial、Wald检验、对数正态先验都是已有工具），而是**巧妙地把"经验贝叶斯收缩"这个技巧，应用到甲基化测序数据的离散度估计上**，解决了"单个CpG位点数据太少、方差估计不稳"这个具体痛点——这也是为什么你之前问"作者是怎么想出这么复杂的模型"时，答案是"每一块都是现成积木，创新在于怎么拼"。

### The Bayesian hierarchical model

To characterize the data, we propose the following Bayesian hierarchical model, based on the beta-binomial distribution. Notation for our model is as follows: at the $i$-th CpG site, $j$-th group and $k$-th replicate, $X_{ijk}$ is the number of reads that show methylation, $N_{ijk}$ is the total number of reads that cover this position and $p_{ijk}$ is the underlying 'true' methylation proportion. Since the process of sequencing involves the random sampling of two kinds of reads—methylated or unmethylated, $X_{ijk}\mid p_{ijk},N_{ijk}$ will follow a binomial distribution:

$$
X_{ijk}\mid p_{ijk},N_{ijk} \sim \text{Binomial}(N_{ijk},p_{ijk})
$$

- $X_{ijk}$：第 $i$ 个CpG位点、第 $j$ 组、第 $k$ 个重复中，显示甲基化的read数
- $N_{ijk}$：覆盖这个位置的总read数
- $p_{ijk}$：潜在的"真实"甲基化比例
- 测序过程涉及对两类read（甲基化/非甲基化）的随机抽样，因此给定 $p_{ijk}$ 和 $N_{ijk}$，$X_{ijk}$ 服从二项分布

给定 $p_{ijk},N_{ijk}$ ，那么 $X_{ijk}$ 的分布服从二项分布

Since the true methylation proportions among replicates can be anywhere between 0 and 1, we **assume** that the proportions for each CpG site within each group of replicates **follow a beta distribution**. The beta distribution has long been a natural choice to model binomial proportions as it is a conjugate distribution of the binomial distribution and is the most flexible distribution with a support interval of $[0,1]$.

$$
p_{ijk} \sim \text{Beta}(\mu_{ij},\phi_{ij})
$$

Here the beta distribution is parameterized by mean (denoted by $\mu_{ij}$) and dispersion (denoted by $\phi_{ij}$). Compared with the traditional parameterization of the Beta$(\alpha,\beta)$ distribution, the parameters have the following relationship:

$$
\mu=\frac{\alpha}{\alpha+\beta},\quad \phi=\frac{1}{\alpha+\beta+1}
$$

- $p_{ijk}$：位点 $i$、组 $j$、重复 $k$ 的真实甲基化比例（前面公式里的那个未知参数）
- **假设它服从 Beta 分布**——因为真实比例只能落在 $[0,1]$ 区间内，而 Beta 分布正是这个区间上最灵活、又是二项分布共轭先验的分布
- $\mu_{ij}$：这里的 Beta 分布不用传统的 $(\alpha,\beta)$ 参数化，改用**均值**参数化——直接表示"这组重复的平均甲基化比例"，更直观
- $\phi_{ij}$：**离散度**参数——衡量"这组重复之间，甲基化比例的分散程度"，$\phi$ 越大，同一组内不同重复之间的差异越大
- $\mu=\dfrac{\alpha}{\alpha+\beta}$：均值参数化和传统 $(\alpha,\beta)$ 参数化之间的换算关系——两者本质是同一个分布，只是换了个更贴合应用场景的参数写法
- $\phi=\dfrac{1}{\alpha+\beta+1}$：$\alpha+\beta$ 越大，$\phi$ 越小——$\alpha+\beta$ 相当于"先验里的虚拟样本量"，虚拟样本量越大代表先验越确信，对应离散度越小（组内重复越一致）

In this hierarchical model, the biological variation among replicates is captured by the beta distribution and the variation due to the random sampling of DNA segments during sequencing is captured by the binomial distribution. The dispersion parameter $\phi_{ij}$ captures the variation of a CpG site's methylation proportion relative to the group mean. We allow each CpG site within a single condition (e.g. within cases, or controls) to have its own dispersion. This is a flexible assumption because it allows either different or common dispersions for both conditions; however, our software also includes an option to assume a common dispersion for cases and controls.

To combine information across all CpG sites, based on the observed distribution of dispersion from a publicly available RRBS dataset on mouse embryogenesis (21), we assumed the following prior on $\phi_{ij}$:

$$
\phi_{ij} \sim \text{log-normal}(m_{0j}, r_{0j}^2)
$$

where $m_{0j}$ and $r_{0j}^2$ are mean and variance parameters that can be estimated from the data. For each CpG site in this dataset, we applied a method of moments (MOM) estimator to estimate the dispersion parameters. As shown in Figure 1, the genome-wide distribution of logarithm dispersion parameter estimates is approximately Gaussian with mean $=-3.39$ and SD $=1.08$, suggesting that the dispersion parameters can be well-described by a log-normal distribution. However, simulations using dispersions from different distributions also show that our proposed method is robust to violations of this log-normal assumption (Supplementary Figure S1).

- **这个层次模型的两层分别捕捉了什么**：Beta分布负责捕捉"生物学重复间的变异"，Binomial分布负责捕捉"测序时对DNA片段随机抽样带来的变异"——两种不同来源的不确定性，被分配到了模型的两个不同层级
- $\phi_{ij}$：离散度参数，衡量"某个CpG位点的甲基化比例，相对于所在组均值的分散程度"
- **每个CpG位点在单个条件下（比如病例组或对照组内）都有自己独立的离散度**——这是一个灵活的设定，既允许两种条件下离散度不同，也允许设成相同（软件里提供了这个选项）
- **为了跨CpG位点整合信息**，作者给 $\phi_{ij}$ 本身也定义了一个先验分布——**这正是"经验贝叶斯"的关键一步**：不让每个位点的离散度完全独立无关，而是让它们共享一个更高层的先验，借力估计
- $\phi_{ij}\sim\text{log-normal}(m_{0j},r_{0j}^2)$：离散度的先验取对数正态分布——因为离散度 $\phi$ 本身必须是正数，取对数之后变成一个可以是任意实数的量，正好可以用你已经很熟悉的正态分布来建模
- $m_{0j}$ 和 $r_{0j}^2$：对数正态分布的均值和方差参数，**从数据本身估计出来**（不是主观设定）——这是经验贝叶斯"用数据本身确定先验超参数"这个核心思想的具体体现
- **矩估计（MOM）**：一种比最大似然更简单的参数估计方法，用样本的矩（均值、方差）直接反推分布参数，不需要迭代优化，计算成本很低——这也是为什么它适合"先给每个位点粗略估计一下离散度，再汇总看整体分布长什么样"这种探索性的第一步
- **全基因组尺度上，$\log(\phi)$ 的分布近似高斯**（均值$-3.39$，标准差$1.08$）——这是作者从真实数据里观察到的经验规律，支撑了"用对数正态作为先验"这个选择是合理的，不是凭空拍脑袋决定的
- **稳健性说明**：即使这个对数正态假设不完全成立，模拟显示方法依然可靠——这是一个负责任的补充说明，告诉读者这个先验假设的选择不是模型能否正常工作的"生死线"

### Parameter estimation

To estimate the parameters of the prior distribution in a general setting, we first use the MOM to estimate the dispersion parameters for all CpG sites, and then estimate $m_{0j}$ and $r_{0j}^2$ as the mean and variance of the logarithm of the dispersion estimates. The mean methylation levels are estimated as $\hat\mu_{ij}=\dfrac{\sum_k X_{ijk}}{\sum_k N_{ijk}}$. Under the hierarchical model, the conditional posterior distribution of $\phi_{ij}$ satisfies:

$$
\log\left(p(\phi_{ij}\mid x_{ijk},N_{ijk},\mu_{ij})\right) \propto \sum_k \varphi\left(x_{ijk}+(\phi_{ij}^{-1}-1)\mu_{ij}\right)
$$

$$
+\sum_k \varphi\left(N_{ijk}-x_{ijk}+(\phi_{ij}^{-1}-1)(1-\mu_{ij})\right)
$$

$$
-\sum_k \varphi\left(N_{ijk}+(\phi_{ij}^{-1}-1)\right) - n\varphi\left((\phi_{ij}^{-1}-1)\mu_{ij}\right)
$$

$$
-n\varphi\left((\phi_{ij}^{-1}-1)(1-\mu_{ij})\right)+n\varphi\left(\phi_{ij}^{-1}-1\right)
$$

$$
-\log(\phi_{ij})-\log(r_{0j})-\frac{(\log(\phi_{ij})-m_{0j})^2}{2r_{0j}^2}
$$

A point estimate of $\phi_{ij}$ can be obtained by maximizing this conditional posterior likelihood. In practice, we use the Newton–Raphson method after plugging in the estimates of $m_{0j}, r_{0j}^2$ and $\mu_{ij}$. Because we estimate $m_{0j}$ and $r_{0j}^2$ from the data, the estimated $\phi_{ij}$ is therefore an empirical Bayes estimate, which shrinks toward the common prior mean. Also notable is that the last line of the above equation includes the penalty function $-\log(\phi_{ij})-\log(r_{0j})-\dfrac{(\log(\phi_{ij})-m_{0j})^2}{2r_{0j}^2}$, which will penalize extremely large $\phi_{ij}$ in our estimation.

**关键参数中文解释**

- $\varphi(\cdot)$：这里是**对数Gamma函数**（$\log\Gamma(\cdot)$）的简写记号——注意跟 $\phi_{ij}$（离散度参数）**不是同一个符号**，只是印刷体长得像，容易搞混；$\varphi$ 出现是因为 Beta-Binomial 分布的密度公式本身包含 Gamma 函数的比值，取对数后自然变成一串 $\log\Gamma$ 相加减
- $\hat\mu_{ij}=\dfrac{\sum_k X_{ijk}}{\sum_k N_{ijk}}$：均值的估计方式——把这一组内所有重复的甲基化read数直接加总，除以所有重复的总覆盖度（不是先算每个重复的比例再取平均，是先合并计数再算比例），这样覆盖度高的重复自然获得更大的权重
- $\log\left(p(\phi_{ij}\mid\ldots)\right)\propto\ldots$：这一长串公式就是给定观测数据后，$\phi_{ij}$ 的条件后验对数密度——前面几行是似然部分（Beta-Binomial分布对数似然展开后的样子），最后一行是先验部分
- **最后一行的惩罚函数** $-\log(\phi_{ij})-\log(r_{0j})-\dfrac{(\log(\phi_{ij})-m_{0j})^2}{2r_{0j}^2}$：这正是对数正态先验 $\phi_{ij}\sim\text{log-normal}(m_{0j},r_{0j}^2)$ 取对数后的样子
- **牛顿-拉夫逊法（Newton-Raphson）**：一种数值优化算法，用二阶导数信息加速找到函数最大值——这里用它来找让"条件后验密度最大"的 $\phi_{ij}$，对应Q4的一种具体实现（MAP估计的数值求解方式）
- **"这是一个经验贝叶斯估计，会向共同先验均值收缩"**：因为 $m_{0j},r_{0j}^2$ 是从全基因组数据估出来的（不是主观拍脑袋），所以每个位点最终估计出的 $\phi_{ij}$，都会被这个"全基因组共享的先验"往中心拉一把——这正是经验贝叶斯"借力"的具体数学体现，也是"伪计数收缩"思想的正式实现
- **惩罚极端大的 $\phi_{ij}$**：先验里那个 $-\log(\phi_{ij})$ 项，会让优化过程天然不喜欢 $\phi_{ij}$ 取太大的值——防止某些位点因数据碰巧很少、很嘈杂，就估计出一个夸张到不合理的离散度，起到"稳健化"的作用

- **轮廓似然/轮廓后验**：这种"先固定一部分参数，只对剩下的参数做完整贝叶斯"的做法，叫profile likelihood/posterior（轮廓似然/轮廓后验）——不是偷懒，是工程上的权衡。$\mu_{ij}$ 也是模型的一个未知参数，但在这个公式里被当成已知的、固定的数摆在了条件里。

- 这是一个条件后验，只对 $\phi_{ij}$ 一个参数，不是联合后验对 $\phi_{ij}, \mu_{ij}$

#### 条件后验公式是如何推导出来的

1. 第一步：还是"后验 $\propto$ 似然 × 先验"

$$
p(\phi_{ij}\mid x_{ijk},N_{ijk},\mu_{ij}) \propto \underbrace{p(x_{ijk}\mid N_{ijk},\mu_{ij},\phi_{ij})}_{\text{似然：Beta-Binomial}} \times \underbrace{p(\phi_{ij})}_{\text{先验：对数正态}}
$$

取对数，乘法变加法：

$$
\log p(\phi_{ij}\mid\ldots) = \log(\text{似然}) + \log(\text{先验}) + \text{常数}
$$

论文公式的前四行是"$\log(\text{似然})$展开后的样子"，最后一行是"$\log(\text{先验})$"。

2. 第二步：似然部分从哪来——Beta-Binomial 的标准公式

$$
p(X_{ijk}\mid N_{ijk},\alpha,\beta) = \binom{N_{ijk}}{X_{ijk}}\frac{B(X_{ijk}+\alpha,\ N_{ijk}-X_{ijk}+\beta)}{B(\alpha,\beta)}
$$

$B(a,b)=\dfrac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$ 是Beta函数。取对数：

$$
\log p(X_{ijk}\mid\ldots) = \log\binom{N_{ijk}}{X_{ijk}} + \log\Gamma(X_{ijk}+\alpha)+\log\Gamma(N_{ijk}-X_{ijk}+\beta)-\log\Gamma(N_{ijk}+\alpha+\beta) - \left[\log\Gamma(\alpha)+\log\Gamma(\beta)-\log\Gamma(\alpha+\beta)\right]
$$

3. 第三步：$\varphi(\cdot)$ 就是 $\log\Gamma(\cdot)$ 的简写

$$
\varphi(\cdot) \equiv \log\Gamma(\cdot)
$$

代入后：

$$
\varphi(X_{ijk}+\alpha)+\varphi(N_{ijk}-X_{ijk}+\beta)-\varphi(N_{ijk}+\alpha+\beta)-\varphi(\alpha)-\varphi(\beta)+\varphi(\alpha+\beta)
$$

4. 第四步：把 $\alpha,\beta$ 换成 $\mu_{ij},\phi_{ij}$

$$
\alpha=(\phi_{ij}^{-1}-1)\mu_{ij},\qquad \beta=(\phi_{ij}^{-1}-1)(1-\mu_{ij}),\qquad \alpha+\beta=\phi_{ij}^{-1}-1
$$

逐项代入替换，精确得到论文公式的六项：

$$
\varphi\left(X_{ijk}+(\phi_{ij}^{-1}-1)\mu_{ij}\right) + \varphi\left(N_{ijk}-X_{ijk}+(\phi_{ij}^{-1}-1)(1-\mu_{ij})\right) - \varphi\left(N_{ijk}+(\phi_{ij}^{-1}-1)\right)
$$
$$
- \varphi\left((\phi_{ij}^{-1}-1)\mu_{ij}\right) - \varphi\left((\phi_{ij}^{-1}-1)(1-\mu_{ij})\right) + \varphi\left(\phi_{ij}^{-1}-1\right)
$$

组合数项 $\log\binom{N}{X}$ 不含 $\phi_{ij}$，是常数，被 $\propto$ 直接吸收扔掉。

5. 第五步：对所有重复 $k$ 求和

一组内有 $n$ 个重复，每个重复贡献一份似然，独立观测取对数后连加 $\sum_k$——前三项每个重复的 $X_{ijk},N_{ijk}$ 不同要单独求和；后三项不含 $k$（只跟共享的 $\mu_{ij},\phi_{ij}$ 有关），$n$ 个重复贡献同一个值 $n$ 次，直接写成 $n\times(\cdot)$。

6. 第六步：最后一行是对数正态先验的对数密度

$$
\phi_{ij}\sim\text{log-normal}(m_{0j},r_{0j}^2) \;\Rightarrow\; \log p(\phi_{ij}) = -\log(\phi_{ij})-\log(r_{0j})-\frac{(\log\phi_{ij}-m_{0j})^2}{2r_{0j}^2} + \text{常数}
$$

标准正态分布密度公式取对数后的形式，只是套用到了 $\log\phi_{ij}$ 上。

7. 一句话总结

$$
\text{整个公式} = \underbrace{\sum_k\log(\text{Beta-Binomial似然})}_{\text{代入}\alpha,\beta\to\mu,\phi\text{后展开}} + \underbrace{\log(\text{对数正态先验})}_{\text{标准正态密度公式套用到}\log\phi\text{上}}
$$

没有新技巧——纯粹是"后验∝似然×先验"取对数后逐项代数展开的结果，$\alpha,\beta\to\mu,\phi$ 的变量代换让每一项变复杂了，但可以在心里折叠回 $\log(\text{似然})+\log(\text{先验})$ 这一句话。



### Statistical test procedure

After estimating the parameters for each group as described above, hypothesis tests can be performed at each CpG site to compare mean methylation levels between two groups, e.g. test $H_0:\mu_{i1}=\mu_{i2}$. We propose to use a **Wald test**. The variance of $\hat\mu_{ij}$ is derived as follows. First, the variance of $X_{ijk}$ is (based on beta-binomial distribution):

$$
\text{var}(X_{ijk}) = N_{ijk}\mu_{ij}(1-\mu_{ij})\left[1+(N_{ijk}-1)\phi_{ij}\right]
$$

So,

$$
\text{var}(\hat\mu_{ij}) = \text{var}\left(\frac{\sum_k X_{ijk}}{\sum_k N_{ijk}}\right)
$$

$$
= \left(\frac{1}{\sum_k N_{ijk}}\right)^2\sum_k\left\{N_{ijk}\mu_{ij}(1-\mu_{ij})\left[1+(N_{ijk}-1)\phi_{ij}\right]\right\} \tag{1.1}
$$

The estimated variance of $\mu_{ij}$ can be obtained by plugging in estimated values of $\mu_{ij}$ and $\phi_{ij}$ to Equation (1.1). For two-group comparisons, a Wald test of the $i$-th CpG site is:

$$
t_i = \frac{\hat\mu_{i1}-\hat\mu_{i2}}{\sqrt{\hat{\text{var}}_{i1}+\hat{\text{var}}_{i2}}} \tag{1.2}
$$

where $\hat{\text{var}}_{ij}$ ($j=1,2$) is the estimated variance for group 1 or 2. It is not trivial to derive the null distribution of the test statistics. However, based on simulation results which suggest that the empirical null distribution of the test statistics is approximately normal (Figure 5), it is possible to calculate approximate $P$-values based on the normal distribution.

**关键参数中文解释**

- $\text{var}(X_{ijk})=N_{ijk}\mu_{ij}(1-\mu_{ij})\left[1+(N_{ijk}-1)\phi_{ij}\right]$：这是 Beta-Binomial 分布的方差公式——跟普通二项分布方差 $N\mu(1-\mu)$ 相比，多了一个 $\left[1+(N-1)\phi\right]$ 因子，这个因子 $>1$（因为 $\phi>0$），说明 Beta-Binomial 的方差永远比普通二项分布更大——这正是"过度离散"（overdispersion）的数学体现，$\phi_{ij}$ 就是控制这个"多出来多少方差"的旋钮

- $\hat\mu_{ij}=\dfrac{\sum_k X_{ijk}}{\sum_k N_{ijk}}$：均值的点估计——"合并计数再算比例"的估计方式

- **公式(1.1)**：把 $\hat\mu_{ij}$ 写成"加权求和"的形式后，用方差的线性变换性质（$\text{Var}(aX)=a^2\text{Var}(X)$，各重复间独立可以拆开求和）推出来的——本质上是把每个重复各自的方差（用上面那条Beta-Binomial方差公式）加权汇总起来

- $t_i=\dfrac{\hat\mu_{i1}-\hat\mu_{i2}}{\sqrt{\hat{\text{var}}_{i1}+\hat{\text{var}}_{i2}}}$：**Wald检验统计量**——分子是"两组估计值的差"，分母是"这个差值的标准误"，跟 $t$ 检验公式 $t=\dfrac{\bar{x}_1-\bar{x}_2}{\text{SE}}$ 结构完全一样，区别只在于这里的方差用的是Beta-Binomial模型算出来的、已经把覆盖度信息考虑进去的方差，而不是普通$t$检验那种忽略覆盖度的样本方差

- **"渐近正态"**：作者没有严格推导出 $t_i$ 精确服从什么分布，而是通过模拟发现它近似服从正态分布，进而用正态分布去近似计算p值——这正是拉普拉斯近似思路的实际应用，也呼应了 Q4 计算策略"渐近正态近似"这一分支

### Defining DMRs

Based on the calculated $P$-value at each CpG site, we implemented a simple procedure in DSS for calling DMRs based on the approximate Wald test $P$-values described above. To call DMRs, the user needs to **specify a $P$-value threshold and a few other parameters**. Called DMRs must exceed a minimum length (100 bp by default) and cover more than a minimum number of CpG sites (three by default), and the percentage of CpG sites in the DMR with $P$-values less than the threshold must exceed a user-specified value (80% by default). Regions satisfying the above criteria will be reported as DMRs. Note that in this procedure, the correlation of the $P$-values for proximal sites is not considered; incorporation of this information into the DMR detection method is a direction of future research.

### Simulations

We used simulation data to test the proposed method and compare the results with existing methods. Simulations are based on mouse embryogenesis data (Gene Expression Omnibus accession GSE34864) from RRBS experiments in a study on mouse embryogenesis (21). For each simulation, we simulated 20,000 CpG sites for replicates from two groups, where the number of replicates per group is taken as 2, 3 or 5. We first computed $\mu_{ij}$ for each of the CpG sites based on the average methylation proportions from a set of 20,000 contiguous CpG sites in the mouse embryogenesis data. For Type I error simulations, we let $\mu_{i1}=\mu_{i2}$ for all CpG sites; for simulations that included DML, we allowed $\mu_{ij}$ to vary between groups for a randomly selected 5% of CpG sites in each simulation. We next simulated the dispersion parameter $\phi_{ij}$ for each CpG site from a log-normal distribution with parameters estimated from the data (mean $=-3.39$, var $=1.08$) as described above. To check the robustness of our model to departures from this distributional assumption, we also performed simulations with $\phi_{ij}$ drawn from a Gamma distribution (with parameters estimated from the data, shape $=1.5$, scale $=0.02$) and empirically sampled from real data estimates. For coverage, we simulated coverage depth ($N_{ijk}$) for each CpG site and replicate by sampling the coverage depth from real RRBS data. Finally, for each replicate at each CpG site, we then used $\mu_i,\phi_{ij}$ and $N_{ijk}$ to simulate methylated counts for each CpG site based on the beta-binomial distribution.

For additional simulations based on a different genome, parameters were estimated from a publicly available whole genome Arabidopsis dataset (Gene Expression Omnibus accession GSE38991). In this situation, a similar approach was used for generating simulation data. Again, we used log-normal (mean $=-4.3$, var $=1.7$), Gamma (shape $=0.43$, scale $=0.06$) and empirical distributions to generate dispersion parameter $\phi$.

We also performed additional simulations using a distribution other than the assumed beta distribution to generate methylation levels. In these simulations, methylation levels of biological replicates within each treatment group and CpG site were generated from a truncated normal distribution. Each CpG site and group had its own truncated normal distribution, with the parameters estimated from the mouse embryogenesis data. Since methylation levels range from 0 to 1, the boundaries of each truncated normal distribution were set to be 0 and 1.