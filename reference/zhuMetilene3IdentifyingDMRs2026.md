> 原文：*Zhu, Z., Bernhart, S.H., Jühling, F., Kretzmer, H. & Hoffmann, S. metilene³: identifying DMRs across multiple conditions with auto-classification. Nature Communications (2026). https://doi.org/10.1038/s41467-026-74931-y*
> 本文为开放获取文章（CC BY 4.0），以下为全文中文翻译，供内部研究使用。参考文献列表保留原文（英文），未逐条翻译。


---

| | metilene³ 对应 | 一句话记忆点 |
|---|---|---|
| **Q1 潜变量** | 每个样本在每个基因组区段的真实甲基化状态（低/高/中间）+ 样本真实分组身份 | 看不见的是"谁跟谁是一类" |
| **Q2 观测模型** | βᵢ,ₑ甲基化值；用Z分数/2D-KS/Kruskal-Wallis**检验统计量代替显式似然** | 噪声大→统计量自然不显著，没有参数化噪声模型 |
| **Q3 先验结构** | 同时用了5种：**马尔可夫**(基因组坐标光滑)、**收缩到常数**(段内视为常数)、**分组**(G/H/I局部分类)、**稀疏**(DMR全局稀疏)、**低秩/层次共享**(DMTree跨DMR汇总模式) | 这是全文最精华的一格，五种先验拼出了整个算法 |
| **Q4 计算策略** | 贪心递归二分(CBS)+确定性阈值分类+贪心层次分裂(DMTree)，**无MCMC/VI，一次性可并行** | 只决定快不快：牺牲了后验，换来笔记本电脑几分钟出结果 |
| **Q5 独立轴：选择性推断** | **先搜索最优切分/分组，再在同一份数据上检验**——经典double dipping；无监督模式下"分组"和"检验"是同一批数据的双重循环 | 报告的p值/显著性是**乐观估计**，metilene³更适合当"假说生成器"而非终审结论(论文自己也用"hypothesis-generating"措辞) |

---

# 创新性（5点精简版）

1. **顺序倒转**：不再"先分组、后测差异"，而是"先局部测差异、再反推分组"——能抓住只在小片区域出现的稀有/过渡信号。
2. **DMTree即解释**：每次分裂天然对应一批具体DMR，聚类结果自带"可查到基因座"的可解释性，不是先聚类后解释的黑箱。
3. **靠重复出现来抗噪**：不显式建模批次效应，而是要求分组模式在大量独立DMR上反复出现(ΣWᵢ>Wmin)，隐式实现跨区段一致性检验。
4. **工程上的完全并行**：区段独立可并行+预分段，把平方级增长的比较量压到笔记本电脑可用的量级。
5. **软肋**：分组发现和显著性检验用同一份数据，存在选择效应——结果适合排序筛选，不适合当无偏统计结论。

> **一句话记忆**：metilene³ = 局部测差异 → 跨区段聚合成树 → 树反哺再检测 的**自举闭环**；统计上朴素(非贝叶斯)，工程上高效，代价是Q5选择效应没有被显式修正。

# metilene³：跨多条件识别差异甲基化区域（DMR），并支持样本自动分类

**收稿日期：** 2025年10月18日　**接受日期：** 2026年6月17日

**作者：** Zhihan Zhu¹˒²、Stephan H. Bernhart³、Frank Jühling⁴、Helene Kretzmer¹˒⁵、Steve Hoffmann²˒⁶

¹马克斯·普朗克分子遗传学研究所（德国柏林）；²柏林自由大学数学与计算机科学系；³莱比锡大学代谢研究中心（LeiCeM）；⁴斯特拉斯堡大学/法国国家健康与医学研究院转化医学研究所；⁵波茨坦大学哈索·普拉特纳数字工程研究所；⁶耶拿弗里茨·利普曼研究所（莱布尼茨老龄化研究所）Hoffmann课题组。通讯作者：helene.kretzmer@hpi.de；steve.hoffmann@leibniz-fli.de

## 摘要

DNA甲基化是众多物种中一种关键的表观遗传标记，识别差异甲基化区域（DMR）对理解基因组调控至关重要。现有的大多数DMR检测方法都要求预先设定样本所属条件，这限制了在群体身份未知或不确定情形下（临床场景中很常见）对新表观遗传模式的发现。此外，能够进行多条件比较的方法也十分有限。为填补这一重要空白，我们提出了metilene³——一种可在监督与无监督两种模式下运行的快速多条件DMR检测方法，既可使用用户提供的标签，也可对未标注样本进行自主聚类。metilene³通过对多个成对甲基化差异信号进行基因组分段，实现了样本分类以及以DMR为锚点的表观遗传关系推断。基于模拟数据与多种人类数据集，我们证明metilene³能够准确检测DMR、稳健地对样本聚类，并有潜力揭示新的调控元件与样本分层。具体而言，在胰腺组织数据集中，metilene³识别出富集了参与胰腺癌发生的关键转录因子的DMR，提示了一个可能改变的NFKB–NFAT调控程序。总之，metilene³为探索异质性甲基化组、发现复杂生物学及临床数据集中的表观遗传模式，提供了一个快速、可解释的分析框架。

## 引言

胞嘧啶5号碳原子的甲基化（5mC）是一种可遗传的表观遗传标记，在细胞分裂过程中得以维持，并在不依赖底层DNA序列的情况下，对塑造细胞身份和发育状态起着关键作用。在哺乳动物中，胞嘧啶甲基化主要发生于CpG位点，是发育、衰老与疾病等多种调控过程的重要指标[1]。例如在癌症中，DNA甲基化模式常常发生改变，包括全局性低甲基化和局部启动子高甲基化[2]。这些异常模式可能参与抑癌基因沉默或致癌通路激活，从而促进疾病的发生和进展，例如胰腺癌中CCND2（细胞周期蛋白D2）的甲基化[3]，以及IDH突变型胶质瘤中的表观遗传改变[4]。因此，研究DNA甲基化差异能够为理解生物学过程提供重要线索，并为各类疾病中的生物标志物发现提供机会。此类研究的一个关键环节，就是识别差异甲基化的CpG位点及差异甲基化区域（DMR）。然而，包括metilene（第5号参考文献）在内的现有DMR调用工具和分析流程，通常都要求样本被划分类别，即被分配到特定条件[5–11]。无论是监督式还是无监督式，在多个（可能未知的）分组中检测DMR，仍然是一项重大挑战。此外，随着大规模甲基化研究数量的增长，系统性地识别具有生物学意义的DMR正变得日益复杂，这凸显了开发更灵活、更高效分析方法的必要性。

诸如SMART2[12]和wgbstools[13]之类的方法已被开发用于支持多条件比较。简而言之，SMART2利用香农熵来识别对特定细胞类型具有特异性的CpG位点和区域[12]；wgbstools通过最小化相邻CpG的甲基化方差对基因组进行分段，随后通过"一对多"比较各分段的平均甲基化水平来识别群组特异性模式[13]。近期，两种不依赖分组信息的差异甲基化检测方法相继被提出：其一为Aclust2.0[14]，基于Illumina甲基化芯片；其二为Methylscore[15]，专为植物设计，基于隐马尔可夫模型提供无监督DMR检测及局部聚类，但不进行全局分组。然而，许多研究需要一个更灵活的分析框架——既能在不完全依赖预定义条件的情况下识别甲基化差异，又能支持对样本的下游分类。这对于识别新的亚组（例如疾病亚型）以改善诊断、预后预测甚至治疗尤为重要[16,17]。迄今为止，用于此类任务的大多数聚类方法都是先选取高变异性的CpG位点[18]，再进行K-均值或层次聚类。这种预先选择过程可能会大幅限制生物学可解释性，例如将CpG位点从其CpG邻域中孤立出来，并剔除那些变化虽小但一致的CpG位点。

metilene³填补了这一方法学空白。在其前身核心算法基础上，metilene³采用新算法实现了无监督、无标签的DMR检测与样本分类、监督式多组比较，以及对所推断DMR的系统性分析。metilene³实现了一种多样本版本的环形二元分割（circular binary segmentation, CBS）算法，可在各个基因组位点独立地对样本进行分类（局部分类）。DMR及其对应的局部分类结果被用于对样本进行分裂式聚类，从而构建出反映输入样本间表观遗传关系的二叉连接树（差异甲基化树，DMTree）。这一步骤对于将样本相似性中反复出现的模式，与噪声或混杂因素驱动的分组区分开来尤为重要。因此，metilene³通过聚焦于稳健模式，自动对无监督样本聚类结果进行自我校正。基于模拟数据和生物学数据，我们证明metilene³具有敏感、准确、快速的特点，即使在多条件场景下，也能在标准笔记本电脑上完成DMR检测。通过对公开数据集的重新分析，我们展示了metilene³在揭示细胞分化和癌症生物学相关洞见方面的能力，凸显了利用先进DNA甲基化分析工具重新审视既有数据的价值。

## 结果

### 算法

**一个用于跨多条件研究甲基化组的工具。** metilene³引入了两项关键进展：（1）除传统的双组分析外，还支持监督式多组比较；（2）具备无监督模式，可在无需用户提供样本分类信息的情况下检测DMR。为实现这一目标，metilene³在DMR检测中增加了局部分类步骤，从而基于DMR间反复出现的"样本—分类"模式构建二叉聚类树（DMTree）。该方法旨在通过利用具有显著甲基化差异的、反复出现的分组组合，降低对噪声和批次效应的敏感性，从而有助于发现稳健的、可能是新颖的样本亚组，并支持对样本或亚组特异性改变的探索——这些结果随后被用于监督式DMR检测。关键参数及其默认值在方法部分描述。此外，metilene³还提供汇总报告，支持数据可视化及下游处理，包括图形生成、将DMR注释到最近基因，以及基因集富集分析（图1a）。

对于两组或多组（在无监督模式下则是多个样本）之间的DMR检测，metilene³采用最初提出的环形二元分割（CBS）算法[5]，对平均甲基化差异信号（MDS）进行分割（图1b）。由于无监督模式代表更一般的情形，我们首先在此背景下描述metilene³的方法。

**无监督模式下的环形二元分割与局部样本分类。** 对于某一区段 [s, t] 内所有样本对，在满足 s ≤ u < v ≤ t 的基因组坐标 [u, v] 条件下，对MDS的Z分数进行最大化。使Z分数最大化的样本对 (O, R) 及区段 [a, b] 将被进一步处理：在 [a, b] 内，其中一个样本相对另一个呈高甲基化（或低甲基化）状态。因此，条件 O 和 R 被分别选作集合 G 和 H 的首个元素。接下来，根据在 [a, b] 内与 O、R 的最小甲基化距离，将其余所有样本分配至 G 或 H。若无法进行明确分类，样本则被归入一个模糊的、中间态的分组（图1b）——这提示该样本可能处于过渡状态，其潜在细胞群体在该区域呈现异质性甲基化模式，导致其与低甲基化组和高甲基化组之间均只表现出细微差异。我们将此过程称为"局部样本分类"，因为样本仅根据其在 [a, b] 内的甲基化值被分配至三个集合之一。随后，metilene³对G组和H组在 [a, b] 内的所有甲基化值进行二维柯尔莫哥洛夫-斯米尔诺夫（KS）检验。基于检验结果及其他附加标准（如CpG数目），最终判定该区段是否被接受为DMR。若被接受，递归过程将在区段 [s, a) 和 (b, t] 上继续；否则，区段 [s, t] 的递归终止。

**DMTree的构建。** 在CBS完成后，所识别出的DMR及其对应的局部样本分类结果被用于构建DMTree（图1c）。从概念上讲，对于给定的一组样本，metilene³会识别在多个DMR之间反复出现的局部分类模式，即样本被分配到G、H、I三个集合中的相似模式。为确定树构建每一步中最具影响力的分类，metilene³计算每种模式的ΣWᵢ，即该模式对应的所有DMR的平均MDS值之和，称为总甲基化差异。

DMTree的构建从考虑全部样本开始。具有最大ΣWᵢ值的分类模式被置于树根，样本随之被划分为低甲基化组和高甲基化组，构成DMTree的最初两个分支。为避免分裂结果被组大小所主导，中间组I中的样本被临时归入较小的一组（详见方法部分）。对于每个分支，ΣWᵢ值会基于该分支所含样本子集重新计算，随后选取ΣWᵢ最高的模式来决定下一次分裂。此过程反复进行，直至仅剩一个样本，即到达叶节点。若某分支的样本总ΣWᵢ值低于Wmin，这些样本将被合并归入metilene³所报告的某一聚类簇。在最终聚类完成后，metilene³会以这些聚类簇作为分组，再次运行DMR检测（此时为监督模式），以识别更多的DMR（图1a）。

**监督式多组DMR检测。** 在监督模式下，对于某区段 [s, t]，metilene³同样寻求在满足 s ≤ u < v ≤ t 的基因组坐标 [u, v] 条件下，最大化MDS的Z分数。但与无监督模式不同的是，此时Z分数是针对所有"组对"（而非"样本对"）计算的。后续步骤与无监督情形下的计算完全相同。

需要指出的是，上述计算面临的一项主要计算挑战，是随着条件数（监督模式下）或样本数（无监督模式，需自动样本分类时）的增加，MDS数量呈二次方增长（补充图1a）。然而在两种模式下，只要两个CpG位点间的距离过大，基因组都会被预先分段。所有此类分段均可被独立分析，从而实现该方法的完全并行化（图1b）。

## 基准测试

**模拟数据上的性能测试。** 为评估metilene³，我们针对人类10号染色体模拟了50个全基因组亚硫酸氢盐测序（WGBS）样本（CpG密集区甲基化水平低，其余区域甲基化水平高），并分别置于两种不同的甲基化背景（全局噪声；详见方法）之下。样本被划分为五组，并根据分组将一组模拟DMR"注入"每种背景中（详见方法）。具体而言，我们引入了五种不同类型的DMR，每种类型200个区域，共计1000个DMR，并与一个未引入任何DMR的对照组（Ctr）进行对比（图2a，补充数据1）。类型I至IV代表各组分别呈现低甲基化、高甲基化或中间甲基化水平的场景。相比之下，类型V反映了更复杂的模式，其中某一组内的DMR在甲基化水平和长度上存在差异，从而在核心DMR周围产生偏离（即"侧翼区域"，图2a）。考虑到在实际生物组织样本中，DMR内部的CpG位点很少完全甲基化或完全未甲基化，我们通过将每种DMR类型分配至四个混合水平之一（混合因子c的取值范围为1至0.6）引入了异质性。总体而言，较低的混合因子会缩小组间平均甲基化差异，从而增加准确检测DMR的难度。为进一步模拟真实生物样本中常见的复杂性，我们还引入了两类批次效应（1000个批次特异性未甲基化位点）、样本特异性的肿瘤纯度因子（DMR内样本与对照甲基化值的混合，因样本而异），以及两个混杂因子C1和C2（分别对应500个与C1或C2相关的混杂特异性DMR，详见方法，补充图2a–d）。

我们首先仅在监督模式下，将metilene³与其他多组DMR调用工具进行比较。总体而言，metilene³在背景1和背景2中分别识别出1116个和957个DMR。相比之下，wgbstools分别检测到920个和255个DMR，而SMART2分别识别出577个和22个DMR（所有方法的错误发现率均约为1%或更低）。值得注意的是，即便在混合因子较低（即模拟难度较大）的情形下，metilene³在所有DMR类型中都保持了较高的灵敏度（图2b）。我们进一步利用Jaccard指数评估模拟DMR与预测DMR之间的一致性。metilene³在背景1中的相似度达到0.96，在背景2中为0.83，而其他DMR调用工具在两种背景下均未达到0.5（图2c）。我们还针对已报告与遗漏DMR的长度进行了分析，结果显示metilene³和SMART2在真阳性方面模拟长度与预测长度之间具有很高的相关性（均>0.8）；相比之下，wgbstools则倾向于低估DMR的长度（图2d左）。

我们还通过将metilene³局部聚类得到的样本标签（低甲基化、高甲基化或中间态）与模拟标签进行比较，评估了其准确性。总体准确率较高（背景1/2分别为0.901/0.882），但不同DMR类型间的准确率存在差异：组间差异更明显的类型I和II被更准确地恢复，而噪声更大的类型III、IV和V准确率有所下降（补充图2e），这很可能是由于组间模拟差异更为细微所致。针对类型V的DMR，我们研究了组特异性DMR长度对单个模拟DMR内预测DMR数量的影响。正如预期，随着组特异性区域长度（"侧翼区域长度"）的增加，metilene³倾向于检测到更多的DMR（补充图2f、g）。

在无标签DMR识别方面，我们将metilene³的无监督模式与methylscore进行了基准比较。methylscore遗漏了大量DMR（在不同混合因子和背景下分别为8.0–91.6%和18.4–98.8%），且预测DMR的长度与模拟长度差异很大（背景1/2的平均偏差分别为341 bp/363 bp，图2d右，补充图2h）。相比之下，metilene³的无监督模式表现出与监督模式相同的准确度，这是因为DMTree正确推断出了分组，这很可能得益于其在聚类时聚焦于强而稳定的DMR。

最后，我们利用模拟数据测量了metilene³、wgbstools和SMART2的时间与内存消耗。在监督式DMR检测中，metilene³在单核上大约需要29分钟，占用不到0.7 GB内存（10核：12分钟，1.26 GB内存）。相比之下，wgbstools更快，单核约需1.6分钟，占用0.48 GB内存，且性能会随核心数增加进一步提升。SMART2则明显更慢，在10核上完成相同任务需要一个多小时（补充图2i）。由于methylscore的工作流需要两个不同的软件版本以及对中间文件的手动处理，其时间与内存消耗无法测量。

**无监督模式——DMTree聚类的评估。** 无监督模式的核心步骤，是利用我们提出的新型聚类方法DMTree，基于在所有样本间（无标签）检测到的DMR，将样本分配到不同的组。DMTree通过在众多DMR间汇集分组证据，使得即便是精细但稳健的信号也能在各次分裂中不断累积，从而提高信噪比。为将本算法与DNA甲基化数据分析中常用的聚类方法进行基准比较，我们使用先前模拟的数据，将所得的分组结果与层次聚类和K-均值聚类的结果进行了比较（详见方法）。在两种背景下，DMTree都能准确再现模拟分组，而层次聚类和基于前1%高变异CpG的K-均值聚类则均受到批次效应的影响，在噪声更大的背景2中尤为明显（图3a–c，补充图3a–c）。

为进一步检验各聚类方法的稳健性，我们进行了重复子抽样实验（抽取80%样本，重复50次，详见方法）。为同时评估小群组检测的灵敏度，每轮实验中G1组的大小被固定为1至10之间的某个数值（补充图3d）。DMTree的聚类结果更为稳健，其平均Rand指数超过0.9，而层次聚类和K-均值聚类分别低于0.7和0.4，并且无论各算法自身参数n、N、K如何设置，DMTree的表现均优于其余两种方法（n为metilene³中的最小组大小，N为层次聚类的分组数，K为K-均值聚类的簇数，图3d）。有趣的是，当以用于构建DMTree的那些DMR作为输入时，层次聚类和K-均值聚类的表现有所改善，但仍未达到DMTree的稳健水平（补充图3e）。重要的是，DMTree对小群组的检测也更为敏感：只要最小聚类样本数n设置得不大于G1组的大小，DMTree便已能在G1样本占比仅为5%时区分出G1与G2（与G1最相似的组，图3e）。而对于层次聚类，需要G1占比超过15%才能使G1和G2的可恢复度（平均Rand指数）达到0.5；对于K-均值聚类，即便设置较大的簇数K，也需要约25%的占比（图3e）。

考虑到真实数据集常常在样本纯度上存在差异（例如由于血液污染或组织混合），我们利用TCGA队列来捕捉常见的纯度水平（最低纯度水平的中位数约为0.6）。据此，我们又模拟了四个纯度递减的数据集，平均纯度范围为0.8至0.5（补充图3f）。结果显示，在纯度不低于0.6的数据集上，DMTree能够准确预测样本分组，而一旦纯度低于该阈值，分组恢复能力则会下降（补充图3g）。因此，常见的纯度水平大体上是可以被良好接受的，但对于极端低纯度的数据集，结果应谨慎解读。

### DMTree再现了细胞谱系与样本来源

为评估metilene³揭示有意义生物学洞见的潜力，我们利用无监督模式分析了一个大型正常人类细胞类型WGBS数据集[13]。该队列包含来自77种原代细胞类型的205个样本，涵盖反映发育过程和组织来源的39个更广泛类别，例如T细胞分化不同阶段的样本。所得的DMTree将样本划分为56个聚类簇，与其已知的细胞类型、功能及组织来源表现出显著的一致性（补充图4a）。例如，metilene³正确地将胰岛细胞类型聚在一起，同时清晰区分了α、β和δ细胞。与此同时，metilene³还报告了每种细胞亚型特异的DMR，其中包括一个位于CREBBP基因座、仅在α细胞中低甲基化的DMR，以及一个位于ASAP1基因座、仅在β细胞中低甲基化的DMR（补充图5a）。

该图谱虽然涵盖了广泛的组织类型，但大多数研究聚焦于特定的细胞类型或器官。为评估metilene³区分高度相似细胞类型的能力，我们还将其应用于36份标注为五种不同细胞类型、涵盖多个分化阶段的血液样本（图4a）。该分析揭示出6个聚类簇，并识别出363,861个DMR（补充数据2）。所得的DMTree再现了关键的谱系关系，共有五次主要分裂：首先将髓系细胞与淋巴细胞分开，随后将粒细胞与单核细胞分开，接着将B细胞与其他淋巴细胞分开，再将NK细胞与T细胞分开，最后将幼稚T细胞与成熟T细胞分开（后者在图谱原始分析中未被区分，图4b）。值得注意的是，当使用metilene³推导出的DMR的甲基化水平时，样本在PCA中按细胞类型清晰分组；而使用全部CpG、前1%高变异CpG，或与构建DMTree数量相当的CpG时，则均未观察到明显的分离（补充图5b）。这凸显了所识别DMR的细胞类型特异性，并引出了一个问题：这些DMR是否可能在造血谱系分化中具有功能意义。

为评估metilene³所识别DMR的功能相关性，我们研究了这些区域内潜在的转录因子结合位点是否与血细胞特异性功能相关。为此，我们对与DMTree分裂相关联的DMR进行了motif富集分析。值得注意的是，我们在分裂I（髓系细胞中低甲基化）中发现ATTGCGCAAC基序显著富集，该基序被C/EBPβ识别；在分裂V（幼稚T细胞中低甲基化）中发现CCTTTGATST基序显著富集，其与淋巴增强因子1（LEF1）的DNA结合基序相似（图4c）。为评估这些转录因子是否也与细胞类型特异性表达模式相关，我们分析了一个已发表的人类外周血单个核细胞单细胞RNA测序数据集（图4d）。该数据集包含淋巴细胞（T、B及NK细胞）和单核细胞，但不包含粒细胞。与既往研究一致[19]，C/EBPβ在单核细胞中的表达显著高于淋巴细胞，支持驱动DMTree首次主要分裂的DMR反映了具有生物学意义的谱系分离这一观点（图4e）。值得注意的是，已有研究表明C/EBPβ能够诱导v-ABL永生化的原代B细胞（通常属于淋巴系）转分化为具有髓系特征的粒细胞-巨噬细胞祖细胞样细胞（GMP样细胞）[20,21]。在谱系进一步下游，LEF1在T细胞发育中发挥关键作用，负责维持基因组结构并确保其身份与功能的稳定性[22]。事实上，正如既往观察所示，LEF1在幼稚CD4⁺ T细胞中的表达最高[23]（图4e）。这些观察结果表明，metilene³识别出的DMR捕获了参与塑造细胞身份的相关调控区域，从而为理解潜在的基因调控网络提供了洞见。

### 在癌症数据集中的应用

癌症样本通常比健康组织具有更高的异质性，使得DMR检测更具挑战性。为评估metilene³在此类复杂背景下用于生物学假设验证的能力，我们将其应用于两个公开的癌症相关WGBS数据集[24,25]，以及一个源自游离DNA的数据集[26]。

**胶质母细胞瘤。** 来自Wu等人[24]的数据集包含64个样本，来自IDH突变型和IDH野生型胶质母细胞瘤肿瘤及非癌脑组织的混合队列（图5a）。无监督DMTree识别出267,811个DMR，并将样本划分为四组（图5b，补充数据3）。簇A包含IDH野生型肿瘤，簇B和簇C包含IDH突变型肿瘤，簇D主要由非癌样本组成，表明该方法能够根据病理生物学特征和IDH突变状态正确地对样本进行分层（图5b）。为评估IDH突变相关的簇B和簇C是否代表具有生物学差异的亚组，我们计算了差异表达基因（DEG），并对DMR和DEG均进行了基因集富集分析（GSEA，补充数据4和5）。GSEA提示存在潜在的分子亚型，它们在癌症相关通路上呈现出差异富集，例如簇B与参与上皮-间质转化（EMT）的基因相关。值得注意的是，基于DMR计算得到的优势比与基于DEG计算的GSEA网络富集分数显著相关（Spearman R = 0.77，p < 0.0001），提示甲基化变化与基因表达变化之间具有一致性（图5c）。

针对这些DMR的motif分析提示，一个Krüppel型C2H2锌指转录因子（ZNF766，p < 10⁻²⁷）可能参与了所检测出的簇B和簇C之间的甲基化差异（补充图6c）。ZNF766在健康脑组织及其他组织中广泛表达[27]，近期研究表明其可能干扰PRDM9介导的重组过程[28]，因而可能与染色体不稳定性相关[29,30]。有趣的是，另一个富集的转录因子CTCFL（BORIS）通常仅在精子发生过程中表达，据报道可被PRDM9激活[28,31]。我们还观察到CDX1（一种与增殖和凋亡相关的转录因子）的显著富集（p < 10⁻²⁷）。在其他转录因子中，我们的数据提示ZNF766、CTCFL和CDX1在IDH突变型胶质母细胞瘤B簇中存在结合，提示其相对于簇C中其余IDH突变型胶质母细胞瘤可能存在潜在的调控差异。

此外，我们还识别出一个明显的IDH突变型肿瘤离群样本AK015，其在DMTree中与正常样本聚在一起。该样本还表现出非典型的转录组特征：增殖标志物MKI67的表达水平与正常样本一样低，提示其增殖活性低于其他肿瘤（图5d）。然而，胶质瘤相关抗原TNC仍保持高表达，这与单纯的肿瘤纯度解释相悖，支持其肿瘤身份（图5d）。此外，特别是在IDH突变特异性DMR方面，AK015偏离了其他正常样本，表现出部分向癌症方向的偏移，而其他癌症特异性区域则表现得更接近正常状态（图5b）。基于全部CpG、高变异CpG或DMR CpG的PCA，以及转录组数据，均证实了这种中间态特征（补充图6a、b）。虽然阐明该样本非典型特征的确切原因超出了本文范围，但metilene³将该样本标记为离群值，为进一步研究其潜在原因提供了线索。

最后，我们测试了metilene³在应用于因样本纯度可能较低而更具挑战性的样本时的敏感性。为此，我们使用了一个源自脑脊液（CSF）循环游离DNA（cfDNA）的脑癌数据集[26]，即源自脑肿瘤或邻近组织、释放到脑脊液中的小DNA片段。我们在9个cfDNA样本上测试了metilene³，其中4个来自非肿瘤供体（脑积水患者），5个来自处于不同治疗阶段的髓母细胞瘤患者；DMTree成功将其聚类为肿瘤组和非肿瘤组（补充图7，补充数据6）。

**胰腺导管腺癌。** 我们重新分析了Lo等人[25]的数据，该数据包含35个来自腺泡细胞、导管细胞、胰腺上皮内瘤变（PanIN）以及胰腺导管腺癌（PDAC）的样本（图6a）。metilene³识别出275,352个DMR，并构建出一棵能够准确区分癌症与非癌样本、瘤变与正常样本、以及腺泡样本与导管样本的DMTree（图6b，补充数据7）。基于DMR甲基化水平的PCA证实了瘤变样本与两种正常类型样本之间的分离（图6c右）。相比之下，基于全部CpG或前1%高变异CpG的PCA并未显示出腺泡样本与瘤变样本之间的清晰分离，而是主要捕捉到了潜在的文库制备批次效应（样本名中含A表示Swift Accel-NGS Methyl-Seq，含B表示NEBNext Ultra DNA，图6c）。

鉴于关于个体PDAC起源细胞（腺泡还是导管）的争论仍在持续[32]，我们利用腺泡-导管DMR，计算了反映瘤变和PDAC样本与腺泡或导管样本之间甲基化相似性的一致性得分。与腺泡-导管化生（ADM）模型[32]一致，PanIN样本始终与腺泡细胞更为相似，而PDAC样本则更接近导管细胞（图6d），支持了PanIN代表肿瘤进展过程中从腺泡向导管身份转变的中间态这一观点。值得注意的是，当限定于腺泡-导管差异更大的DMR时，PanIN和PDAC样本与导管细胞的相似性进一步增加。这表明，在正常腺泡细胞和导管细胞差异最大的区域，PanIN获得了类导管的甲基化特征，凸显了其在PDAC发展中的过渡状态。对腺泡-导管DMR的进一步分层显示，在腺泡细胞中低甲基化的区域甲基化水平逐渐升高，而在导管细胞中低甲基化的区域则更为骤然地去甲基化（图6e）。在所有分析中，PanIN的甲基化水平始终介于腺泡和PDAC之间，支持了在腺泡-导管化生（ADM）[33]框架内存在一条腺泡-PanIN-PDAC的进展轨迹。

接下来，我们研究了DMR的组成和基因组分布。得益于其局部聚类策略，metilene³检测到了位置相邻但在各组间表现不同的DMR。例如，位于TRAK1转录起始位点下游的一个DMR在PDAC和导管样本中呈低甲基化，但在IDH突变亚簇B和C中并非如此（原文如此，此处指PDAC样本中该DMR在PanIN和腺泡样本中未呈低甲基化），而相邻的一个DMR则仅在腺泡细胞中呈高甲基化（图6f）。在TXNRD1的一个内含子中也观察到类似模式：三个DMR均显示PDAC中低甲基化，而两种正常组织均呈高甲基化（图6g）。在瘤变样本中，其中一个DMR也呈低甲基化，而其余DMR则保持接近正常的甲基化水平。在这两个位点，这些表观遗传变化均伴随着表达水平的逐步上调，与ADM模型一致（补充图8a）。TRAK1编码一种线粒体运输驱动蛋白，已被证明与其他癌症的侵袭行为有关[34]，可能与PDAC的侵袭性相关[35]。同样，TXNRD1编码硫氧还蛋白还原酶1，是一种与癌症相关的酶，也是获FDA批准的治疗靶点，参与胰腺β细胞的抗氧化防御[36]，其过表达通过AKT/mTOR激活与肝细胞癌的不良预后相关[37]。尽管这两个基因均未被直接与PDAC相关联，但它们表观遗传上的"双重开关"模式提示其可能在疾病机制中发挥作用，值得进一步研究。

为表征DMTree的首次分裂，我们检索了可能参与或受PDAC样本中差异甲基化影响的转录因子结合位点（图6b）。对PDAC-DMR的motif富集分析显示，在PDAC强烈高甲基化区域（平均甲基化增加>0.5，n=125）中未发现显著富集。相比之下，在PDAC强烈低甲基化DMR（平均甲基化丢失>0.5，n=524，补充图9a）中，若干motif表现出显著富集。以其余所有DMR（n=8661）为背景，我们观察到核因子κB2（NFKB2；OR=3.35，p=8.3×10⁻²³）和活化T细胞核因子1（NFATc1；OR=2.26，p=2.7×10⁻¹¹）基序在PDAC低甲基化DMR中显著富集（补充图9b）。已知NF-κB信号通路通过抑制凋亡、促进增殖、血管生成和转移来促进癌变[38]，NFAT家族因子（NFATc1–4）则调控细胞周期、存活、迁移及血管生成，并被认为与PDAC相关[39,40]。值得注意的是，34%的PDAC低甲基化DMR（179/524）含有NFKB2或NFATc1基序，其中8%（44/524）同时含有两种基序，相对背景呈现显著富集（OR：9.83；Fisher精确检验p值=3.8×10⁻²⁴，补充图9b，补充数据8）。在全基因组范围内，仅1.6%的NFATC1基序位于NFKB2基序30 bp范围内，而在PDAC低甲基化DMR中这一比例为38.7%。值得注意的是，这两个结合位点都属于一组表现出显著共现的转录因子结合位点，这种共现无法用DMR长度来解释，因为PDAC低甲基化DMR的平均长度反而短于其他DMR（补充图9c、d）。

我们接下来探究这些转录因子结合位点是否与差异甲基化更普遍地相关。在所有DMR（不仅是PDAC低甲基化DMR）中，我们观察到以NFKB2和NFATc1基序为中心的对称低甲基化（补充图9e左）。重要的是，Arlt等人的研究已将NF-κB和NFAT（连同Nrf2）确定为促进PDAC发展的主要核因子[41]，而NFKB或NFAT基序的低甲基化在PDAC中最为显著，在PanIN中程度较轻（在腺泡和导管样本中则明显更弱）。此外，携带这些基序的DMR邻近富集于WNT/β-catenin及通过NF-κB介导的TNFα信号通路的基因（补充图9e右）。在同时含有两种基序的DMR中，PanIN的甲基化水平几乎接近正常组织（图7a），提示双基序低甲基化比单基序位点更具PDAC特异性。与此一致，公开的RNA-seq数据[42,43]显示，NFKB2和NFATC1（及其家族成员）的表达沿ADM连续谱逐渐上调，并在PDAC中达到峰值（图7b；补充图10a）。将metilene³检测到的DMR与TCGA-PAAD样本和配对正常组织之间的DEG（来自GEPIA2[44]）进行整合，我们发现midkine基因（MDK）是与双基序DMR相关联的、最显著上调的DEG（图7c及补充图10b、c，补充数据9）。既往研究表明，MDK参与促进EMT[45]、在胰腺癌细胞胞外囊泡中富集[46]，以及促进胰腺癌的增殖和迁移[47]。

综合来看，这些数据支持了一个具有假说生成价值的模型：PDAC特异性的双基序低甲基化标记了MDK及其他癌症进展相关基因附近的调控元件，并可能与NF-κB/NFAT活性存在耦合（可能通过Akt/mTORC1/NF-κB通路[45]；补充图10d）。除MDK外，metilene³还揭示了具有复杂甲基化模式的其他位点（如TRAK1、TXNRD1），可作为生物标志物开发及机制研究的候选靶点（补充数据9）。

## 讨论

在本研究中，我们提出了metilene³，用于分析全基因组甲基化数据，实现了多条件无监督DMR识别，并构建了基于DMR的树（DMTree）用于样本聚类。基于改进的环形二元分割算法，metilene³在模拟数据中取得了出色的性能（背景1/2中灵敏度分别为96.8%/83.2%，精确度分别为98.9%/99.0%）。此外，DMTree准确聚类了所有模拟样本，以及来自正常人体组织和肿瘤的样本。通过将无标签DMR检测与以DMR为锚点的聚类相结合，甲基化分析变得更少依赖先验知识、更少偏倚。重要的是，DMTree提供了直接的可解释性，因为单个DMR被明确关联到树的分裂节点，使聚类结果更加直观、可解释和可验证。最后，metilene³运行快速、内存高效：使用10核处理35个样本的人类癌症WGBS数据集约需10分钟，占用内存不到1 GB，可支持探索性数据分析中的参数调优。

我们证明，将metilene³应用于生物医学数据集有潜力提供有关基因组调控相互作用的重要洞见。DMTree再现了血细胞谱系关系，并在一个人类细胞图谱[13]中重新识别出细胞类型特异性的调控区域、因子和靶点。同样，在癌症和cfDNA数据集等更具挑战性的场景中，metilene³识别出了合理的DMR和DMTree，支持了探索性分析和生物学假设的形成。需要指出的是，与任何聚类方法一样，连续的轨迹会被表示为离散的分支；然而，其潜在的相似性在树的拓扑结构中得以保留，未来的扩展可以引入基于扩散的方法来显式建模连续数据。在胰腺组织中，metilene³清晰区分了胰腺腺癌、瘤变及对照组织样本。其中，metilene³报告了高度具有区分力的癌症相关DMR，这些DMR富集于提示NFKB和NFAT共结合的motif。结合RNA表达数据的结果，让我们不禁推测存在一条潜在的midkine–NF-κB–NFAT信号轴，参与胰腺肿瘤的发生。除通路层面的信号外，metilene³还揭示了具有复杂甲基化模式的、能够区分健康、瘤变和恶性组织的单个位点。这些位点是极具前景的机制研究后续靶点和潜在诊断生物标志物候选。

未来工作将聚焦于将metilene³的应用范围扩展至其他技术平台。由于获取高覆盖度的牛津纳米孔技术（Oxford Nanopore）或PacBio直接DNA测序数据存在困难，我们尚未对该场景进行专门测试。然而，已有若干研究表明，这些技术之间的甲基化率具有可比性，因此我们推测不会存在显著偏倚。随着单细胞亚硫酸氢盐测序（scBS）技术的发展，得益于其标准化的数据输入格式以及处理大样本量的能力，metilene³已具备处理scBS数据的准备；然而，迄今为止，该技术中较高的脱失率（dropout rate）对先验数据插补构成了重大挑战。

综上所述，metilene³为在先验知识有限的情况下分析甲基化数据，提供了一种敏感、准确且高效的解决方案，有望促进生物医学研究中新假说与新发现的产生，并有待在实验或转化研究场景中作进一步验证。

## Methods

### Segmentation

We modified the circular binary segmentation described in Jühling et al.⁵ Here, we first introduce our basic approach to segmentation and subsequently explain how this is done in the unsupervised mode.

#### Basic pair-wise segmentation strategy

Let $E$ denote the set of samples. For some group of samples, i.e., $G \subset E$, we calculate the mean methylation level at the cytosine $i$ with

$$
M_i(G) = \frac{1}{|G|}\sum_{e \in G} \beta_{i,e} \tag{1}
$$

- **这个公式在做什么：** 它定义的是"某一组样本，在某一个CpG位点上的平均甲基化水平"——说白了就是对一组样本，在同一个位置求个平均值。这是后续所有比较（组间差异、Z分数、DMR判定）的最基础的统计量。

- $M_i(G) = \text{组 } G \text{ 中所有样本在位点 } i \text{ 处甲基化值的平均数}$

- 这一步本质上就是"对一组人、在同一个基因位置，求个平均甲基化率"，是后面公式(2)计算"两组之间差多少"（$\Delta_i(G,H) = M_i(G) - M_i(H)$）的直接输入。

- **$E$**：全部样本组成的集合（比如你测序的所有个体/细胞样本）。

- **$G$**：$E$ 的一个子集（$G \subset E$），代表"某一组样本"——比如"肿瘤组"或无监督模式下"局部分类出来的某个组"。

- **$i$**：基因组上某一个具体的胞嘧啶（CpG）位点的索引，也就是"你正在看第几个CpG"。

- **$e$**：集合 $G$ 中的某一个具体样本（$e \in G$），是求和时的遍历变量，代表"组内的第几个个体"。

- **$\beta_{i,e}$**：样本 $e$ 在位点 $i$ 处的甲基化值（β值）。这是最原始的观测数据——通常是一个0到1之间的比例，表示这个位点在这个样本里"被甲基化的读段占比"（0=完全未甲基化，1=完全甲基化）。

- **$|G|$**：集合 $G$ 的样本数量（即组内有多少个样本）。

- **$\sum_{e \in G} \beta_{i,e}$**：把组 $G$ 里**每一个样本**在位点 $i$ 处的甲基化值 $\beta_{i,e}$ **全部加起来**。

- **$\frac{1}{|G|}\sum_{e \in G} \beta_{i,e}$**：把上面的总和**除以样本数**，也就是求算术平均值。

- **$M_i(G)$**：最终结果——组 $G$ 在位点 $i$ 处的**平均甲基化水平**。

---

where $\beta_{i,e}$ is methylation value of CpG $i$ in the sample $s$. Consequently, the mean difference signal (MDS), $\Delta$, for cytosine $i$ between two groups $G, H \subset E$ with $G \cap H = \emptyset$ is calculated by

$$
\Delta_i(G,H) = M_i(G) - M_i(H) \tag{2}
$$

- **这个公式在做什么：** 它定义的是"平均差异信号"（Mean Difference Signal, MDS）——在同一个CpG位点上，把两组样本各自的平均甲基化水平相减，得到这两组在这个位置上"差了多少"。这是整篇论文里所有分段（segmentation）、Z分数计算的核心原料：先在每个位点算出一个"差异值"，后面的算法本质上都是在这条差异信号曲线上找"哪一段差异最显著、最持续"。

- **$G, H \subset E$**：两个都属于全部样本集合 $E$ 的子集，代表要拿来做比较的"两组样本"（比如无监督模式里种子样本对应扩展出的两组，或监督模式里用户指定的两个条件组）。

- **$G \cap H = \varnothing$**：$G$ 和 $H$ 没有交集，也就是同一个样本不能同时属于这两组——这是"两组"这个概念成立的前提条件。

- **$i$**：基因组上某一个具体的CpG位点索引，和公式(1)里的 $i$ 是同一个意思。

- **$\beta_{i,e}$**：样本 $e$ 在位点 $i$ 处的甲基化值（延续公式1的定义；原文这里写的"sample $s$"是笔误，应该是任意样本 $e$）。

- **$M_i(G)$**：由公式(1)算出的、组 $G$ 在位点 $i$ 处的**平均**甲基化水平。

- **$M_i(H)$**：同理，组 $H$ 在位点 $i$ 处的**平均**甲基化水平。

- **$M_i(G) - M_i(H)$**：两组平均值直接相减——**做减法**，得到"$G$组比$H$组在这个位点甲基化程度高多少（或低多少）"。

- **$\Delta_i(G,H)$**：相减后得到的结果，即位点 $i$ 处的**平均差异信号（MDS）**。

---

and the sum of such differences over sum interval $[u,v]$ is

$$
S_{u,v}(G,H) = \sum_{i=u}^{v} \Delta_i(G,H) = \sum_{i=0}^{v} \Delta_i(G,H) - \sum_{j=0}^{u-1} \Delta_j(G,H) \tag{3}
$$

- **这个公式在做什么：** 而公式(3)是把公式(2)这些逐点的差异值，在一段基因组区间 $[u,v]$ 内**累加起来**，得到"这一整段区域上，两组样本的总差异信号" $S_{u,v}$。这个累加和，正是后面公式(4)(5)里Z分数要去最大化的核心统计量——算法要找的，就是使这个累加差异"相对区间长度而言最大"的那一段 $[u,v]$。

- **$[u,v]$**：一段基因组坐标区间，$u$ 是起点，$v$ 是终点（满足 $u \le v$）。这是算法要去搜索、试探的"候选区段"。

- **$i$**：求和用的遍历变量，代表区间内的每一个CpG位点，从 $u$ 一直取到 $v$。

- **$\Delta_i(G,H)$**：公式(2)定义好的、位点 $i$ 处两组的平均差异信号。

- **$\sum_{i=u}^{v} \Delta_i(G,H)$**：把区间 $[u,v]$ 内**每一个**位点的差异值 $\Delta_i$ **依次相加**，得到该区间的总差异。这是公式最直接、最容易理解的写法。

- **$S_{u,v}(G,H)$**：上面求和的结果，即区间 $[u,v]$ 上的**累积差异信号（cumulative MDS）**。

- **$\sum_{i=0}^{v} \Delta_i(G,H)$**：从基因组（或该染色体片段）最开始的位点 $0$，一直累加到 $v$ 的**前缀和**。

- **$\sum_{j=0}^{u-1} \Delta_j(G,H)$**：同样从 $0$ 开始，累加到 $u-1$（也就是 $u$ 的前一个位点）的**前缀和**；这里用 $j$ 而不是 $i$，只是为了避免和左边的求和变量混淆，本质上是同一个意思。

- **两者相减**：一个"到 $v$ 为止的总和"减去"到 $u-1$ 为止的总和"，剩下的正好就是 $[u,v]$ 这一段区间的和——这是一个经典的**前缀和技巧（prefix sum）**。

### 为什么要写成这种"前缀和相减"的形式

这不是数学上必须的，而是**工程/计算效率上的考量**：如果每次要算某个区间 $[u,v]$ 的和都从头逐个位点相加，当算法要反复尝试成千上万种不同的 $u, v$ 组合时会非常慢。而如果预先算好从头到每个位置的前缀和，之后任意区间的和只需要"两个前缀和相减"就能瞬间得到——这正好呼应了论文里反复强调的"metilene³计算高效、可全基因组并行"的设计目标。

```text
M_i(G)  = [0.82  0.79  0.847 0.113 0.103 0.877 0.897 0.15  0.203 0.193]
M_i(H)  = [0.21  0.17  0.16  0.76  0.79  0.21  0.2   0.785 0.8   0.75 ]
Delta_i = [0.61  0.62  0.687 -0.647 -0.687 0.667 0.697 -0.635 -0.597 -0.557]

前缀和 prefix = [0.61  1.23  1.917 1.27  0.583 1.25  1.947 1.312 0.715 0.158]

[u, v]        直接求和      前缀和相减    是否一致
[ 0, 2]      1.9167        1.9167       True
[ 3, 4]     -1.3333       -1.3333       True
[ 0, 9]      0.1583        0.1583       True
[ 5, 8]      0.1317        0.1317       True
[ 3, 9]     -1.7583       -1.7583       True

换 [5, 8] 举例
```

---

Subsequently, the MDS in some interval $[s,t]$ is segmented using a Z-score. Specifically, the algorithm determines the coordinates $[u,v]$, $s \le u < v \le t$, maximizing the score for two given groups $G, H$ with

$$
\underset{s \le u < v \le t}{\operatorname{argmax}} \; Z_{s,t}(u,v \mid G,H) \tag{4}
$$

- **这个公式在做什么：** 前面公式(3)算出的 $S_{u,v}(G,H)$ 只是"某个候选区间内的差异总和"——但一个大区间天然容易累积出更大的和，这并不代表它就是真正有意义的差异区域。公式(4)要做的是：在一个更大的搜索范围 $[s,t]$ 内，**尝试所有可能的子区间 $[u,v]$**，找出使"Z分数"（公式5会定义）**最大**的那一对坐标 $u,v$——这一步就是CBS分割算法里"找断点"的核心：不是找和最大的区间，而是找**统计显著性最强**的区间。

- **$[s,t]$**：当前正在处理的、较大的一段基因组搜索区间（"母区间"）。算法是递归进行的：一开始 $[s,t]$ 可能是整条染色体，每次分割出DMR后，剩下的两侧子区间又会分别作为新的 $[s,t]$ 继续搜索。

- **$[u,v]$**：$[s,t]$ 内部的一个**候选子区间**，也就是算法想要确定的"最优分割结果"。

- **$s \le u < v \le t$**：对 $u,v$ 取值范围的约束——$u,v$ 必须落在 $[s,t]$ 内部（不能超出母区间边界），并且 $u$ 严格小于 $v$（保证区间至少包含从 $u$ 到 $v$ 之间的位点，是个"有长度"的合法区间，不能是空区间或反向区间）。

- **$G, H$**：要比较的两组样本（延续前面公式(2)(3)的定义），这里是"两个已经给定的组"，公式(4)是在这两组固定的前提下去搜索最优区间。

- **$Z_{s,t}(u,v \mid G,H)$**：Z分数函数，输入是候选区间 $[u,v]$（以及它所处的母区间 $[s,t]$）和两组样本 $G,H$，输出一个衡量"这个区间差异有多显著"的数值。竖线 "$\mid$" 是数学里常见的写法，表示"在给定 $G,H$ 的条件下"，也就是这个函数是关于 $u,v$ 的，$G,H$ 是已知的固定参数。它的具体计算公式会在公式(5)里给出。

- **$\operatorname{argmax}$**：这是"取使函数值最大的自变量"的意思，和"$\max$"（取最大值本身）不同——$\max$ 返回的是最大的**数值**，而 $\operatorname{argmax}$ 返回的是达到这个最大值时对应的**坐标** $(u,v)$。

- **$\underset{s \le u < v \le t}{\operatorname{argmax}}$**：合起来读，就是"在满足 $s \le u < v \le t$ 这个范围约束的所有可能的 $(u,v)$ 取值里，找出让 $Z_{s,t}(u,v\mid G,H)$ 取得最大值的那一组 $(u,v)$"。


- 公式(4)描述的是一个**穷举搜索/优化过程**：算法在母区间 $[s,t]$ 内，遍历所有合法的子区间 $[u,v]$，计算每一个候选区间的Z分数，最终**保留Z分数最高的那一个区间**作为这一轮分割的结果——这正是CBS递归分割里"每一步该在哪里切一刀"的判定依据，切出来的这个区间接下来会被送去做2D-KS检验，通过检验才会被最终报告为一个DMR。
---

where

$$
Z_{s,t}(u,v \mid G,H) = \frac{\left[\,|S_{u,v}(G,H)| - \dfrac{v-u}{t-s}\cdot|S_{s,t}(G,H)|\,\right]^2}{(v-u)\left[1-\dfrac{v-u}{t-s}\right]} \tag{5}
$$

- **这个公式在做什么：** 这是公式(4)里要最大化的那个"Z分数"的具体计算方式。直觉上讲，它衡量的是——**候选子区间 $[u,v]$ 实际累积的差异，比"如果差异是均匀分布在整个 $[s,t]$ 上、按长度比例该分到的差异"多出了多少，并且用一个考虑区间长度的分母把这个"多出来的量"标准化成一个可比较的分数**。这正是经典CBS（环形二元分割）算法里用来"找断点"的统计量的推广形式。

- **$[s,t]$**：当前正在搜索的母区间（延续公式4的定义）。
- **$[u,v]$**：$[s,t]$ 内部正在测试的候选子区间。
- **$G,H$**：要比较的两组样本。
- **$S_{u,v}(G,H)$**：公式(3)定义的、子区间 $[u,v]$ 上的**累积差异信号**。
- **$S_{s,t}(G,H)$**：同样用公式(3)计算，但换成整个母区间 $[s,t]$ 上的**累积差异信号**——代表"这一大段区域总共积累了多少差异"。
- **$|S_{u,v}(G,H)|$、$|S_{s,t}(G,H)|$**：取绝对值——因为差异信号可正可负（$G$比$H$高还是低），这里只关心"差异有多大"，不关心方向。
- **$v-u$**：候选子区间 $[u,v]$ 的长度（含多少个位点）。
- **$t-s$**：母区间 $[s,t]$ 的长度。
- **$\dfrac{v-u}{t-s}$**：子区间长度占母区间长度的**比例**，记作一个介于0和1之间的分数，可以理解成"如果差异是均匀分布的，这段子区间理论上应该占多大份额"。


- **$\dfrac{v-u}{t-s}\cdot|S_{s,t}(G,H)|$**：这是一个"期望值"——假设整个母区间 $[s,t]$ 的总差异是按长度**均匀摊派**的，那么长度占比为 $\frac{v-u}{t-s}$ 的子区间，"应该"分到这么多差异。
- **$|S_{u,v}(G,H)| - \dfrac{v-u}{t-s}\cdot|S_{s,t}(G,H)|$**：用"子区间实际观测到的差异"减去"按比例算出的期望差异"，得到一个**偏离量**——如果子区间的差异明显超出了它按长度该分到的份额，说明这段区间"异常集中"地承载了差异，很可能就是真正的DMR所在。
- **平方 $(\cdot)^2$**：把偏离量平方，一是消除正负号（不管偏离方向，只看偏离程度），二是让统计量的形式更接近经典的（近似卡方分布的）检验统计量，方便后续判断显著性。

- **$(v-u)$**：子区间的长度。
- **$1-\dfrac{v-u}{t-s}$**：子区间长度占比的"补集"，即母区间中**不属于**子区间的那部分所占的比例。
- **两者相乘**：这一项的形式类似统计学中"比例的方差"（$n \cdot p \cdot (1-p)$ 的结构），作用是**归一化**——当子区间占母区间的比例接近0或接近1（也就是子区间特别短，或者几乎和母区间一样长）时，分子中的偏离量本来就容易因为随机波动而显得很大或很小；用这个分母去除，能让不同长度、不同位置的候选区间之间的Z分数**具有可比性**，不会因为区间长度不同而产生偏差。

公式(5)本质上是在问："子区间 $[u,v]$ 实际携带的差异，相对于它按长度理应占的份额，**偏离得有多显著**？" 分子衡量偏离的绝对大小，分母根据区间长度做归一化。公式(4)遍历所有可能的 $[u,v]$，找出让这个"显著性分数"最大的那一段——分数越大，说明这段区间的差异越不像是随机波动，越可能是真正的差异甲基化区域（DMR）候选。

---

The segmentation algorithm proceeds recursively until no further increase of a two-dimensional Kolmogorov-Smirnov (KS) test statistic is achieved. To be computationally efficient, MDS with

$$
\sum_{i=u}^{v} \max\left[0, \operatorname{sgn}\!\left(|\Delta_i(G,H)| - \delta\right)\right] < r_{\min}
$$

where $\delta$ is a user-defined methylation difference threshold (default: 0.5) and $r_{\min}$ is a user-defined minimum CpG with difference exceeding $\delta$ (default: 5), will be omitted and no 2D-KS test will be performed. Once a segment is found that satisfies the minimal methylation difference $\delta$ and has no sub-intervals with better test statistics, it is reported as DMR.

- **这个公式在做什么：** 这不是一个"分数"公式，而是一个**提前剪枝/跳过判断**——在正式对某个候选区间 $[u,v]$ 做代价较高的二维KS检验之前，先用这个简单规则快速判断"这段区间有没有必要认真检验"。如果区间内差异明显的CpG数量太少，就直接跳过，不做2D-KS检验，从而节省大量计算量。这也是metilene³能做到"全基因组扫描依然很快"的一个关键工程手段。

- **$[u,v]$**：正在评估的候选区间（延续前面的定义）。
- **$i$**：区间内遍历的CpG位点索引，从 $u$ 到 $v$。
- **$G, H$**：正在比较的两组样本。
- **$\Delta_i(G,H)$**：延续公式(2)的定义，位点 $i$ 处 $G,H$ 两组的平均差异信号。
- **$\delta$**：用户设定的**甲基化差异阈值**（默认 $0.5$）——判断"这个位点差异算不算大"的分界线，作用和公式(7)里的 $\epsilon$ 是同一类东西，只是用在不同的场景（这里是基础两组比较策略，$\epsilon$ 用在无监督模式的距离函数里）。
- **$|\Delta_i(G,H)| - \delta$**：实际差异减去阈值，正数表示"超过阈值"，负数表示"没超过"。
- **$\operatorname{sgn}(\cdot)$**：取符号——超过阈值为 $+1$，没超过为 $-1$，正好相等为 $0$。
- **$\max[0, \operatorname{sgn}(\cdot)]$**：把 $-1$ 截断成 $0$，只保留 $+1$（或 $0$），等价于一个"是否超标"的指示器：超标记1分，不超标记0分。
- **$\sum_{i=u}^{v}[\cdot]$**：把区间内所有位点的这个0/1得分加起来，得到**区间 $[u,v]$ 内差异超过 $\delta$ 的CpG位点总数**——注意这里**没有**像公式(7)那样除以区间长度 $(v-u)$，算出来的是一个原始计数，不是比例。
- **$r_{\min}$**：用户设定的**最小CpG数量阈值**（默认 $5$）——要求"差异超标的位点数"至少要达到这个数量。
- **$< r_{\min}$**：判断条件——如果超标位点数量**小于** $r_{\min}$，说明这段区间里"差异足够大的位点"太少，不太可能是一个有意义的DMR，于是这段区间会被**直接跳过**，不再耗费算力去做2D-KS检验。

- 这条规则做的事情是：**先数一数区间内"差异明显超过 $\delta$ 的CpG有几个"，如果这个数量连 $r_{\min}$ 都不到，就认为这段区间大概率不是DMR，直接跳过、省下一次昂贵的统计检验**——是一个纯粹为了计算效率服务的"预筛选"步骤，和公式(4)(5)里"找最优切分点"的Z分数、以及公式(7)里"局部样本分类"的伪距离，是三套服务于不同目的、但结构上很相似（都用了 $\operatorname{sgn}$ 和阈值比较）的计算工具。

---

### Segmentation in the unsupervised mode

In the unsupervised mode, for some interval $[s,t]$, metilene³ assigns the two samples maximizing the Z-score in $[s,t]$, $O$ and $R$, to groups $G^0$ and $H^0$, respectively. After determining the Z-score maximizing interval $[u,v]$, metilene³ proceeds by clustering the remaining samples in $E \setminus \{G^0 \cup H^0\}$. The cluster index $C$ assigning the samples to one or none of the two groups is calculated by

$$
C_{u,v}(e) = \operatorname{sgn}\!\big(d_{u,v}(e, G^0)\big) - \operatorname{sgn}\!\big(d_{u,v}(e, H^0)\big) \tag{6}
$$

- **这个公式在做什么：** 这是无监督模式下"局部样本分类"的核心判定规则。前面通过Z分数搜索，已经找到了一对最具代表性的种子样本 $O$（归入组 $G^0$）和 $R$（归入组 $H^0$），以及它们所在的最优区间 $[u,v]$。现在要处理**剩下的其他样本** $e$：判断每一个 $e$ 到底该归入"像 $G^0$"的一边、"像 $H^0$"的一边，还是两边都不像（归入中间组）。公式(6)给出的 $C_{u,v}(e)$ 就是这个判定结果的数值编码——它本身还不是最终分类，而是一个**中间指标**，后面再根据它的正负号来决定分组。

- **$e$**：待分类的某一个样本（不是种子样本 $O$、$R$ 本身，而是集合 $E \setminus \{G^0 \cup H^0\}$ 里剩下的样本）。

- **$[u,v]$**：前面Z分数最大化搜索出的最优区间，分类判断就是在这个区间范围内进行的。

- **$G^0$**：以种子样本 $O$ 为代表的组（"高甲基化组"或"低甲基化组"其中之一，取决于哪个方向）。

- **$H^0$**：以种子样本 $R$ 为代表的另一组。

- **$d_{u,v}(e, G^0)$**：一个"伪距离函数"（具体定义在公式7），衡量样本 $e$ 在区间 $[u,v]$ 内与组 $G^0$ 有多"接近"或多"不同"。这个函数会在下一条消息里详细拆解，这里先把它当作一个已知的距离度量来用。

- **$d_{u,v}(e, H^0)$**：同理，衡量样本 $e$ 与组 $H^0$ 的距离。

- **$\operatorname{sgn}(\cdot)$**：符号函数（sign function），只保留输入值的**正负号**，不保留具体数值大小：
  - 输入为正 → 输出 $+1$
  - 输入为负 → 输出 $-1$
  - 输入为 $0$ → 输出 $0$

- **$\operatorname{sgn}\big(d_{u,v}(e,G^0)\big)$**：把"$e$与$G^0$的距离"转换成一个只表示方向的信号（$e$相对$G^0$是"偏离"还是"接近"）。

- **$\operatorname{sgn}\big(d_{u,v}(e,H^0)\big)$**：同理，把"$e$与$H^0$的距离"也转换成方向信号。

- **两者相减**：$\operatorname{sgn}(d_{u,v}(e,G^0)) - \operatorname{sgn}(d_{u,v}(e,H^0))$——比较这两个方向信号，得到样本 $e$ 更偏向哪一组的最终指标。

- **$C_{u,v}(e)$**：最终结果，一个取值范围在 $\{-2,-1,0,1,2\}$ 之间的整数（两个 $\operatorname{sgn}$ 值各自是 $-1,0,1$ 之一，相减后的可能组合），用来编码样本 $e$ 在区间 $[u,v]$ 内的归属倾向。

- $C_{u,v}(e) < 0$ → 样本 $e$ 被归入新的组 $G$（更像 $G^0$ 这一边）
- $C_{u,v}(e) > 0$ → 样本 $e$ 被归入新的组 $H$（更像 $H^0$ 这一边）
- $C_{u,v}(e) = 0$ → 样本 $e$ 归入**中间组 $I$**（既不明显偏向 $G^0$ 也不明显偏向 $H^0$，可能处于过渡状态或本身该区域甲基化模式不清晰）


$$C_{u,v}(e) = \text{"} e \text{ 偏向 } G^0 \text{ 的方向信号"} - \text{"} e \text{ 偏向 } H^0 \text{ 的方向信号"}$$

- 它把"样本 $e$ 到两个种子组的距离"简化成一个正负号相减的整数指标，用来判断这个样本更靠近哪一边——这一步就是"局部样本分类"（把每个样本分到 $G$/$H$/中间组 $I$）的具体计算规则,而距离函数 $d_{u,v}$ 本身的定义,是下一步公式(7)要讲的内容。

---

using the pseudo-distance function

$$
d_{u,v}(A,B) = \frac{\displaystyle\sum_{i=u}^{v} \max\!\big[0, \operatorname{sgn}(|\Delta_i(A,B)| - \epsilon)\big]}{v-u} - \rho \tag{7}
$$

- **这个公式在做什么：** 这是公式(6)里用到的"伪距离函数"$d_{u,v}$ 的具体定义。它衡量的是"$A$ 和 $B$ 在区间 $[u,v]$ 内，到底有多不一样"——但它不是常见的欧氏距离之类的连续度量，而是**先数出'差异明显超标的CpG位点个数'占区间总长度的比例，再减去一个偏置量**，得到一个可正可负的"伪距离"分数。之所以叫"伪距离"，是因为它不满足严格数学意义上距离函数的所有性质，只是借用了"距离越大越不像"这个直觉。

- **$A, B$**：要比较的两个对象——可以是"一个样本 vs 一个组"（如公式6里的 $d_{u,v}(e,G^0)$），也可以是"一个组 vs 一个组"（后面监督模式公式8里会用到 $d_{u,v}(T,G^0)$）。

- **$[u,v]$**：正在评估的区间（延续前面的定义）。

- **$i$**：区间内遍历的CpG位点索引，从 $u$ 到 $v$。

- **$\Delta_i(A,B)$**：延续公式(2)的定义，位点 $i$ 处 $A$ 和 $B$ 的平均差异信号，只是这里把公式(2)里的 $G,H$ 换成了更一般的 $A,B$。

- **$|\Delta_i(A,B)|$**：取绝对值，只关心差异的大小，不关心方向。

- **$\epsilon$**：一个用户设定的**甲基化差异阈值**（默认0.5）——用来判断"这个位点的差异算不算大"的分界线。

- **$|\Delta_i(A,B)| - \epsilon$**：拿"实际差异"减去"阈值"，结果为正就说明这个位点差异**超过了**阈值，为负就说明**没超过**。

- **$\operatorname{sgn}(\cdot)$**：符号函数，只取上面这个差值的正负号：超过阈值 → $+1$；没超过 → $-1$；正好等于 → $0$。

- **$\max[0, \operatorname{sgn}(\cdot)]$**：把符号函数的结果和 $0$ 取最大值——这一步的作用是把 $-1$ 全部**截断成 $0$**，只保留 $+1$（或 $0$）。效果上，这一整项等价于一个"指示函数"：**只要这个位点差异超过阈值，就记1分；没超过，就记0分**。

- **$\sum_{i=u}^{v}[\cdot]$**：把区间 $[u,v]$ 内所有位点的这个"0/1得分"加起来，等于**区间内差异明显超过阈值 $\epsilon$ 的CpG位点总数**。

- **$v-u$**：区间长度（位点数）。

- **$\dfrac{\sum(\cdot)}{v-u}$**：用"超标位点数"除以"区间总长度"，得到**超标位点占比**——一个介于0和1之间的比例，表示"这段区域里有多大比例的CpG，$A$和$B$差异明显"。

- **$\rho$**：一个用户设定的**偏置项**（默认0.5），起"平移/门槛"的作用。

- **减去 $\rho$**：最终结果 $d_{u,v}(A,B)$ = 超标位点占比 $- \rho$。

where $\epsilon$ is a user-defined methylation threshold (default: 0.5) and $\rho$ (default: 0.5) is a user-defined bias. In brief, the function $d$ sums up the number of cytosines with an absolute mean difference between $A$, e.g., the sample $e$, and $B$, e.g., the set $G$, that is strictly larger than $\epsilon$. This number is normalized to the size of the interval. The bias term $\rho \ge 0$ shifts this value, making it harder to achieve positive distances. In turn, for the calculation of the cluster index $C$, higher values of $\rho$ shrink the decision boundary, making the classification more conservative and assigning more samples to none of the two groups ($C_{u,v} = 0$). Consequently, for $C_{u,v}(e) < 0$, $e$ will be added to a new group $G$; for $C_{u,v}(e) > 0$, $e$ will be added to a new group $H$. For $C_{u,v}(e) = 0$, samples are assigned to the intermediate group $I$. Within the segment $[u,v]$, the two-dimensional KS (2D-KS) test for equality is applied to $G$ and $H$. The determined groups also play a critical role for the DMTree explained below.

### Supervised multi-group segmentation

In the supervised multi-group segmentation mode, metilene³ performs the steps outlined in the basic segmentation strategy, i.e., the Z-score calculations for all possible combinations of groups. Once the pair of groups $G^0$ and $H^0$, as well as the region $[u,v]$, that maximize the Z-score are found, metilene³ clusters the remaining groups with the group-version cluster index $C_{u,v}(T)$:

$$
C_{u,v}(T) = \operatorname{sgn}\!\big(d_{u,v}(T, G^0)\big) - \operatorname{sgn}\!\big(d_{u,v}(T, H^0)\big) \tag{8}
$$

Three new groups will be first initialized: $G \leftarrow G^0$, $H \leftarrow H^0$, $I \leftarrow \emptyset$. For $C_{u,v}(T) < 0$, $T$ will be merged with the group $G$. For $C_{u,v}(T) > 0$, $T$ will be merged with the group $H$. For $C_{u,v}(T) = 0$, $T$ will be merged with the group $I$. The new groups $G$ and $H$ will be used for 2D-KS test within $[u,v]$.

In addition to the 2D-KS test, Kruskal–Wallis test on the mean methylation level of DMRs will be performed across all groups using scipy⁴⁸, with a default $p$-value threshold 0.01.

### DMTree construction

We recall that metilene³ has locally assigned samples to one of three sets during each recursion step of the unsupervised segmentation. Thus, at the end of the segmentation, for some DMR $k$, there is set $G_k$, $H_k$ and $I_k$. In the following, we will use these local assignments to build the differential methylation tree (DMTree) by performing recursive binary splits. In the following we assume that some set $K$ holds all DMRs. For the DMR $k \in K$, the set $G_k$ holds the hypo-methylated samples and the set $H_k$ holds the hyper-methylated samples. Also, let $\sigma: \{1, \dots, |E|\} \to E$ be an indexing function that returns the $i$-th element of the sample set $E$.

The algorithm recursively splits the set of samples into smaller group (Box 1). In each step of the recursion and separately for each DMR, the group assignments made for individual samples are encoded in a sample assignment vector $\mathbf{v}$. A sample is encoded with $-1$ or $1$ if it belongs to the largest of the two sets $G_k$ or $H_k$; with $0$, otherwise. If a sufficient number of samples ($N_{\min}$) are encoded with either $-1$ (or $1$) or $0$, the methylation difference $S(G_k, H_k)$ is added to a hash map using the assignment pattern $\mathbf{v}$ as a key. In this way, the methylation differences associated with a particular sample assignment pattern are summed up over all available DMRs. Subsequently, the pattern with the largest cumulative difference $\mathbf{u}$ is sought. If the total difference underscores the threshold $W_{\min}$, the samples are reported as leaves of the current branch. Otherwise, the samples are split using the pattern encoded in $\mathbf{u}$ and the recursion continues. During the computation, metilene³ stores the assignment patterns and the associated DMRs to allow the investigation of all DMRs that gave rise to the respective splits in the DMTree.

### Simulation

In this study, we modified the methylome simulation framework from our previous work⁵ to extend its applicability to multi-condition. We simulated two CpG methylation backgrounds, each consisting of 50 samples, on human Chromosome 10 as previously described in ref. 5. These samples were classified into five groups: G1 ($n=10$), G2 ($n=10$), G3 ($n=10$), G4 ($n=10$), G5 ($n=10$). Based on the five aforementioned groups, we introduced five types of DMRs, totaling 1000 DMRs, and the methylation rates $p$ within the DMRs were re-simulated under four conditions (Fig. 2a):

condition A:

$$
p \sim \operatorname{Beta}(\beta,\alpha)\cdot c + \operatorname{Beta}(\alpha,\beta)\cdot(1-c) \tag{9}
$$

condition B:

$$
p \sim \operatorname{Beta}(\alpha,\beta)\cdot c + \operatorname{Beta}(\beta,\alpha)\cdot(1-c) \tag{10}
$$

condition C:

$$
p \sim \operatorname{Beta}\!\left(\frac{\lceil \alpha+\beta \rceil}{2}, \frac{\lceil \alpha+\beta \rceil}{2}\right) \tag{11}
$$

condition D:

$$
p \sim \operatorname{Beta}(\beta,\alpha)\cdot[r\cdot(c-0.5)+0.5] + \operatorname{Beta}(\alpha,\beta)\cdot[0.5-r\cdot(c-0.5)] \tag{12}
$$

here, $\alpha$ and $\beta$ are the parameters of beta distribution, $r$ is the scaling factor and $c$ is the mixing factor. We used $\alpha=40$ and $\beta=3$ (or $\alpha=3$ and $\beta=40$, randomly with a 50% probability) for background 1 and $\alpha=15$ and $\beta=5$ (or $\alpha=5$ and $\beta=15$, randomly with a 50% probability) for background 2, and randomly assign $c$ from $[0.6, 0.73, 0.87, 1]$ to a simulated DMR.

In each DMR type, each group uses the conditions in Table 1. Unlike the other DMR types, for DMR type V, the re-simulation of the CpG methylation rates in G1 started from the $fl^{\text{th}}$ CpG in the DMR and stopped at the $(N-fl)^{\text{th}}$ CpG where $N$ is the number of CpGs in the DMR and $fl$ is a random number from 6 to 15.

### Purity simulation

To mimic the tumour microenvironment, we assigned each sample a simulated purity value drawn from a normal distribution with mean 0.86 and variance 0.11, derived from the TCGA-LGG + GBM cohort⁴⁹. The final methylation rate is calculated as:

$$
\beta_i = \beta_{tumor_i} + purity_i \cdot \beta_{ctr_i} \tag{13}
$$

where $\beta_{ctr_i}$ is the background (without DMR) methylation values of sample $i$. To identify the minimal purity required for robust DMTree clustering, we further simulated datasets with purity mean of 0.8 to 0.5 and variance of 0.1, and evaluated the average Rand index of DMTree clustering results on these datasets.

### Confounding factor simulation

We introduced two confounding factors, $C_1$ and $C_2$ to increase the heterogeneity within subgroups. $C_1$ is drawn from $\operatorname{unif}(0,1)$ and $C_2$ is drawn from $N(0.5, 0.1)$. For each confounding factor, 500 regions were simulated to be correlated with $C_1$ or $C_2$. For each region $j$ correlated with $C$ ($C \in \{C_1, C_2\}$), the methylation rate $p_{i,j}$ for sample $i$ is:

$$
p_{i,j} \sim \operatorname{Beta}(\alpha,\beta)\cdot s_{i,j} + \operatorname{Beta}(\beta,\alpha)\cdot(1-s_{i,j}) \tag{14}
$$

where

$$
s_{i,j} \sim N\big(\rho\cdot(C_i-0.5)+0.5,\ 0.3\big) \tag{15}
$$

where $\rho$ controls the correlation coefficient between the methylation of the region and the confounding factor:

$$
\rho \sim \operatorname{unif}(0.1, 1) \tag{16}
$$

Similar to the DMR simulation, $\alpha=40$ and $\beta=3$ (or $\alpha=3$ and $\beta=40$, randomly with a 50% probability) were used for background 1 and $\alpha=15$ and $\beta=5$ (or $\alpha=5$ and $\beta=15$, randomly with a 50% probability) were used for background 2.

### Batch effect simulation

We also acknowledged that there might be a batch effect in bisulfite sequencing. In the simulation, 25 samples were assigned to batch 1, and the others were assigned to batch 2. 500 CpGs were randomly selected to be unmethylated in batch 1 and another 500 CpGs to be unmethylated in batch 2.

### Parameter of the metilene³

In the supervised mode, we keep minimal mean difference of a DMR $\delta$ to 0.1 as the original metilene⁵. In the unsupervised mode, we tested different combinations of the parameters: minimal mean difference of a DMR $\delta \in \{0.1, 0.25, 0.5\}$, minimal number of samples in a cluster $N_{\min}$ from 2 to 6, and minimal cumulative difference $W_{\min} \in \{1, 2, 5, 10, 20\}$. To identify an optimal parameter set, we performed grid-search of these parameters on background 1 and background 2, and kept parameter sets that use the same minimal number of samples and led to consistent clustering results between two backgrounds. Five parameter sets were kept for the background 1 and one parameter set was kept for the background 2. The parameter sets were then ranked by $\delta$ and $W_{\min}$ in descending order and the first one was set as the default parameter of metilene³. For background 1, $\delta = 0.5$, $N_{\min} = 6$, $W_{\min} = 10$ were used, and $\delta = 0.25$, $N_{\min} = 6$, $W_{\min} = 1$ were used for background 2.

Considering that real-world datasets vary in number of samples and number of available CpGs $n_{CpG}$ (normalized to the number of simulated CpGs $n_{simulatedCpG}$), we proposed a formula to estimate the optimal $W_{\min}$ with:

$$
W' = \alpha \cdot \frac{n_{CpG}}{n_{simulatedCpG}} \tag{17}
$$

where $\alpha$ is a coefficient to take the number of samples in the real-world dataset $n$ into account by comparing to the number of samples simulated $n_{sim}$ ($n_{sim} = 50$):

$$
\alpha = 1 + \left[\frac{n-1}{n_{sim}-1}\right]\cdot (W_{default} - 1) \tag{18}
$$

Here $W_{default}$ would be 10 for the background 1 and 1 for the background 2. In practice, metilene³ will first run with the first parameter set (optimal for background 1) and if no cluster can be found, the second parameter set (optimal for background 2) will be used. In addition, to simplify the parameter, $W_{\min}$ will be set to 1 if $W' < 5$ and will be set to 100 if $W' \ge 50$, otherwise, $W_{\min}$ will be set to 10.

### Benchmarking on DMR identification

We tested metilene³, wgbstools (version 0.2.2)¹³, SMART (version 2.2.8)¹² and methylscore¹⁵ using simulated datasets on a Linux computing cluster. All DMR callers were executed with default parameters, except for the number of threads. Sensitivity, precision and Jaccard index were calculated as follows:

$$
\text{sensitivity} = \frac{\sum Length_{overlapped}}{\sum Length_{simulated}} \tag{19}
$$

$$
\text{precision} = \frac{\sum Length_{overlapped}}{\sum Length_{predicted}} \tag{20}
$$

$$
\text{Jaccard} = \frac{\sum Length_{overlapped}}{\sum Length_{predicted} + \sum Length_{simulated} - \sum Length_{overlapped}} \tag{21}
$$

The overlaps between predicted DMRs and simulated DMRs were computed using pybedtools⁵⁰.

For local clustering, accuracy is defined as the proportion of DMRs for which the predicted group (from the predicted DMR with the largest overlap with the simulated DMR) matches the simulated group.

For time and memory usage benchmarking, DMR calling was performed using both a single core and ten cores. Since SMART does not offer a parameter to control the number of threads, it was only tested with ten cores. The workflow of methylscore could not be executed in a single run, and thus the time and memory usage could not be benchmarked.

The unsupervised mode of metilene³ was tested on the simulated datasets using default parameters.

### Benchmarking on clustering

We compared DMTree to unsupervised clustering methods commonly used in DNA methylation data analysis, including hierarchical clustering and K-means clustering⁵¹. We first applied DMTree, hierarchical clustering, and K-means clustering based on the first two principal components to the simulated dataset. To further measure the robustness, and the sensitivity of small group detection, we stochastically sampled $n$ ($n$ from one to ten) samples from group G1 and $40-n$ samples from other groups, and repeated for 50 times. For hierarchical clustering and K-means clustering, we evaluated the results under three to seven clusters. For metilene³, as the number of clusters is automatically determined, we evaluated the results under different minimal number of samples in a cluster, $N_{\min}$, from two to six. The performance of the clustering methods was measured by the average Rand index using scikit-learn⁵², compared to the simulated labels. The sensitivity of small group detection was measured by the average Rand index among samples from G1 and G2, as G2 is the most similar group to G1.

---

# 内容总结

metilene³是原有DMR（差异甲基化区域）检测工具metilene的第三代版本，由马克斯·普朗克分子遗传学研究所、莱布尼茨老龄化研究所等团队开发，2026年发表于*Nature Communications*。它要解决的核心问题是：现有的DMR检测方法几乎都要求预先知道样本属于哪个分组（例如"肿瘤 vs. 正常"），但在很多真实场景（尤其是临床数据）中，分组信息是未知或不可靠的，而且大多数工具只能做两组比较，难以处理三组及以上的多条件比较。

metilene³给出的方案是一套统一的算法框架，同时支持两种运行模式：

1. **监督模式**：与传统DMR分析类似，用户提供样本的分组标签，metilene³在此基础上把原有metilene的环形二元分割（CBS）算法从两组推广到任意多组之间的比较。
2. **无监督（自动分类）模式**：不需要任何分组标签，metilene³在基因组的每个局部区段独立地对样本进行"局部分类"——先找出该区段甲基化差异最大的一对样本作为两个分组的种子，再根据其余样本与这两个种子的甲基化距离，把它们分别归入"低甲基化组G"、"高甲基化组H"，或无法明确归类的"中间组I"。这一步在全基因组的每个候选区段上重复进行，因此每个DMR都会产生一套局部的样本分组结果。

metilene³的第二个核心创新，是把这些成千上万个"局部分组结果"汇聚起来，构建一棵被称为**DMTree（差异甲基化树）**的二叉聚类树：算法寻找在众多DMR中反复出现、且累积甲基化差异（ΣWᵢ）最大的分组模式，将其作为树的第一次分裂，然后对每个分支递归地重复这一过程，直至无法再分裂（样本量或累积差异低于阈值）为止。由于只有"多次重复出现、差异幅度大"的分组模式才会主导树的构建，DMTree天然具有抗噪声、抗批次效应的能力，能够把真正稳定的生物学分组同随机噪声或混杂因素区分开来。聚类完成后，metilene³还会以这些聚类簇为分组，再跑一次监督模式的DMR检测，得到更完整的DMR集合。

在方法学验证方面，作者在模拟数据（人类10号染色体、50个样本、5种分组、5类DMR、并引入批次效应/肿瘤纯度/混杂因子等干扰因素）上，将metilene³与wgbstools、SMART2（监督模式对比）和methylscore（无监督模式对比）进行了系统比较，结果显示metilene³在灵敏度、Jaccard相似度（与真实模拟DMR的一致性最高达0.96）、DMR长度准确性等指标上均显著优于对比工具；在聚类稳健性上，DMTree的平均Rand指数持续超过0.9，明显优于层次聚类（<0.7）和K-均值聚类（<0.4），且对小群体（低至5%的样本占比）的识别也远比传统聚类敏感；计算效率上，metilene³在单核/多核条件下均能以标准笔记本电脑量级的时间和内存（如10核处理35个样本约10分钟、不到1GB内存）完成分析。

在真实生物学数据上的应用则展示了该工具的价值：在正常人体细胞图谱（205个样本、77种细胞类型）中，DMTree自动重现了已知的细胞类型和谱系关系（如准确区分胰岛的α/β/δ细胞）；在血液细胞谱系数据中，DMTree自动重建出髓系-淋系分化树，并通过motif富集与单细胞RNA-seq数据交叉验证，确认了C/EBPβ与LEF1分别驱动髓系和T细胞谱系分化节点；在胶质母细胞瘤数据中，DMTree正确按IDH突变状态和肿瘤/正常状态分层，并识别出一例转录组特征异常的"离群"肿瘤样本；在胰腺导管腺癌（PDAC）数据中，DMTree清晰重现了从正常腺泡/导管细胞、经癌前病变（PanIN）到PDAC的转化轨迹，并发现PDAC特异性低甲基化区域显著富集NFKB2和NFATc1转录因子结合基序，二者常共同出现在同一DMR中，其邻近基因（如MDK）在肿瘤中表达逐步上调，提示一条潜在的"NF-κB–NFAT–MDK"信号轴可能参与胰腺癌的发生发展。

总体而言，这篇论文提供的是一个"检测DMR + 无监督分类样本 + 可解释聚类树 + 下游功能注释（motif、GSEA）"一体化的分析框架，其方法学贡献主要体现在：把经典的两组CBS分割算法推广到多组/无标签场景，并首创性地把"逐区段局部分类结果的累积模式"作为构建可解释聚类树的依据，从而在不依赖预先分组、且计算高效的前提下，实现了稳健的样本分层和DMR发现。

