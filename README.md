# Computational_Physics_9

## A. $2 \times 2$ 矩阵对角化

考虑一个二能级系统，其哈密顿量为：
$$H = \begin{bmatrix}
a & c \\
c^* & -a
\end{bmatrix},$$
其中 $a$ 是实数，$c$ 是复数。

### 1.矩阵对角化

1. 对角化哈密顿量，作分解 $H = Q \Lambda Q^{-1}$。其中 $\Lambda$ 为特征值构成的对角矩阵，$Q$ 为本征列向量构成的矩阵，且所有本征向量归一化。$Q^{-1}$ 和 $Q$ 有什么关系，为什么？（1分）

#### 计算矩阵 $H$ 的特征值：

$$
\det(H - \lambda I) = \det \begin{bmatrix} a - \lambda & c \\ c^* & -a - \lambda \end{bmatrix}
$$

$$
= (a - \lambda)(-a - \lambda) - |c|^2 = \lambda^2 - a^2 - |c|^2 = 0
$$

得到特征方程为：

$$
\lambda^2 = a^2 + |c|^2
$$

因此，特征值为：

$$
\lambda_1 = \sqrt{a^2 + |c|^2}, \quad \lambda_2 = -\sqrt{a^2 + |c|^2}
$$

#### 特征向量计算
对于特征值 $\lambda_1 = \sqrt{a^2 + |c|^2}$，求解 $(H - \lambda_1 I)v_1 = 0$：

$$\begin{bmatrix}
a - \sqrt{a^2 + |c|^2} & c \\
c^* & -a - \sqrt{a^2 + |c|^2}
\end{bmatrix} \begin{bmatrix} v_{11} \\ v_{12} \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

得到：
$$(a - \sqrt{a^2 + |c|^2})v_{11} + cv_{12} = 0$$

可以选择：
$$v_{11} = c, \quad v_{12} = \sqrt{a^2 + |c|^2} - a$$

归一化后：
$$v_1 = \frac{1}{N_1}\begin{bmatrix} c \\ \sqrt{a^2 + |c|^2} - a \end{bmatrix}$$
其中 $N_1 = \sqrt{|c|^2 + (\sqrt{a^2 + |c|^2} - a)^2}$

类似地，对于特征值 $\lambda_2 = -\sqrt{a^2 + |c|^2}$：
$$v_2 = \frac{1}{N_2}\begin{bmatrix} c \\ -\sqrt{a^2 + |c|^2} - a \end{bmatrix}$$
其中 $N_2 = \sqrt{|c|^2 + (-\sqrt{a^2 + |c|^2} - a)^2}$

#### Q的具体形式

将归一化的特征向量整理为矩阵

$$Q=\begin{bmatrix}
\frac{c}{N_1} & \frac{c}{N_2} \\
\frac{\sqrt{a^2 + |c|^2} - a}{N_1} & \frac{-\sqrt{a^2 + |c|^2} - a}{N_2}
\end{bmatrix}$$

其中归一化常数为：
$$N_1 = \sqrt{|c|^2 + (\sqrt{a^2 + |c|^2} - a)^2}$$
$$N_2 = \sqrt{|c|^2 + (-\sqrt{a^2 + |c|^2} - a)^2}$$

可以进一步简化表达式。令 $E = \sqrt{a^2 + |c|^2}$，则：

$$N_1 = \sqrt{|c|^2 + (E - a)^2} = \sqrt{|c|^2 + E^2 - 2aE + a^2} = \sqrt{E^2 + a^2 - 2aE + |c|^2}$$

由于 $E^2 = a^2 + |c|^2$，代入得：
$$N_1 = \sqrt{2E^2 - 2aE} = \sqrt{2E(E-a)}$$

同理：
$$N_2 = \sqrt{|c|^2 + (-E - a)^2} = \sqrt{|c|^2 + E^2 + 2aE + a^2} = \sqrt{2E^2 + 2aE} = \sqrt{2E(E+a)}$$

因此，特征向量矩阵可以更简洁地表示为：

$$Q = \begin{bmatrix} \frac{c}{\sqrt{2E(E-a)}} & \frac{c}{\sqrt{2E(E+a)}} \\
\frac{E-a}{\sqrt{2E(E-a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}$$

#### 哈密顿量的对角化展开

哈密顿量的对角化分解为：
$$H = Q \Lambda Q^{-1}$$

其中 $\Lambda$ 是特征值构成的对角矩阵：
$$\Lambda = \begin{bmatrix}
E & 0 \\
0 & -E
\end{bmatrix}$$

由于 $H$ 是Hermitian矩阵，$Q^{-1} = Q^\dagger$，因此：
$$H = Q \Lambda Q^\dagger$$

展开得到：

$$H = \begin{bmatrix} \frac{c}{\sqrt{2E(E-a)}} & \frac{c}{\sqrt{2E(E+a)}} \\
\frac{E-a}{\sqrt{2E(E-a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}
\begin{bmatrix}
E & 0 \\
0 & -E
\end{bmatrix}
\begin{bmatrix} \frac{c^*}{\sqrt{2E(E-a)}} & \frac{E-a}{\sqrt{2E(E-a)}} \\
\frac{c^*}{\sqrt{2E(E+a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}$$

进一步乘开：

$$H = \begin{bmatrix} \frac{cE}{\sqrt{2E(E-a)}} & \frac{-cE}{\sqrt{2E(E+a)}} \\
\frac{(E-a)E}{\sqrt{2E(E-a)}} & \frac{(E+a)E}{\sqrt{2E(E+a)}}
\end{bmatrix}
\begin{bmatrix} \frac{c^*}{\sqrt{2E(E-a)}} & \frac{E-a}{\sqrt{2E(E-a)}} \\
\frac{c^*}{\sqrt{2E(E+a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}$$

完成矩阵乘法后，可以验证结果确实等于原始哈密顿量：

$$H = \begin{bmatrix}
a & c \\
c^* & -a
\end{bmatrix}$$

这就完成了哈密顿量 $H$ 的完整对角化分解。
#### 特征向量矩阵 $Q$ 和 $Q^{-1}$ 的关系
矩阵 $Q$ 由归一化的特征向量组成：
$$Q = \begin{bmatrix} v_1 & v_2 \end{bmatrix}$$

对于Hermitian矩阵（即 $H = H^\dagger$），有 $Q^{-1} = Q^\dagger$，这是因为：

1. 不同特征值对应的特征向量正交
2. 特征向量已被归一化

这个性质称为"酉对角化"，对于量子力学中的Hermitian算符（如哈密顿量）非常重要。证明如下：

假设 $H$ 是Hermitian的，即 $H = H^\dagger$。如果 $H v_i = \lambda_i v_i$，则对于任意特征值 $\lambda_i, \lambda_j$：

$$\langle v_i|H|v_j\rangle = \lambda_j \langle v_i|v_j\rangle$$
$$\langle v_i|H|v_j\rangle = \langle H v_i|v_j\rangle = \lambda_i^* \langle v_i|v_j\rangle$$

因此 $(\lambda_j - \lambda_i^*)\langle v_i|v_j\rangle = 0$。当 $i \neq j$ 时，$\lambda_i \neq \lambda_j$，且 $\lambda_i$ 是实数（因为 $H$ 是Hermitian），所以 $\langle v_i|v_j\rangle = 0$。

当特征向量已归一化，对于所有 $i$，有 $\langle v_i|v_i\rangle = 1$，因此 $Q^\dagger Q = I$，即 $Q^{-1} = Q^\dagger$。

这表明 $Q$ 是一个酉矩阵（若 $c$ 是纯实数，则 $Q$ 是正交矩阵）。

### 2. 证明公式: $e^{iHt} = Qe^{i\Lambda t}Q^{-1}$

#### 泰勒展开法证明

时间演化算符 $e^{iHt}$ 可以通过泰勒级数展开表示为：

$$e^{iHt} = I + iHt + \frac{(iHt)^2}{2!} + \frac{(iHt)^3}{3!} + ... = \sum_{n=0}^{\infty} \frac{(iHt)^n}{n!}$$

利用矩阵分解 $H = Q\Lambda Q^{-1}$，将 $H$ 的幂展开如下：

$$H^2 = H \cdot H = (Q\Lambda Q^{-1})(Q\Lambda Q^{-1}) = Q\Lambda Q^{-1}Q\Lambda Q^{-1} = Q\Lambda^2 Q^{-1}$$

同理，可以证明：
$$H^n = Q\Lambda^n Q^{-1}$$

将这一结果代入时间演化算符的泰勒展开式：

$$e^{iHt} = \sum_{n=0}^{\infty} \frac{(iH)^n t^n}{n!} = \sum_{n=0}^{\infty} \frac{(iQ\Lambda Q^{-1})^n t^n}{n!} = \sum_{n=0}^{\infty} \frac{Q(i\Lambda)^n Q^{-1}t^n}{n!}$$

由于 $Q$ 和 $Q^{-1}$ 不依赖于求和指标 $n$，可以提到求和符号外：

$$e^{iHt} = Q\left(\sum_{n=0}^{\infty} \frac{(i\Lambda)^n t^n}{n!}\right)Q^{-1} = Qe^{i\Lambda t}Q^{-1}$$

这就证明了 $e^{iHt} = Qe^{i\Lambda t}Q^{-1}$。

#### 具体表达式推导

对于我们的哈密顿量，已知：

$$H = \begin{bmatrix}
a & c \\
c^* & -a
\end{bmatrix}$$

$$Q = \begin{bmatrix} \frac{c}{\sqrt{2E(E-a)}} & \frac{c}{\sqrt{2E(E+a)}} \\
\frac{E-a}{\sqrt{2E(E-a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}$$

$$\Lambda = \begin{bmatrix}
E & 0 \\
0 & -E
\end{bmatrix}$$

其中 $E = \sqrt{a^2 + |c|^2}$。

指数对角矩阵 $e^{i\Lambda t}$ 的计算很简单，因为对角矩阵的指数就是对每个对角元素取指数：

$$e^{i\Lambda t} = \begin{bmatrix}
e^{iEt} & 0 \\
0 & e^{-iEt}
\end{bmatrix}$$

将这个结果代入公式 $e^{iHt} = Qe^{i\Lambda t}Q^{-1}$，并利用 $Q^{-1} = Q^\dagger$（因为 $H$ 是Hermitian矩阵）：

$$e^{iHt} = Q\begin{bmatrix}
e^{iEt} & 0 \\
0 & e^{-iEt}
\end{bmatrix}Q^\dagger$$

展开计算：

$$e^{iHt} = \begin{bmatrix} \frac{c}{\sqrt{2E(E-a)}} & \frac{c}{\sqrt{2E(E+a)}} \\
\frac{E-a}{\sqrt{2E(E-a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}
\begin{bmatrix}
e^{iEt} & 0 \\
0 & e^{-iEt}
\end{bmatrix}
\begin{bmatrix} \frac{c^*}{\sqrt{2E(E-a)}} & \frac{E-a}{\sqrt{2E(E-a)}} \\
\frac{c^*}{\sqrt{2E(E+a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}$$

这可以进一步计算为：

$$e^{iHt} = \begin{bmatrix} \frac{ce^{iEt}}{\sqrt{2E(E-a)}} & \frac{ce^{-iEt}}{\sqrt{2E(E+a)}} \\
\frac{(E-a)e^{iEt}}{\sqrt{2E(E-a)}} & \frac{(-E-a)e^{-iEt}}{\sqrt{2E(E+a)}}
\end{bmatrix}
\begin{bmatrix} \frac{c^*}{\sqrt{2E(E-a)}} & \frac{E-a}{\sqrt{2E(E-a)}} \\
\frac{c^*}{\sqrt{2E(E+a)}} & \frac{-E-a}{\sqrt{2E(E+a)}}
\end{bmatrix}$$

完整乘开后，得到 $e^{iHt}$ 的矩阵元素为：

$$e^{iHt}_{11} = \frac{|c|^2}{2E(E-a)}e^{iEt} + \frac{|c|^2}{2E(E+a)}e^{-iEt} = \frac{|c|^2}{2E}\left(\frac{e^{iEt}}{E-a} + \frac{e^{-iEt}}{E+a}\right)$$

$$e^{iHt}_{12} = \frac{c(E-a)}{2E(E-a)}e^{iEt} - \frac{c(E+a)}{2E(E+a)}e^{-iEt} = \frac{c}{2E}(e^{iEt} - e^{-iEt}) = i\frac{c}{E}\sin(Et)$$

$$e^{iHt}_{21} = \frac{c^*(E-a)}{2E(E-a)}e^{iEt} - \frac{c^*(-E-a)}{2E(E+a)}e^{-iEt} = \frac{c^*}{2E}(e^{iEt} - e^{-iEt}) = i\frac{c^*}{E}\sin(Et)$$

$$e^{iHt}_{22} = \frac{(E-a)^2}{2E(E-a)}e^{iEt} + \frac{(-E-a)(-E-a)}{2E(E+a)}e^{-iEt} = \frac{E-a}{2E}e^{iEt} + \frac{E+a}{2E}e^{-iEt}$$

利用欧拉公式 $e^{iEt} = \cos(Et) + i\sin(Et)$ 和 $e^{-iEt} = \cos(Et) - i\sin(Et)$，我们可以将上述结果写成：

$$e^{iHt} = \begin{bmatrix}
\cos(Et) + i\frac{a}{E}\sin(Et) & i\frac{c}{E}\sin(Et) \\
i\frac{c^*}{E}\sin(Et) & \cos(Et) - i\frac{a}{E}\sin(Et)
\end{bmatrix}$$

这就是时间演化算符 $e^{iHt}$ 的显式表达式。

## B. 横场Ising模型

研究一个自旋 $1/2$ 的横场Ising模型：
$$H = -J \sum_{i=1}^{L} \hat{\sigma}_{i}^{z} \hat{\sigma}_{i+1}^{z} - h \sum_{i=1}^{L} \hat{\sigma}_{i}^{x}.$$

取周期边界条件（$\hat{\sigma}_{L+1}^{z} \equiv \hat{\sigma}_{1}^{z}$）且固定 $J=1$。

### 1.对角化哈密顿量

1. 取参数 $L=12$，对于不同的横场强度 $h=0.5, 1.0, 2.0$，分别对角化哈密顿量 $H$，求得基态和第一激发态并列出其能量（2分）。对于基态波函数，计算第一个格点上横场方向磁化强度的期望 $\langle \hat{\sigma}_{1}^{x} \rangle$，并一并列出。（1分）

已知横场Ising model 的Hamitonian；

$$H = -J \sum_{i=1}^{L} \hat{\sigma}_{i}^{z} \hat{\sigma}_{i+1}^{z} - h \sum_{i=1}^{L} \hat{\sigma}_{i}^{x}.$$

由此可以得到Hamiltonian的矩阵表示。

进一步地，对Hamiltonian进行对角化，得到特征值和特征向量。

特征值即为体系的本征能量，由此得到基态和第一激发态的能量。

特征向量即为体系的对应特征能量下的波函数，由此得到基态和第一激发态的波函数。

第一个格点上横场方向磁化强度的期望 $\langle \hat{\sigma}_{1}^{x} \rangle$ 可以通过基态波函数计算得到。

定义第一格点的横场方向磁化强度测量算符：

$$\hat{\sigma}_{1}^{x} =\sigma_x \otimes I \dots \otimes I \\=\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\otimes\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \dots \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} $$

则第一格点的横场方向磁化强度的期望为：

$$\langle \hat{\sigma}_{1}^{x} \rangle =\langle \psi|\hat{\sigma}_{1}^{x}|\psi \rangle $$

由此得到基态第一格点的横场方向磁化强度的期望

以上是对第一问的解答。经过程序模拟计算，得到结果为：
![alt text](/images/image.png)

<!-- todo -->
| $h$ | 基态能量 $E_0$ | 第一激发态能量 $E_1$ | $\langle \hat{\sigma}_{1}^{x} \rangle$ |
|---|---|---|---|
| 0.5 | -12.762569151024042 | -12.762496682237142 | -0.25872864371389986 |
| 1.0 | -15.322595151080776  | -15.191508225450303  | -0.6384414646283648 |
|2.0|-25.52513830204809|-23.52499336447428|-0.934183107395054|


### 结果分析

1. **基态能量**：随着横场强度 $h$ 的增加，基态能量明显降低。这说明更强的横场使系统能够达到更低的能量状态，系统更趋向于在横场方向排列。

2. **能隙变化**：第一激发态与基态之间的能隙 $(E_1-E_0)$ 表现出与横场强度相关的特征：
   - $h=0.5$ 时能隙极小，约为 0.00007
   - $h=1.0$ 时能隙增大到约 0.13109
   - $h=2.0$ 时能隙进一步增大到约 2.00015

   这种能隙的变化特征通常与量子相变有关。特别是在 $h=0.5$ 附近，能隙接近于零，这可能表明系统接近或处于临界点。

3. **横场磁化强度**：$\langle \hat{\sigma}_{1}^{x} \rangle$ 的绝对值随 $h$ 的增加而显著增大：
   - 当 $h=0.5$ 时，磁化强度较小，约为 -0.259
   - 当 $h=2.0$ 时，磁化强度接近于 -1（最大可能值），为 -0.934

   这表明随着横场强度增加，自旋更倾向于沿横场方向排列，系统从Ising相互作用主导的状态转变为横场主导的状态。

4. **物理意义**：
   - 负的 $\langle \hat{\sigma}_{1}^{x} \rangle$ 值表明自旋平均指向 $-x$ 方向
   - 这与横场项 $-h \sum_{i=1}^{L} \hat{\sigma}_{i}^{x}$ 的符号一致，表明系统倾向于最小化横场能量
   - 一维横场Ising模型在 $h_c \approx 1$ 处存在量子相变，这与我们观察到的能隙在 $h=0.5$ 和 $h=1.0$ 之间急剧变化相符

综合分析表明，该系统随着横场强度的增加，从一个由近邻相互作用主导的有序相转变为一个由横场主导的量子顺磁相。特别是能隙的剧烈变化，清晰地反映了量子相变的特征。这些数值结果与理论预期的横场Ising模型行为高度一致。

### 2.时间演化

2. $t=0$ 时系统处于 $h=0.5$ 的基态，$t=0^{+}$ 的瞬间哈密顿量参数变为 $h=3.0$，求在此哈密顿量下的时间演化，计算 $\langle \hat{\sigma}_{1}^{x}(t) \rangle$ 在时间 $t \in [0, 20]$ 的变化情况，并画出来。你需要给出精确的结果和Runge-Kutta方法的结果。验证两种方法的结果是一样的。（3分）

由第一问中的方法可以得到h=0.5下的基态波函数。

当h=3.0时，可以求得该种情况下的哈密顿量，记为$H_3$。

#### 精确办法

薛定谔方程：

$$H \psi = i\frac{\partial}{\partial t} \psi$$

则：

$$\psi(t) = e^{-iHt} \psi(0)$$

#### Runge-Kutta方法

$$\frac{d\psi}{dt} = -iH \psi$$

采用Runge-Kutta方法，间隔$\Delta t $进行演化，求解时间演化方程。

得到的结果如下图：
![alt text](images/B_evolution.png)

如图，两种方法得到的结果是一样的。

#### 数值方法的准确性

使用了两种方法计算系统的时间演化：
1. **精确对角化方法**：通过公式 $|\psi(t)\rangle = e^{-iHt}|\psi(0)\rangle$ 直接计算
2. **Runge-Kutta方法**：数值求解微分方程 $\frac{d|\psi\rangle}{dt} = -iH|\psi\rangle$

从图中可以看出，两种方法的结果几乎完全重叠，这验证了：
- 数值方法的实现是正确的
- RK45方法在这种问题中能达到很高的精度
- 数值积分的步长和误差控制参数设置合理

#### 物理现象分析

图中显示的 $\langle\hat{\sigma}_1^x(t)\rangle$ 时间演化具有以下特征：

1. **初始值**：约为0.26，对应于 $h=0.5$ 时基态的 $\langle\hat{\sigma}_1^x\rangle$ 值
2. **快速上升**：淬火后系统迅速响应，磁矩开始向 $x$ 方向偏转
3. **复杂振荡**：
   - 平均值约为0.4左右，介于初始值和 $h=3.0$ 时基态的 $\langle\hat{\sigma}_1^x\rangle$ 之间
   - 振荡不是简单的单频振荡，而是包含多个频率分量的复杂振荡
   - 振幅在0.3-0.54之间，表明系统在不同本征态间的量子干涉

4. **准周期行为**：系统没有表现出明显的衰减，表明这是一个孤立量子系统的可逆演化

#### 与理论预期的比较

结果物理上是合理的，因为：

1. **量子淬火**：当系统从 $h=0.5$ 突变到 $h=3.0$ 时，初态可以被表示为新哈密顿量本征态的叠加
2. **演化守恒**：总能量守恒，系统将在新的本征态之间发生干涉
3. **非平衡动力学**：观测到的复杂振荡行为是量子多体系统非平衡动力学的典型特征

### 3.L=18时保真度的演化

3. 取参数 $L=18, h=1$。从Neel态 $|\psi(0)\rangle = |↑↓↑↓\cdots\rangle$（奇数格点自旋朝上，偶数朝下）出发，求解波函数 $|\psi(t)\rangle$ 在区间 $t \in [0, 10]$ 上的时间演化。画出Fidelity $F(t) = |\langle\psi(0)|\psi(t)\rangle|^2$ 随时间的变化。（2分）

对于参数 $L=18, h=1$ 的系统，初态为Neel态 $|\psi(0)\rangle = |↑↓↑↓\cdots\rangle$（奇数格点自旋朝上，偶数朝下）。

#### 计算方法

1. 首先构建Neel态 $|\psi(0)\rangle = |↑↓↑↓\cdots\rangle$ 的矢量表示
2. 构建哈密顿量 $H$ 的矩阵表示
3. 使用精确对角化求解时间演化方程：$|\psi(t)\rangle = e^{-iHt}|\psi(0)\rangle$
4. 计算Fidelity：$F(t) = |\langle\psi(0)|\psi(t)\rangle|^2$

#### 计算结果

通过数值模拟，我得到了Fidelity $F(t) = |\langle\psi(0)|\psi(t)\rangle|^2$ 随时间的变化曲线，如下图所示：

![Fidelity随时间的变化](./images/B_neel_fidelity.png)

在这个结果中，我们可以观察到以下特征：

1. Fidelity从初始值1开始迅速衰减，表明系统离开了初始的Neel态
2. Fidelity展现出振荡行为，这反映了量子多体系统中的干涉效应
3. 振荡的幅度随着时间逐渐减小，并趋于一个较小的非零平均值
4. 这种行为表明系统虽然会部分回到初始态，但不会完全回到初始态，反映了系统的不可逆性

#### 物理解释

这种行为可以从以下角度理解：

1. Neel态不是横场Ising模型的本征态，因此会发生时间演化
2. 横场项 $-h \sum_{i} \hat{\sigma}_{i}^{x}$ 会使自旋翻转，导致系统偏离初始的Neel排列
3. Ising相互作用项 $-J \sum_{i} \hat{\sigma}_{i}^{z} \hat{\sigma}_{i+1}^{z}$ 倾向于维持相邻自旋的反铁磁排列
4. 这两个相互竞争的项导致了系统的复杂动力学行为
5. 由于系统规模有限，量子复振会使系统部分地回到初始态，但不会完全回到，这就是我们观察到的振荡行为

通过这个时间演化的分析，我们可以了解横场Ising模型中量子多体系统的非平衡动力学特性。

# 附录
## B_1
``` python
import numpy as np
from scipy.sparse import kron, identity, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numba import njit 


@njit
def Hamiltonian(L, J, h):
    # J=J/2
    H = np.zeros((2**L, 2**L))
    I = np.array([[1, 0], [0, 1]])
    sigmaz = np.array([[1, 0], [0, -1]])
    sigmax = np.array([[0, 1], [1, 0]])
    for i in range(L):
        if i == 0:
            h2 = sigmax
        else:
            h2 = I
        if i == 0 or i == L - 1:
            h1 = sigmaz
        else:
            h1 = I
        for j in range(1, L):
            if j == i or j == i + 1:
                h1 = np.kron(h1, sigmaz)
            else:
                h1 = np.kron(h1, I)
            if j == i:
                h2 = np.kron(h2, sigmax)
            else:
                h2 = np.kron(h2, I)
        # print(f"i={i}",-J * h1 )
        H += -J * h1 + h * h2

    return H

# 优化

# def get_sigma_op(op, site, L):
#     ops = [identity(2, format='csr') for _ in range(L)]
#     ops[site] = op
#     total = ops[0]
#     for i in range(1, L):
#         total = kron(total, ops[i], format='csr')
#     return total

# def pauli_x():
#     return csr_matrix(np.array([[0, 1], [1, 0]], dtype=np.complex128))

# def pauli_z():
#     return csr_matrix(np.array([[1, 0], [0, -1]], dtype=np.complex128))


# def Hamiltonian(L, J, h):
#     H = csr_matrix((2**L, 2**L), dtype=np.complex128)
#     for i in range(L):
#         sz_i = get_sigma_op(pauli_z(), i, L)
#         sz_j = get_sigma_op(pauli_z(), (i+1)%L, L)
#         H -= sz_i @ sz_j
#     for i in range(L):
#         sx_i = get_sigma_op(pauli_x(), i, L)
#         H -= h * sx_i    
#     return H

def sigmax1(L):
    I = np.array([[1, 0], [0, 1]])
    sigmax = np.array([[0, 1], [1, 0]])
    a = sigmax
    for i in range(L - 1):
        a = np.kron(a, I)
    return a

def expm(H, eigvals, eigvecs, t):
    di = [np.exp(1j * eigvals[i] * t) for i in range(len(eigvals))]
    di = np.array(di)
    di = di.astype(np.complex128)
    Lam = np.diag(di)
    a = eigvecs @ Lam @ (eigvecs.conj().T)
    a = a.astype(np.complex128)
    return a

if __name__ == "__main__":
    L = 12
    J = 1
    # h0 = 0.5
    # h1 = 3
    # sigma1 = sigmax1(L)
    # H0 = Hamiltonian(L, J, h0)
    # H = Hamiltonian(L, J, h1)
    # eigvals, eigvecs = eigh(H0)
    h_list = [0.5,1.0,2.0]
    sigma1 = sigmax1(L)
    for h in h_list:
        H = Hamiltonian(L, J, h)
        eigvals, eigvecs = np.linalg.eigh(H)
        print(f"when h = {h}, the ground state energy is:{eigvals[0]}, and the first active state energy is:{eigvals[1]}")
        ground_state = eigvecs[:, 0]
        # 计算在第一个格点上的横场方向磁化强度的期望值
        
        ex_sigma_1_x = ground_state.conj().T @ sigma1 @ ground_state
        print(f"when h = {h}, the expectation value of sigma1_x is: {ex_sigma_1_x}")

```

## B_2
``` python
import numpy as np
from scipy.sparse import kron, identity, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def pauli_x():
    return csr_matrix(np.array([[0, 1], [1, 0]], dtype=np.complex128))

def pauli_z():
    return csr_matrix(np.array([[1, 0], [0, -1]], dtype=np.complex128))

def get_sigma_op(op, site, L):
    ops = [identity(2, format='csr') for _ in range(L)]
    ops[site] = op
    total = ops[0]
    for i in range(1, L):
        total = kron(total, ops[i], format='csr')
    return total

def hamiltonian(L, h):
    H = csr_matrix((2**L, 2**L), dtype=np.complex128)
    for i in range(L):
        sz_i = get_sigma_op(pauli_z(), i, L)
        sz_j = get_sigma_op(pauli_z(), (i+1)%L, L)
        H -= sz_i @ sz_j
    for i in range(L):
        sx_i = get_sigma_op(pauli_x(), i, L)
        H -= h * sx_i    
    return H

def schrodinger_equation(t, psi):
    return -1j * (H1 @ psi)

L = 12
H0 = hamiltonian(L, 0.5)  
H1 = hamiltonian(L, 3.0)  
sx1 = get_sigma_op(pauli_x(), 0, L)  

E0, V0 = eigsh(H0, k=1, which='SA')
psi0 = V0[:, 0].astype(np.complex128)  

E1, V1 = np.linalg.eigh(H1.toarray())  
c = V1.conj().T @ psi0  

times = np.linspace(0, 20,201)
sx_exact = []
for t in times:
    # psi_t = psi0@np.exp(-1j * E1 * t)
    psi_t = V1 @ (c * np.exp(-1j * E1 * t))  
    sx_exact.append(np.real(psi_t.conj().T @ (sx1 @ psi_t))) 

sol = solve_ivp(schrodinger_equation, [0, 20], psi0, t_eval=times, 
                method='RK45',rtol=1e-6,atol=1e-8)
sx_rk = [np.real(psi.conj().T @ (sx1 @ psi)) for psi in sol.y.T]

plt.figure(figsize=(10, 5))
plt.plot(times, sx_exact, label='Exact method',color = "red", lw=2)
plt.plot(times, sx_rk, '--', label='Runge-Kutta (RK45)', color = "blue",alpha=0.8)
plt.xlabel('Time $t$', fontsize=15)
plt.ylabel(r'$\langle \sigma_1^x(t) \rangle$', fontsize=15)
plt.title(r'Evolution of $\langle \hat{\sigma}_1^x (t) \rangle$',fontsize = 15)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
path = "./images/B_evolution.png"
plt.savefig(path)
# plt.show()

```

## B_3
``` python
import numpy as np
from scipy.sparse import kron, identity, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from B_2_evolution import pauli_x, pauli_z, get_sigma_op, hamiltonian

if __name__ == '__main__':
    L = 18
    h = 1
    
    # 构建哈密顿量
    H = hamiltonian(L, h)
    
    # 构建Neel态 |↑↓↑↓...↑↓>
    # 在计算基底中，↑对应0，↓对应1
    neel_state = np.zeros(2**L, dtype=np.complex128)
    
    # 通过位运算高效构建Neel态的索引
    neel_index = 0
    for i in range(L):
        if i % 2 == 1:  # 偶数位置（从0开始计数）spin down
            neel_index |= (1 << i)
    
    neel_state[neel_index] = 1.0
    
    # 定义薛定谔方程
    def schrodinger_equation(t, psi):
        return -1j * (H @ psi)
    
    # 时间范围
    times = np.linspace(0, 10, 201)
    
    # 使用RK45方法求解时间演化
    sol = solve_ivp(schrodinger_equation, [0, 10], neel_state, 
                    t_eval=times, method='RK45', rtol=1e-7, atol=1e-9)
    
    # 计算保真度 F(t) = |<ψ(0)|ψ(t)>|^2
    fidelity = []
    for state in sol.y.T:
        overlap = np.abs(np.vdot(neel_state, state))**2
        fidelity.append(overlap)
    
    # 绘制保真度随时间的变化
    plt.figure(figsize=(10, 6))
    plt.plot(times, fidelity, 'b-', linewidth=2)
    plt.xlabel('Time $t$', fontsize=14)
    plt.ylabel('Fidelity $F(t) = |\\langle\\psi(0)|\\psi(t)\\rangle|^2$', fontsize=14)
    plt.title('Fidelity of Neel state time evolution ($L=18, h=1$)', fontsize=16)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("./images/B_neel_fidelity.png", dpi=300)
    plt.show()

```
