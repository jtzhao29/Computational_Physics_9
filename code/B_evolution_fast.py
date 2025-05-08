import numpy as np
from duijiao import Hamiltonian,sigmax1,expm
import pandas as pd
import matplotlib.pyplot as plt
from numba import njit as jit
import scipy.sparse as sp
from scipy.sparse.linalg import expm

# 2. $t=0$ 时系统处于 $h=0.5$ 的基态，$t=0^{+}$ 的瞬间哈密顿量参数变为 $h=3.0$，求在此哈密顿量下的时间演化，计算 $\langle \hat{\sigma}_{1}^{x}(t) \rangle$ 在时间 $t \in [0, 20]$ 的变化情况，并画出来。你需要给出精确的结果和Runge-Kutta方法的结果。验证两种方法的结果是一样的。（3分）
@jit
def evolution_exact(psi0: np.ndarray, H: sp.csr_matrix, t: float) -> np.ndarray:
    """
    使用稀疏矩阵计算在哈密顿量 H 下，初始态 psi0 的时间演化
    :param psi0: 初始态
    :param H: 稀疏哈密顿量矩阵
    :param t: 时间
    :return: 演化后的态
    """
    # 计算矩阵指数 e^(-iHt) 直接使用稀疏矩阵
    expm_Ht = expm(-1j * H * t)
    psi_t = expm_Ht @ psi0  # 演化后的态
    return psi_t

@ jit
def runge_kutta(psi0:np.ndarray,H:np.ndarray,t:np.ndarray)->np.ndarray:
    """
    使用Runge-Kutta方法计算时间演化
    :param psi0: 初始态
    :param H: 哈密顿量
    :param t0: 初始时间
    :param tf: 结束时间
    :param t: 时间点
    :return: 演化后的态
    """
    # $\frac{d\psi}{dt} = -iH \psi$$
    
    dt = (tf - t0) / len(t)
    psi = np.zeros((len(t), len(psi0)), dtype=np.complex128)
    psi[0] = psi0

    for i in range(1, len(t)):
        print(i)
        k1 = -1j * H @ psi[i - 1] * dt
        k2 = -1j * H @ (psi[i - 1] + k1 / 2) * dt
        k3 = -1j * H @ (psi[i - 1] + k2 / 2) * dt
        k4 = -1j * H @ (psi[i - 1] + k3) * dt
        psi[i] = psi[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return psi

def plot_evolution(psi:np.ndarray, t:np.ndarray, L=12, J=1, h0=0, h1=3):
    plt.figure(figsize=(10, 6))
    plt.plot(t, np.abs(psi[:, 0])**2, label=f"h = {h0}")




if __name__ == "__main__":
    L = 12
    J = 1
    h0 = 0.5
    h1 = 3.0
    sigma1 = sigmax1(L)
    H0 = Hamiltonian(L, J, h0)
    H1 = Hamiltonian(L, J, h1)
    eigvals, eigvecs = np.linalg.eigh(H0)
    # ground_state = eigvecs[:, 0]
    ground_state = eigvecs[:, 0].astype(np.complex128)
    # 计算随时间的演化，在H1下的演化
    t0=0
    tf = 20
    dt = 1
    # 写出t
    t = np.arange(t0, tf, dt)
    psi_t_exact = np.zeros((len(t), len(ground_state)), dtype=complex)
    for i in range(len(t)):
        print(i)
        psi_t_exact[i] = evolution_exact(ground_state,H1,t[i])

    psi_t_runge = runge_kutta(ground_state, H1, t)

    pd_psi_t_exact = pd.DataFrame(psi_t_exact)
    pd_psi_t_runge = pd.DataFrame(psi_t_runge)

    pd_psi_t_exact.to_csv("./data/psi_t_exact.csv")
    pd_psi_t_runge.to_csv("./data/psi_t_runge.csv")



