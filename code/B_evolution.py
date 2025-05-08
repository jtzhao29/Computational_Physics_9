import numpy as np
from duijiao import Hamiltonian,sigmax1,expm
import pandas as pd
import matplotlib.pyplot as plt
from numba import njit 
import numpy as np
from duijiao import Hamiltonian, sigmax1
from numba import njit, prange
import pandas as pd
import matplotlib.pyplot as plt
# 2. $t=0$ 时系统处于 $h=0.5$ 的基态，$t=0^{+}$ 的瞬间哈密顿量参数变为 $h=3.0$，求在此哈密顿量下的时间演化，计算 $\langle \hat{\sigma}_{1}^{x}(t) \rangle$ 在时间 $t \in [0, 20]$ 的变化情况，并画出来。你需要给出精确的结果和Runge-Kutta方法的结果。验证两种方法的结果是一样的。（3分）
global sigma1
L=12
sigma1 = sigmax1(L)

# def calculate_sigma_1_x(psi:np.ndarray)->float:
    
#     ex_sigma_1_x = psi.conj().T @ sigma1 @ psi
#     return ex_sigma_1_x

# 优化计算
@njit
def calculate_sigma_1_x(psi: np.ndarray) -> float:
    ex_sigma_1_x = np.vdot(psi, sigma1 @ psi)
    return ex_sigma_1_x


@njit
# def evolution_exact(psi0:np.ndarray,H:np.ndarray,t:float)->tuple[np.ndarray,float]:
#     """
#     计算在哈密顿量H下，初始态psi0的时间演化
#     :param psi0: 初始态
#     :param H: 哈密顿量
#     :param t: 时间
#     :return: 演化后的态
#     """
    
#     # $$\psi(t) = e^{-iHt} \psi(0)$$
#     psi_t = np.exp(-1j * H * t) @ psi0
#     ex_sigma_1_x = calculate_sigma_1_x
#     return psi_t,ex_sigma_1_x

def evolution_exact(psi0: np.ndarray, H: np.ndarray, t: float) -> np.ndarray:
    """
    精确时间演化计算
    """
    # 使用矩阵指数进行演化
    expHt = np.exp(-1j * H * t)
    psi_t = expHt @ psi0
    ex_sigma_1_x = calculate_sigma_1_x(psi_t)
    return psi_t, ex_sigma_1_x



# @ njit
# def runge_kutta(psi0:np.ndarray,H:np.ndarray,t:np.ndarray)->tuple[np.ndarray,np.ndarray]:
#     """
#     使用Runge-Kutta方法计算时间演化
#     :param psi0: 初始态
#     :param H: 哈密顿量
#     :param t0: 初始时间
#     :param tf: 结束时间
#     :param t: 时间点
#     :return: 演化后的态
#     """
#     # $\frac{d\psi}{dt} = -iH \psi$$
    
#     dt = (tf - t0) / len(t)
#     psi = np.zeros((len(t), len(psi0)), dtype=np.complex128)
#     psi[0] = psi0

#     sigma1_x_t = np.zeros(len(t))
#     sigma1_x_t[0] = calculate_sigma_1_x(psi0)
#     for i in range(1, len(t)):
#         print(i)
#         k1 = -1j * H @ psi[i - 1] * dt
#         k2 = -1j * H @ (psi[i - 1] + k1 / 2) * dt
#         k3 = -1j * H @ (psi[i - 1] + k2 / 2) * dt
#         k4 = -1j * H @ (psi[i - 1] + k3) * dt
#         psi[i] = psi[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
#         sigma1_x_t[i] = calculate_sigma_1_x(psi[i])
#     return psi,sigma1_x_t

#优化后
@njit(parallel=True)
def runge_kutta(psi0: np.ndarray, H: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    使用Runge-Kutta方法进行时间演化
    """
    dt = t[1] - t[0]  # 计算时间步长
    psi = np.zeros((len(t), len(psi0)), dtype=np.complex128)
    sigma1_x_t = np.zeros(len(t))

    psi[0] = psi0
    sigma1_x_t[0] = calculate_sigma_1_x(psi0)
    
    for i in prange(1, len(t)):
        k1 = -1j * H @ psi[i - 1] * dt
        k2 = -1j * H @ (psi[i - 1] + k1 / 2) * dt
        k3 = -1j * H @ (psi[i - 1] + k2 / 2) * dt
        k4 = -1j * H @ (psi[i - 1] + k3) * dt
        psi[i] = psi[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        sigma1_x_t[i] = calculate_sigma_1_x(psi[i])

    return psi, sigma1_x_t

def plot_evolution(sigma1_x_t_exact, sigma1_x_t_runge, t):
    plt.figure(figsize=(10, 6))
    plt.plot(t, sigma1_x_t_exact, label="Exact Evolution")
    plt.plot(t, sigma1_x_t_runge, label="Runge-Kutta Evolution", linestyle='dashed')
    plt.xlabel('Time (t)',fontsize = 15)
    plt.ylabel(r'$\langle \hat{\sigma}_1^x (t) \rangle$',fontsize = 15)
    plt.title('Evolution of $\langle \hat{\sigma}_1^x (t) \rangle$',fontsize = 15)
    plt.legend()
    plt.grid(True)
    path = r"./images/B_evolution.png"
    plt.show()

# if __name__ == "__main__":
#     L = 12
#     J = 1
#     h0 = 0.5
#     h1 = 3.0
#     sigma1 = sigmax1(L)
#     H0 = Hamiltonian(L, J, h0)
#     H1 = Hamiltonian(L, J, h1)
#     eigvals, eigvecs = np.linalg.eigh(H0)
#     # ground_state = eigvecs[:, 0]
#     ground_state = eigvecs[:, 0].astype(np.complex128)
#     # 计算随时间的演化，在H1下的演化
#     t0=0
#     tf = 20
#     dt = 1
#     # 写出t
#     t = np.arange(t0, tf, dt)
#     psi_t_exact = np.zeros((len(t), len(ground_state)), dtype=complex)
#     sigma1_x_t_exact = np.zeros((len(t)))
#     for i in range(len(t)):
#         print(i)
#         psi_t_exact[i],sigma1_x_t_exact[i]= evolution_exact(ground_state,H1,t[i])

#     psi_t_runge,sigma1_x_t_runge = runge_kutta(ground_state, H1, t)

    # pd_psi_t_exact = pd.DataFrame(psi_t_exact)
    # pd_psi_t_runge = pd.DataFrame(psi_t_runge)

    # pd_psi_t_exact.to_csv("./data/psi_t_exact.csv")
    # pd_psi_t_runge.to_csv("./data/psi_t_runge.csv")

    # pd_sigma1_x_t_exact = pd.DataFrame(sigma1_x_t_exact)
    # pd_sigma1_x_t_runge = pd.DataFrame(sigma1_x_t_runge)
    # pd_sigma1_x_t_exact.to_csv("./data/sigma1_x_t_exact.csv")
    # pd_sigma1_x_t_runge.to_csv("./data/sigma1_x_t_runge.csv")


if __name__ == "__main__":
    # 设置参数
    h0 = 0.5
    h1 = 3.0
    J = 1
    H0 = Hamiltonian(L, J, h0)
    H1 = Hamiltonian(L, J, h1)

    # 计算基态
    eigvals, eigvecs = np.linalg.eigh(H0)
    ground_state = eigvecs[:, 0].astype(np.complex128)

    # 设置时间范围
    t0 = 0
    tf = 20
    dt = 1
    t = np.arange(t0, tf, dt)

    # 精确演化
    sigma1_x_t_exact = np.zeros(len(t))
    for i in range(len(t)):
        _, sigma1_x_t_exact[i] = evolution_exact(ground_state, H1, t[i])

    # Runge-Kutta演化
    _, sigma1_x_t_runge = runge_kutta(ground_state, H1, t)

    # 保存数据
    pd_sigma1_x_t_exact = pd.DataFrame(sigma1_x_t_exact)
    pd_sigma1_x_t_runge = pd.DataFrame(sigma1_x_t_runge)
    pd_sigma1_x_t_exact.to_csv("./data/sigma1_x_t_exact.csv")
    pd_sigma1_x_t_runge.to_csv("./data/sigma1_x_t_runge.csv")

    # 绘制结果
    plot_evolution(sigma1_x_t_exact, sigma1_x_t_runge, t)




