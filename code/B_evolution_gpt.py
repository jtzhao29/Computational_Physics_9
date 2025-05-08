import numpy as np
from duijiao import Hamiltonian, sigmax1,get_sigma_op,pauli_x
from numba import njit, prange
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import kron, identity, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

L = 12
sigma1 = sigmax1(L)
# Modify the calculate_sigma_1_x function to use np.conj() and np.transpose()

def calculate_sigma_1_x(psi: np.ndarray) -> float:
    ex_sigma_1_x = np.real(np.conj(psi).T @ sigma1 @ psi ) # Use np.conj() and np.transpose() instead
    return ex_sigma_1_x


# Main evolution function (ensure both psi and H are complex arrays)

def evolution_exact(psi0: np.ndarray, H: np.ndarray, t: float) -> np.ndarray:
    expHt = np.exp(-1j * H * t)  # exp(-iHt) for time evolution
    psi_t = expHt @ psi0
    ex_sigma_1_x = calculate_sigma_1_x(psi_t)
    return psi_t, ex_sigma_1_x

# Runge-Kutta 方法
# @njit(parallel=True)
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

# 绘制时间演化结果
def plot_evolution(sigma1_x_t_exact, sigma1_x_t_runge, t):
    plt.figure(figsize=(10, 6))
    plt.plot(t, sigma1_x_t_exact, label="Exact Evolution")
    plt.plot(t, sigma1_x_t_runge, label="Runge-Kutta Evolution", linestyle='dashed')
    plt.xlabel('Time (t)')
    plt.ylabel(r'$\langle \hat{\sigma}_1^x (t) \rangle$')
    plt.legend()
    plt.show()

def schrodinger_equation(t, psi):
    return -1j * (H1 @ psi)


# 主函数
if __name__ == "__main__":
    # 设置参数
    h0 = 0.5
    h1 = 3.0
    J = 1
    L=12
    H0 = Hamiltonian(L, J, h0)
    H1 = Hamiltonian(L, J, h1)
    sx1 = get_sigma_op(pauli_x(), 0, L)  

    # 计算基态
    eigvals, eigvecs =eigsh(H0, k=1, which='SA')
    ground_state = eigvecs[:, 0].astype(np.complex128)
    E1, V1 = np.linalg.eigh(H1.toarray())  
    c = V1.conj().T @ ground_state  

    # 设置时间范围
    t0 = 0
    tf = 20
    dt = 1
    t = np.arange(t0, tf, dt)

    sx_exact = []
    for t in t:
        # psi_t = psi0@np.exp(-1j * E1 * t)
        psi_t = V1 @ (c * np.exp(-1j * E1 * t))  
        sx_exact.append(np.real(psi_t.conj().T @ (sx1 @ psi_t))) 

    sol = solve_ivp(schrodinger_equation, [0, 20], ground_state, t_eval=t, 
                    method='RK45',rtol=1e-6,atol=1e-8)
    sx_rk = [np.real(psi.conj().T @ (sx1 @ psi)) for psi in sol.y.T]

    # # 精确演化
    # sigma1_x_t_exact = np.zeros(len(t))
    # for i in range(len(t)):
    #     # print(t[i])
    #     _, sigma1_x_t_exact[i] = evolution_exact(ground_state, H1, t[i])

    # # Runge-Kutta演化
    # _, sigma1_x_t_runge = runge_kutta(ground_state, H1, t)

    # 保存数据
    # pd_sigma1_x_t_exact = pd.DataFrame(sigma1_x_t_exact)
    # pd_sigma1_x_t_runge = pd.DataFrame(sigma1_x_t_runge)
    # pd_sigma1_x_t_exact.to_csv("./data/sigma1_x_t_exact.csv")
    # pd_sigma1_x_t_runge.to_csv("./data/sigma1_x_t_runge.csv")

    # # 绘制结果
    # plot_evolution(sigma1_x_t_exact, sigma1_x_t_runge, t)

    plot_evolution(sx_exact, sx_rk, t)