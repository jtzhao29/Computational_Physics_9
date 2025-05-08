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
    # plt.figure(figsize=(10, 6))
    # plt.plot(times, fidelity, 'b-', linewidth=2)
    # plt.xlabel('Time $t$', fontsize=14)
    # plt.ylabel('Fidelity $F(t) = |\\langle\\psi(0)|\\psi(t)\\rangle|^2$', fontsize=14)
    # plt.title('Fidelity of Neel state time evolution ($L=18, h=1$)', fontsize=16)
    # plt.grid(alpha=0.3)
    # plt.tight_layout()
    # plt.savefig("./images/B_neel_fidelity.png", dpi=300)
    # plt.show()


    plt.figure(figsize=(10, 6))
    plt.plot(times, fidelity, 'b-',color = "b", linewidth=2)
    plt.xlabel('Time $t$')
    plt.ylabel('Fidelity ')
    plt.title('Fidelity ')
    # plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("./images/B_neel_fidelity_2.png", dpi=300)
    plt.show()


