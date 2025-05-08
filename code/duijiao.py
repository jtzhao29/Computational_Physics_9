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


# for test
# if __name__ == "__main__":
#     H = Hamiltonian(3,1,2)
#     print(H)
# import numpy as np
# from scipy.linalg import eigh
# from numba import njit
# import matplotlib.pyplot as plt
# from joblib import Parallel, delayed

# @njit
# def Hamiltonian(L, J, h):
#     H = np.zeros((2**L, 2**L), dtype=np.complex128)  # 更高效的初始化
#     I = np.array([[1, 0], [0, 1]], dtype=np.complex128)
#     sigmaz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
#     sigmax = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    
#     for i in range(L):
#         h2 = sigmax if i == 0 else I
#         h1 = sigmaz if i == 0 and i == L - 1 else I
        
#         for j in range(1, L):
#             if j == i and j == i + 1:
#                 h1 = np.kron(h1, sigmaz)
#             else:
#                 h1 = np.kron(h1, I)
#             if j == i:
#                 h2 = np.kron(h2, sigmax)
#             else:
#                 h2 = np.kron(h2, I)
#         H += -J * h1 - h * h2
#     return H

# @njit
# def sigmax1(L):
#     I = np.array([[1, 0], [0, 1]], dtype=np.complex128)
#     sigmax = np.array([[0, 1], [1, 0]], dtype=np.complex128)
#     a = sigmax
#     for i in range(L - 1):
#         a = np.kron(a, I)
#     return a

# @njit
# def expm(H, eigvals, eigvecs, t):
#     di = np.exp(1j * eigvals * t)
#     Lam = np.diag(di)
#     return eigvecs @ Lam @ eigvecs.conj().T

# # 用并行处理不同 h 的结果
# def compute_for_h(h, L, J):
#     H = Hamiltonian(L, J, h)
#     eigvals, eigvecs = eigh(H)
#     print(f"when h = {h}, the ground state energy is:{eigvals[0]}, and the first active state energy is:{eigvals[1]}")
    
#     ground_state = eigvecs[:, 0]
    
#     # 计算第一个格点的横场方向磁化强度的期望值
#     sigma1 = sigmax1(L)
#     ex_sigma_1_x = np.abs(ground_state.conj().T @ sigma1 @ ground_state)[0]
#     print(f"when h = {h}, the expectation value of sigma1_x is: {ex_sigma_1_x}")
    
#     return eigvals[0], eigvals[1], ex_sigma_1_x

# if __name__ == "__main__":
#     L = 12
#     J = 1
#     h_list = [0.5, 1.0, 2.0]  # 可以增加多个 h 值进行并行计算
    
#     # 并行计算不同 h 值的结果
#     results = Parallel(n_jobs=-1)(delayed(compute_for_h)(h, L, J) for h in h_list)
    
#     # 打印所有结果
#     for h, result in zip(h_list, results):
#         print(f"Results for h={h}: {result}")
