import torch
import numpy as np
import matplotlib.pyplot as plt
from torch.linalg import eig

# PyTorch 版本的 Hamiltonian 构造
def Hamiltonian(L, J, h, device):
    I = torch.eye(2, dtype=torch.complex128, device=device)
    sigmaz = torch.tensor([[1, 0], [0, -1]], dtype=torch.complex128, device=device)
    sigmax = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128, device=device)
    
    H = torch.zeros((2**L, 2**L), dtype=torch.complex128, device=device)

    for i in range(L):
        h2 = sigmax if i == 0 else I
        h1 = sigmaz if i == 0 and i == L - 1 else I
        
        for j in range(1, L):
            if j == i and j == i + 1:
                h1 = torch.kron(h1, sigmaz)
            else:
                h1 = torch.kron(h1, I)
            if j == i:
                h2 = torch.kron(h2, sigmax)
            else:
                h2 = torch.kron(h2, I)
        H += -J * h1 - h * h2
    return H

# 创建 sigmax 张量
def sigmax1(L, device):
    I = torch.eye(2, dtype=torch.complex128, device=device)
    sigmax = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128, device=device)
    
    a = sigmax
    for i in range(L - 1):
        a = torch.kron(a, I)
    return a

# 计算时间演化
def time_evolution(H, eigvals, eigvecs, t, device):
    di = torch.exp(1j * eigvals * t).to(device)
    Lam = torch.diag(di)
    result = torch.matmul(torch.matmul(eigvecs, Lam), eigvecs.conj().T)
    return result

def compute_for_h(h, L, J, device):
    # 构造哈密顿量
    H = Hamiltonian(L, J, h, device)
    
    # 计算特征值和特征向量
    eigvals, eigvecs = eig(H)
    eigvals = eigvals.real  # 只取实部
    
    print(f"when h = {h}, the ground state energy is:{eigvals[0]}, and the first active state energy is:{eigvals[1]}")
    
    ground_state = eigvecs[:, 0]
    
    # 计算第一个格点的横场方向磁化强度的期望值
    sigma1 = sigmax1(L, device)
    ex_sigma_1_x = torch.abs(torch.matmul(ground_state.conj().T, torch.matmul(sigma1, ground_state)))[0]
    
    print(f"when h = {h}, the expectation value of sigma1_x is: {ex_sigma_1_x}")
    
    return eigvals[0], eigvals[1], ex_sigma_1_x

if __name__ == "__main__":
    L = 12
    J = 1
    h_list = [0.5, 1.0, 2.0]  # 可以增加多个 h 值进行并行计算
    
    # 将数据移到 GPU 上
    device = torch.device('cuda:0')  # 选择 A100 GPU
    
    # 并行计算不同 h 值的结果
    results = []
    for h in h_list:
        results.append(compute_for_h(h, L, J, device))
    
    # 打印所有结果
    for h, result in zip(h_list, results):
        print(f"Results for h={h}: {result}")
