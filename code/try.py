import numpy as np
from scipy.sparse import kron, identity

# 定义 Pauli 矩阵和单位矩阵
sigma_z = np.array([[1, 0], [0, -1]])
I = np.identity(2)

def H_12():
    H = kron(kron(sigma_z, sigma_z), I).toarray()  # 转换为稠密矩阵
    return H

def H_13():
    H = kron(kron(sigma_z, I), sigma_z).toarray()  # 转换为稠密矩阵
    return H        



if __name__ == "__main__":
    # 打印构造的哈密顿量矩阵
    print("12",H_12())
    print("13", H_13())