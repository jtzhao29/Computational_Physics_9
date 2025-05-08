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