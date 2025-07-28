import numpy as np
N = 64
np.random.seed(0)
# random diagonals with non-zero main diag
b = np.random.randn(N) + 1j*np.random.randn(N) + 5
c = np.random.randn(N) + 1j*np.random.randn(N)
a = np.random.randn(N) + 1j*np.random.randn(N)
# ensure first a[0], last c[N-1] not used
c[-1] = 0
a[0] = 0
# RHS
d = np.random.randn(N) + 1j*np.random.randn(N)

# solve using numpy's linear solver
A = np.zeros((N,N),dtype=np.complex64)
for i in range(N):
    A[i,i] = b[i]
    if i>0:
        A[i,i-1] = a[i]
    if i<N-1:
        A[i,i+1] = c[i]

x = np.linalg.solve(A,d)

def save(fname, arr):
    with open(fname,'w') as f:
        for v in arr:
            f.write(f"{v.real} {v.imag}\n")

save('test/a.dat', a)
save('test/b.dat', b)
save('test/c.dat', c)
save('test/d.dat', d)
save('test/golden.dat', x)
print('data generated')
