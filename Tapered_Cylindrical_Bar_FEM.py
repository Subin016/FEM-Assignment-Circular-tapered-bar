import numpy as np
import matplotlib.pyplot as plt

L = 1000      # mm
d1 = 50       # mm
d2 = 25       # mm
E = 200000    # N/mm²
P = 10000     # N

# User input
n = int(input("Enter the number of elements (n): "))
# Length of an element
Le = L/n
#No. of Nodes
nn = n+1
#Creates Global Stiffness Matrix with Zero Value.
K = np.zeros((nn,nn))

for e in range(n):

    xmid = (e+0.5)*Le

    dmid = d1 + (d2-d1)*xmid/L

    A = np.pi*dmid**2/4

    ke = (E*A/Le)*np.array([[1,-1],
                            [-1,1]])

    K[e:e+2,e:e+2] += ke

F = np.zeros(nn)
F[-1] = P

Kr = K[1:,1:]
Fr = F[1:]

ur = np.linalg.solve(Kr,Fr)

u = np.zeros(nn)
u[1:] = ur

print("\nFree-end displacement =", u[-1])

delta_exact = 4*P*L/(np.pi*E*d1*d2)

print("Analytical displacement =", delta_exact)

error = abs(delta_exact-u[-1])/delta_exact*100

print("Percentage error =", error,"%")


# -----------------------------
# Plot FEM Mesh
# -----------------------------

x = np.linspace(0, L, 200)

# Actual tapered profile
d = d1 + (d2 - d1) * x / L

plt.figure(figsize=(10,4))

# Upper and lower boundaries
plt.plot(x, d/2)
plt.plot(x, -d/2)

# Element boundaries
for i in range(n+1):
    xi = i * Le
    di = d1 + (d2-d1) * xi / L

    plt.plot([xi, xi],
             [-di/2, di/2],
             'k--',
             linewidth=0.8)

# Node labels
for i in range(n+1):
    xi = i * Le
    plt.text(xi, 0, f'N{i+1}',
             ha='center',
             va='center')

plt.title(f'FEM Discretization of Tapered Bar (n = {n} elements)')
plt.xlabel('Length')
plt.ylabel('Diameter')
plt.axis('equal')
plt.grid(True)

plt.show()
