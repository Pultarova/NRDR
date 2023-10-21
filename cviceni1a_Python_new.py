import math # balicek math potrebujeme kvuli cislu e
import numpy as np
import matplotlib.pyplot as plt

# definujeme nasi funkci F(x,y) = -3*y
F = lambda x, y: -3 * y

N = 30 # pocet kroku; zkuste 10, 50, 100, ...
T = 5 # zajima nas y(T)
h = T / N # delka kroku
t = 0 # pocatecni cas
y_poc = 1 # pocatecni podminka y(0)

# vypocet jen y(T):
y = y_poc
for k in range(N):
    y1 = y + h * F(t, y)
    y = y1

print("Vypoctene reseni:", y)
print("Presne reseni:", np.exp(-3 * T))

# vypocet a vykresleni reseni na celem intervalu (0,T):
t = 0
x = np.linspace(0, T, N + 1)
y = np.zeros(N + 1)
y_pres = np.zeros(N + 1)
y[0] = y_poc;

for k in range(N): # priblizne reseni     
    y[k + 1] = y[k] + h * F(x[k], y[k])
    y_pres[k] = np.exp(-3 * x[k]) # presne reseni pro porovnani

plt.plot(x, y, "b")
plt.plot(x, y_pres, "r--")
plt.ylabel('Vypoctene a presne reseni')
plt.legend(['vypoctene reseni', 'presne reseni'])
plt.show()
