import math # balicek math potrebujeme kvuli cislu e
import numpy as np
import matplotlib.pyplot as plt

# definujeme nasi funkci F(x,y) = -3*y:
F = lambda x, y: -3 * y

# definujeme resic - Eulerovu metodu:
def EM(y_poc,T,N): 
    x = np.linspace(0, T, N + 1)
    y = np.zeros(N + 1)
    y_pres = np.zeros(N + 1)
    y[0] = y_poc;
    for k in range(N): # priblizne reseni     
        y[k + 1] = y[k] + h * F(x[k], y[k])
        y_pres[k] = np.exp(-3 * x[k]) # presne reseni pro porovnani  
    return y

N = 30 # pocet kroku; zkuste 10, 50, 100, ...
T = 5 # zajima nas y(T)
h = T / N # delka kroku
t = 0 # pocatecni cas
y_poc = 1 # pocatecni podminka y(0)
x = np.linspace(0, T, N + 1)
    
y = EM(y_poc,T,N) # volani Eulerovy metody

plt.plot(x, y, "b")
plt.ylabel('Vypoctene a presne reseni')
plt.legend(['vypoctene reseni', 'presne reseni'])
plt.show()
