# reseni DR y'=-3*y, y(0)=1 na intervalu (0,T)
# (presne reseni je y(x)=e^(-3*x))
import math # balicek math potrebujeme kvuli cisu e
N = 100 # pocet kroku; zkuste 5, 10, 50
T = 5 # zajima nas y(T)
h = T/N # delka kroku
t = 0 # pocatecni cas
y0 = 1 # pocatecni podminka y(0)
for k in range(N):    
    y1 = y0 + h * (-3)* y0
    y0 = y1 
print("vypoctene reseni",y1)
print("presne reseni",math.exp(-3*T))

# vykresleni reseni na celem intervalu (0,T):
import array as arr
#from array import *
t = 0
x = arr.array('f',[0])
y = arr.array('f',[y0])
for k in range(N):    
    y1 = y0 + h * (-3)* y0
    y0 = y1 
    y.append(y1)
    x.append(k*h)
import matplotlib.pyplot as plt
plt.plot(x,y)
plt.ylabel('vypoctene reseni')
plt.show()    