import numpy as np
import matplotlib.pyplot as plt

# Definir la función 
def f(t, y):
    return -2 * y + 1

# Método RK4
def runge_kutta_4(f, y0, t0, tf, h):
    n = int((tf - t0) / h)
    t_values = np.linspace(t0, tf, n+1)
    y_values = np.zeros(n+1)
    y_values[0] = y0

    for i in range(n):
        t = t_values[i]
        y = y_values[i]
        
        k1 = h * f(t, y)
        k2 = h * f(t + h/2, y + k1/2)
        k3 = h * f(t + h/2, y + k2/2)
        k4 = h * f(t + h, y + k3)
        
        y_values[i+1] = y + (k1 + 2*k2 + 2*k3 + k4)/6

    return t_values, y_values

# Parámetros
y0 = 0
t0 = 0
tf = 5
h = 0.1

# RK4
t_vals, y_vals = runge_kutta_4(f, y0, t0, tf, h)

# Visualizar
plt.plot(t_vals, y_vals, label="RK4")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Solución de dy/dt = -2y + 1 usando RK4")
plt.grid(True)
plt.legend()
plt.show()