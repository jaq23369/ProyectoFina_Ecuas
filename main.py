import numpy as np
import matplotlib.pyplot as plt

# Función 1
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


# 2x2
def sistema(t, u):
    x, y = u
    dxdt = 3 * x + 4 * y
    dydt = -4 * x + 3 * y
    return np.array([dxdt, dydt])

def runge_kutta_4_sistema(f, u0, t0, tf, h):
    n = int((tf - t0) / h)
    t_values = np.linspace(t0, tf, n+1)
    u_values = np.zeros((n+1, len(u0)))
    u_values[0] = u0

    for i in range(n):
        t = t_values[i]
        u = u_values[i]
        k1 = h * f(t, u)
        k2 = h * f(t + h/2, u + k1/2)
        k3 = h * f(t + h/2, u + k2/2)
        k4 = h * f(t + h, u + k3)
        u_values[i+1] = u + (k1 + 2*k2 + 2*k3 + k4)/6

    return t_values, u_values

# Parámetros iniciales
u0 = np.array([1.0, 0.0])  # x(0) = 1, y(0) = 0
t0_sistema = 0
tf_sistema = 5
h_sistema = 0.1

# RK4
t_sys, u_sys = runge_kutta_4_sistema(sistema, u0, t0_sistema, tf_sistema, h_sistema)

# Graficar resultados
plt.plot(t_sys, u_sys[:, 0], 'r-', label="x(t) RK4")
plt.plot(t_sys, u_sys[:, 1], 'b-', label="y(t) RK4")
plt.xlabel("t")
plt.ylabel("Valores")
plt.title("Sistema dx/dt = 3x + 4y, dy/dt = -4x + 3y usando RK4")
plt.grid(True)
plt.legend()
plt.show()

