import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# PROYECTO FINAL: ECUACIONES DIFERENCIALES1 1
# Joel Jaquez - 23369 y Samuel Mejía - 23442
# Fecha: 4 de junio del 2025
# Alternativa escogida: Simulación a través de métodos numéricos
# ============================================================================

# ============================================================================
# MÉTODOS NUMÉRICOS GENERALES
# ============================================================================

def heun(f, y0, t0, tf, h):
    """Método de Heun para EDO de primer orden"""
    n = int((tf - t0) / h)
    t_values = np.linspace(t0, tf, n+1)
    y_values = np.zeros(n+1)
    y_values[0] = y0

    for i in range(n):
        t = t_values[i]
        y = y_values[i]
        k1 = f(t, y)
        k2 = f(t + h, y + h * k1)
        y_values[i+1] = y + (h / 2) * (k1 + k2)

    return t_values, y_values

def runge_kutta_4(f, y0, t0, tf, h):
    """Método RK4 para EDO de primer orden"""
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

def heun_sistema(f, u0, t0, tf, h):
    """Método de Heun para sistemas de EDOs"""
    n = int((tf - t0) / h)
    t_values = np.linspace(t0, tf, n+1)
    u_values = np.zeros((n+1, len(u0)))
    u_values[0] = u0

    for i in range(n):
        t = t_values[i]
        u = u_values[i]
        k1 = f(t, u)
        k2 = f(t + h, u + h * k1)
        u_values[i+1] = u + (h / 2) * (k1 + k2)

    return t_values, u_values

def runge_kutta_4_sistema(f, u0, t0, tf, h):
    """Método RK4 para sistemas de EDOs"""
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

# ============================================================================
# PROBLEMA 1: EDO DE PRIMER ORDEN
# dy/dt = -2y + 1, y(0) = 0
# Solución analítica: y(t) = 1/2(1 - e^(-2t))
# ============================================================================

def f1(t, y):
    return -2 * y + 1

def solucion_analitica_1(t):
    return 0.5 * (1 - np.exp(-2 * t))

# Parámetros
y0_1 = 0
t0_1 = 0
tf_1 = 1
h_1 = 0.05

# Métodos numéricos
t_heun_1, y_heun_1 = heun(f1, y0_1, t0_1, tf_1, h_1)
t_rk4_1, y_rk4_1 = runge_kutta_4(f1, y0_1, t0_1, tf_1, h_1)

# Solución analítica
t_analitica_1 = np.linspace(t0_1, tf_1, 1000)
y_analitica_1 = solucion_analitica_1(t_analitica_1)

# Gráfica 1
plt.figure(figsize=(10, 6))
plt.plot(t_heun_1, y_heun_1, 'ro-', markersize=6, linewidth=2, markeredgewidth=1, markeredgecolor='darkred', label='Método de Heun')
plt.plot(t_rk4_1, y_rk4_1, 'bs-', markersize=5, linewidth=2, markeredgewidth=1, markeredgecolor='darkblue', label='Método RK4')
plt.plot(t_analitica_1, y_analitica_1, 'k-', linewidth=3, label='Solución Analítica')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('EDO de Primer Orden: dy/dt = -2y + 1, y(0) = 0')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================================
# PROBLEMA 2: EDO DE SEGUNDO ORDEN NO HOMOGÉNEA
# y'' + 3y' + 2y = sin(t), y(0) = 1, y'(0) = 0
# Solución analítica: y(t) = (5/2)e^(-t) - (3/5)e^(-2t) - (3/10)cos(t) + (1/10)sin(t)
# ============================================================================

def sistema_segundo_orden(t, Y):
    """Convierte la EDO de segundo orden en sistema de primer orden"""
    y1, y2 = Y  # y1 = y, y2 = y'
    dy1dt = y2
    dy2dt = -3 * y2 - 2 * y1 + np.sin(t)
    return np.array([dy1dt, dy2dt])

def solucion_analitica_2(t):
    return (5/2) * np.exp(-t) - (3/5) * np.exp(-2*t) - (3/10) * np.cos(t) + (1/10) * np.sin(t)

# Parámetros
Y0_2 = np.array([1.0, 0.0])  # y(0) = 1, y'(0) = 0
t0_2 = 0
tf_2 = 1
h_2 = 0.05

# Métodos numéricos
t_heun_2, Y_heun_2 = heun_sistema(sistema_segundo_orden, Y0_2, t0_2, tf_2, h_2)
t_rk4_2, Y_rk4_2 = runge_kutta_4_sistema(sistema_segundo_orden, Y0_2, t0_2, tf_2, h_2)

# Solución analítica
t_analitica_2 = np.linspace(t0_2, tf_2, 1000)
y_analitica_2 = solucion_analitica_2(t_analitica_2)

# Gráfica 2
plt.figure(figsize=(10, 6))
plt.plot(t_heun_2, Y_heun_2[:, 0], 'ro-', markersize=6, linewidth=2, markeredgewidth=1, markeredgecolor='darkred', label='Método de Heun')
plt.plot(t_rk4_2, Y_rk4_2[:, 0], 'bs-', markersize=5, linewidth=2, markeredgewidth=1, markeredgecolor='darkblue', label='Método RK4')
plt.plot(t_analitica_2, y_analitica_2, 'k-', linewidth=3, label='Solución Analítica')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('EDO de Segundo Orden: y\'\' + 3y\' + 2y = sin(t), y(0) = 1, y\'(0) = 0')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================================
# PROBLEMA 3: SISTEMA 2x2 CON AUTOVALORES COMPLEJOS
# x' = 3x + 4y, y' = -4x + 3y, con x(0) = 1, y(0) = 0
# Solución analítica: x(t) = e^(3t)cos(4t), y(t) = -e^(3t)sin(4t)
# ============================================================================

def sistema_2x2(t, u):
    """Sistema 2x2"""
    x, y = u
    dxdt = 3 * x + 4 * y
    dydt = -4 * x + 3 * y
    return np.array([dxdt, dydt])

def solucion_analitica_3(t):
    """Solución analítica del sistema 2x2"""
    x_analitica = np.exp(3 * t) * np.cos(4 * t)
    y_analitica = -np.exp(3 * t) * np.sin(4 * t)
    return x_analitica, y_analitica

# Parámetros
u0_3 = np.array([1.0, 0.0])  # x(0) = 1, y(0) = 0
t0_3 = 0
tf_3 = 1  # Tiempo entre 0 y 1
h_3 = 0.01

# Métodos numéricos
t_heun_3, u_heun_3 = heun_sistema(sistema_2x2, u0_3, t0_3, tf_3, h_3)
t_rk4_3, u_rk4_3 = runge_kutta_4_sistema(sistema_2x2, u0_3, t0_3, tf_3, h_3)

# Solución analítica
t_analitica_3 = np.linspace(t0_3, tf_3, 1000)
x_analitica_3, y_analitica_3 = solucion_analitica_3(t_analitica_3)

# Gráfica 3a - Componente x(t)
plt.figure(figsize=(10, 6))
plt.plot(t_heun_3, u_heun_3[:, 0], 'ro-', markersize=6, linewidth=2, markeredgewidth=1, markeredgecolor='darkred', label='Método de Heun')
plt.plot(t_rk4_3, u_rk4_3[:, 0], 'bs-', markersize=5, linewidth=2, markeredgewidth=1, markeredgecolor='darkblue', label='Método RK4')
plt.plot(t_analitica_3, x_analitica_3, 'k-', linewidth=3, label='Solución Analítica')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Sistema 2x2: dx/dt = 3x + 4y, dy/dt = -4x + 3y - Componente x(t)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# Gráfica 3b - Componente y(t)
plt.figure(figsize=(10, 6))
plt.plot(t_heun_3, u_heun_3[:, 1], 'ro-', markersize=6, linewidth=2, markeredgewidth=1, markeredgecolor='darkred', label='Método de Heun')
plt.plot(t_rk4_3, u_rk4_3[:, 1], 'bs-', markersize=5, linewidth=2, markeredgewidth=1, markeredgecolor='darkblue', label='Método RK4')
plt.plot(t_analitica_3, y_analitica_3, 'k-', linewidth=3, label='Solución Analítica')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('Sistema 2x2: dx/dt = 3x + 4y, dy/dt = -4x + 3y - Componente y(t)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================================
# ANÁLISIS DETALLADO DE ERRORES
# ============================================================================

print("=" * 60)
print("           ANÁLISIS DETALLADO DE PRECISIÓN")
print("=" * 60)

print("\n1. EDO DE PRIMER ORDEN: dy/dt = -2y + 1, y(0) = 0")
print("-" * 50)
valor_analitico_1 = solucion_analitica_1(tf_1)
valor_heun_1 = y_heun_1[-1]
valor_rk4_1 = y_rk4_1[-1]
error_heun_1 = np.abs(valor_heun_1 - valor_analitico_1)
error_rk4_1 = np.abs(valor_rk4_1 - valor_analitico_1)

print(f"   Solución Analítica en t={tf_1}: {valor_analitico_1:.8f}")
print(f"   Método de Heun en t={tf_1}:     {valor_heun_1:.8f}")
print(f"   Método RK4 en t={tf_1}:        {valor_rk4_1:.8f}")
print(f"   Error Absoluto Heun:           {error_heun_1:.2e}")
print(f"   Error Absoluto RK4:            {error_rk4_1:.2e}")
print(f"   Mejora RK4 vs Heun:            {error_heun_1/error_rk4_1:.1f}x más preciso")

print("\n2. EDO DE SEGUNDO ORDEN: y'' + 3y' + 2y = sin(t), y(0) = 1, y'(0) = 0")
print("-" * 65)
valor_analitico_2 = solucion_analitica_2(tf_2)
valor_heun_2 = Y_heun_2[-1, 0]
valor_rk4_2 = Y_rk4_2[-1, 0]
error_heun_2 = np.abs(valor_heun_2 - valor_analitico_2)
error_rk4_2 = np.abs(valor_rk4_2 - valor_analitico_2)

print(f"   Solución Analítica en t={tf_2}: {valor_analitico_2:.8f}")
print(f"   Método de Heun en t={tf_2}:     {valor_heun_2:.8f}")
print(f"   Método RK4 en t={tf_2}:        {valor_rk4_2:.8f}")
print(f"   Error Absoluto Heun:           {error_heun_2:.2e}")
print(f"   Error Absoluto RK4:            {error_rk4_2:.2e}")
print(f"   Mejora RK4 vs Heun:            {error_heun_2/error_rk4_2:.1f}x más preciso")

print("\n3. SISTEMA 2x2: dx/dt = 3x + 4y, dy/dt = -4x + 3y, x(0) = 1, y(0) = 0")
print("-" * 68)
x_analitico_final, y_analitico_final = solucion_analitica_3(tf_3)
x_heun_final = u_heun_3[-1, 0]
y_heun_final = u_heun_3[-1, 1]
x_rk4_final = u_rk4_3[-1, 0]
y_rk4_final = u_rk4_3[-1, 1]

error_heun_3x = np.abs(x_heun_final - x_analitico_final)
error_rk4_3x = np.abs(x_rk4_final - x_analitico_final)
error_heun_3y = np.abs(y_heun_final - y_analitico_final)
error_rk4_3y = np.abs(y_rk4_final - y_analitico_final)

print(f"   COMPONENTE X(t):")
print(f"   Solución Analítica en t={tf_3}: {x_analitico_final:.8f}")
print(f"   Método de Heun en t={tf_3}:     {x_heun_final:.8f}")
print(f"   Método RK4 en t={tf_3}:        {x_rk4_final:.8f}")
print(f"   Error Absoluto Heun:           {error_heun_3x:.2e}")
print(f"   Error Absoluto RK4:            {error_rk4_3x:.2e}")
print(f"   Mejora RK4 vs Heun:            {error_heun_3x/error_rk4_3x:.1f}x más preciso")

print(f"\n   COMPONENTE Y(t):")
print(f"   Solución Analítica en t={tf_3}: {y_analitico_final:.8f}")
print(f"   Método de Heun en t={tf_3}:     {y_heun_final:.8f}")
print(f"   Método RK4 en t={tf_3}:        {y_rk4_final:.8f}")
print(f"   Error Absoluto Heun:           {error_heun_3y:.2e}")
print(f"   Error Absoluto RK4:            {error_rk4_3y:.2e}")
print(f"   Mejora RK4 vs Heun:            {error_heun_3y/error_rk4_3y:.1f}x más preciso")