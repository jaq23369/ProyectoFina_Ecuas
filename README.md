# Proyecto Final - Ecuaciones Diferenciales 🧮

Este proyecto implementa métodos numéricos para resolver ecuaciones diferenciales ordinarias (EDOs) de primer y segundo orden, así como sistemas de EDOs 2x2. El objetivo es simular soluciones usando los métodos de **Heun** y **Runge-Kutta de 4to orden (RK4)** en Python.

---

## 📌 Ecuaciones resueltas

### 1. EDO de primer orden

```math
\frac{dy}{dt} = -2y + 1,\quad y(0) = 0
```

### 2. Sistema de EDOs 2x2

```math
\begin{cases}
\frac{dx}{dt} = 0.3x + 0.4y \\
\frac{dy}{dt} = -0.4x + 0.3y
\end{cases},\quad x(0) = 1,\ y(0) = 0
```

### 3. EDO de segundo orden no homogénea

```math
\frac{d^2y}{dt^2} + 3\frac{dy}{dt} + 2y = sen(t),\quad y(0) = 1,\ y'(0) = 0
```

---

## 🧠 Métodos numéricos utilizados

- Método de Heun (Euler mejorado)
- Método de Runge-Kutta de cuarto orden (RK4)

Ambos métodos fueron aplicados a cada ecuación para comparar la precisión y comportamiento de las soluciones.

---

## 📦 Requisitos

Este proyecto requiere Python 3.8+ y utiliza las siguientes bibliotecas:

- NumPy
- Matplotlib

Estas librerías ya vienen preinstaladas en algunos entornos como Anaconda o Google Colab.  
Si estás usando una instalación limpia de Python, seguí estas instrucciones según tu sistema operativo:

### 🔧 Instalación en Windows

1. Hay que tener Python y `pip` instalados.
2. Abrí CMD o PowerShell y ejecuta:

```bash
pip install numpy matplotlib
```

### 🔧 Instalación en macOS

1. Asegurate de tener Python 3 instalado (puedes usar Homebrew).
2. Abrí la Terminal y ejecuta:

```bash
brew install python
pip3 install numpy matplotlib
```

---

## ▶️ Cómo ejecutar

Desde la terminal:

```bash
python3 main.py
```

Se generarán varias gráficas mostrando las soluciones aproximadas con los métodos de Heun y RK4.

---

## 📈 Salidas esperadas

- Gráficas comparativas de las soluciones para cada ecuación.
- Visualización del comportamiento dinámico de cada sistema.

---

## 👨‍💻 Autores

- Joel Jaquez - 23369  
- Samuel Mejía - 23442 
