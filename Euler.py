import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Entry, Button, StringVar, messagebox, Frame

# Método de Euler para una EDO
def euler_method(f, y0, t0, t_end, h):
    t_values = np.arange(t0, t_end + h, h)  # Puntos en el tiempo
    y_values = [y0]  # Solución inicial
    
    for t in t_values[:-1]:
        y_next = y_values[-1] + h * f(t, y_values[-1])
        y_values.append(y_next)
    
    return t_values, np.array(y_values)

# EDO de primer orden: dy/dt = -k * y
def edo_first_order(t, y, k=1):
    return -k * y

# Sistema Lotka-Volterra
def lotka_volterra(t, state, alpha=0.1, beta=0.02, delta=0.01, gamma=0.1):
    x, y = state
    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y
    return np.array([dxdt, dydt])

# Método de Euler para sistemas
def euler_system(f, state0, t0, t_end, h):
    t_values = np.arange(t0, t_end + h, h)
    state_values = [state0]
    
    for t in t_values[:-1]:
        state_next = state_values[-1] + h * f(t, state_values[-1])
        state_values.append(state_next)
    
    return t_values, np.array(state_values)

# Función para mostrar pantalla de EDO
def show_edo_screen():
    main_frame.pack_forget()
    edo_frame.pack()

# Función para mostrar pantalla de Lotka-Volterra
def show_lotka_volterra_screen():
    main_frame.pack_forget()
    lv_frame.pack()

# Función para resolver la EDO
def solve_edo():
    try:
        k = float(k_var.get())
        y0_edo = float(y0_var.get())

        if not (0.1 <= k <= 5):
            raise ValueError("k fuera del rango [0.1, 5]")
        if not (0 < y0_edo <= 100):
            raise ValueError("y₀ fuera del rango (0, 100]")

        t_values, y_values = euler_method(lambda t, y: edo_first_order(t, y, k), y0_edo, 0, 10, 0.1)

        # Solución analítica
        y_analytical = y0_edo * np.exp(-k * t_values)

        # Graficar resultados EDO de primer orden
        plt.figure(figsize=(10, 6))
        plt.plot(t_values, y_values, label="Euler (Numérico)", marker='o')
        plt.plot(t_values, y_analytical, label="Analítico", linestyle='--')
        plt.title("EDO de Primer Orden: $\\frac{dy}{dt} = -ky$")
        plt.xlabel("Tiempo (t)")
        plt.ylabel("y(t)")
        plt.legend()
        plt.grid()
        plt.show()
    except ValueError as e:
        messagebox.showerror("Error", f"Error en los datos: {e}")

# Función para resolver el sistema Lotka-Volterra
def solve_lotka_volterra():
    try:
        alpha = float(alpha_var.get())
        beta = float(beta_var.get())
        delta = float(delta_var.get())
        gamma = float(gamma_var.get())
        x0 = float(x0_var.get())
        y0_lv = float(y0_lv_var.get())

        if not (0 <= alpha <= 1):
            raise ValueError("α fuera del rango [0, 1]")
        if not (0 <= beta <= 0.1):
            raise ValueError("β fuera del rango [0, 0.1]")
        if not (0 <= delta <= 0.1):
            raise ValueError("δ fuera del rango [0, 0.1]")
        if not (0 <= gamma <= 1):
            raise ValueError("γ fuera del rango [0, 1]")
        if not (0 < x0 <= 1000):
            raise ValueError("x₀ fuera del rango (0, 1000]")
        if not (0 < y0_lv <= 500):
            raise ValueError("y₀ fuera del rango (0, 500]")

        t_values, system_values = euler_system(
            lambda t, state: lotka_volterra(t, state, alpha, beta, delta, gamma),
            np.array([x0, y0_lv]), 0, 10, 0.1
        )

        # Extraer valores de presas (x) y depredadores (y)
        x_values, y_values = system_values[:, 0], system_values[:, 1]

        # Graficar resultados del sistema Lotka-Volterra
        plt.figure(figsize=(12, 6))
        plt.plot(t_values, x_values, label="Presas (Numérico)", marker='o')
        plt.plot(t_values, y_values, label="Depredadores (Numérico)", marker='x')
        plt.title("Sistema Lotka-Volterra: Dinámica de Depredador-Presa")
        plt.xlabel("Tiempo (t)")
        plt.ylabel("Población")
        plt.legend()
        plt.grid()
        plt.show()

        # Gráfica fase (presas vs depredadores)
        plt.figure(figsize=(8, 6))
        plt.plot(x_values, y_values, label="Trayectoria en Fase", color="purple")
        plt.title("Espacio de Fase: Lotka-Volterra")
        plt.xlabel("Presas (x)")
        plt.ylabel("Depredadores (y)")
        plt.grid()
        plt.show()
    except ValueError as e:
        messagebox.showerror("Error", f"Error en los datos: {e}")

# Crear interfaz gráfica
root = Tk()
root.title("Resolución de EDO y Sistema Lotka-Volterra")

# Pantalla principal
main_frame = Frame(root)
main_frame.pack()

Label(main_frame, text="Seleccione qué desea resolver:").pack()
Button(main_frame, text="EDO de Primer Orden", command=show_edo_screen).pack(pady=5)
Button(main_frame, text="Sistema Lotka-Volterra", command=show_lotka_volterra_screen).pack(pady=5)

# Pantalla de EDO
edo_frame = Frame(root)
Label(edo_frame, text="Resolución de EDO de Primer Orden").grid(row=0, column=0, columnspan=2)
Label(edo_frame, text="Tasa de decaimiento (k) [0.1, 5]:").grid(row=1, column=0)
k_var = StringVar()
Entry(edo_frame, textvariable=k_var).grid(row=1, column=1)
Label(edo_frame, text="Valor inicial (y₀) (0, 100]:").grid(row=2, column=0)
y0_var = StringVar()
Entry(edo_frame, textvariable=y0_var).grid(row=2, column=1)
Button(edo_frame, text="Resolver EDO", command=solve_edo).grid(row=3, column=0, columnspan=2)

# Pantalla de Lotka-Volterra
lv_frame = Frame(root)
Label(lv_frame, text="Sistema Lotka-Volterra").grid(row=0, column=0, columnspan=2)
Label(lv_frame, text="α (crecimiento de presas) [0, 1]:").grid(row=1, column=0)
alpha_var = StringVar()
Entry(lv_frame, textvariable=alpha_var).grid(row=1, column=1)
Label(lv_frame, text="β (tasa de depredación) [0, 0.1]:").grid(row=2, column=0)
beta_var = StringVar()
Entry(lv_frame, textvariable=beta_var).grid(row=2, column=1)
Label(lv_frame, text="δ (reproducción de depredadores) [0, 0.1]:").grid(row=3, column=0)
delta_var = StringVar()
Entry(lv_frame, textvariable=delta_var).grid(row=3, column=1)
Label(lv_frame, text="γ (mortalidad de depredadores) [0, 1]:").grid(row=4, column=0)
gamma_var = StringVar()
Entry(lv_frame, textvariable=gamma_var).grid(row=4, column=1)
Label(lv_frame, text="Población inicial de presas (x₀) (0, 1000]:").grid(row=5, column=0)
x0_var = StringVar()
Entry(lv_frame, textvariable=x0_var).grid(row=5, column=1)
Label(lv_frame, text="Población inicial de depredadores (y₀) (0, 500]:").grid(row=6, column=0)
y0_lv_var = StringVar()
Entry(lv_frame, textvariable=y0_lv_var).grid(row=6, column=1)
Button(lv_frame, text="Resolver Lotka-Volterra", command=solve_lotka_volterra).grid(row=7, column=0, columnspan=2)

root.mainloop()
