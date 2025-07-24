import numpy as np
import matplotlib.pyplot as plt

# =========================
# DADOS FIXOS DO EXPERIMENTO
# =========================
m = 30 / 1000            # massa suspensa em kg
g = 9.81                 # gravidade
diametro_carretel_mm = 50
r = (diametro_carretel_mm / 1000) / 2  # raio do carretel em metros

# =========================
# DADOS EXPERIMENTAIS
# =========================
h_cm = np.array([60, 50, 45, 40, 35, 30])       # altura em cm
tF_s = np.array([7.38, 6.53, 5.81, 5.38, 5.28, 4.88])  # tempo em s

# Conversão para unidades SI
h = h_cm / 100.0        # metros
t2 = tF_s ** 2          # segundos²

# ==============================================================
# (1) MÉTODO DA ACELERAÇÃO (equação 6)
# ==============================================================
a_vals = (2 * h) / t2   # a = 2h / t²
IP_a = m * ((g / a_vals) - 1) * r**2

# ==============================================================
# (2) MÉTODO DOS MÍNIMOS QUADRADOS – Equação (12)
# Ajuste linear de t² = a*h + b
# ==============================================================
n = len(h)
Sx = np.sum(h)
Sy = np.sum(t2)
Sxx = np.sum(h**2)
Sxy = np.sum(h * t2)

slope = (n * Sxy - Sx * Sy) / (n * Sxx - Sx**2)
intercept = (Sy - slope * Sx) / n

IP_ls = (slope * m * g * r**2 / 2) - m * r**2

# ==============================================================
# (3) MÉTODO DA ENERGIA (equação 9)
# Usando v² = 2 a h (Torricelli com a = 2h / t²)
# ==============================================================
v = (2 * h) / tF_s      # v = 2h / tF (Torricelli)
IP_energia = ((2 * g * h) / (v**2) - 1) * m * r**2

# =========================
# SAÍDAS NO CONSOLE
# =========================
print("===== MÉTODO DA ACELERAÇÃO (Equação 6) =====")
for i in range(len(h)):
    print(f"h = {h[i]:.2f} m, tF = {tF_s[i]:.2f} s -> a = {a_vals[i]:.4f} m/s² -> I_P = {IP_a[i]:.6f} kg·m²")

print("\n===== MÉTODO DO AJUSTE LINEAR (Equação 12) =====")
print(f"Massa m = {m:.3f} kg")
print(f"Raio do carretel r = {r:.4f} m")
print(f"Inclinação = {slope:.3f} s²/m")
print(f"Intercepto = {intercept:.3f} s²")
print(f"Momento de inércia I_P (ajuste) = {IP_ls:.6f} kg·m²")

print("\n===== VELOCIDADE USANDO TORRICELLI =====")
for i in range(len(h)):
    print(f"h = {h[i]:.2f} m, tF = {tF_s[i]:.2f} s -> v = {v[i]:.4f} m/s")

print("\n===== MÉTODO DA ENERGIA (Equação 9) =====")
for i, Ip in enumerate(IP_energia):
    print(f"h = {h[i]:.2f} m -> I_P = {Ip:.6f} kg·m²")

# =========================
# GRÁFICO DO AJUSTE
# =========================
x_fit = np.linspace(min(h), max(h), 100)
y_fit = slope * x_fit + intercept

plt.figure(figsize=(6,4))
plt.scatter(h, t2, label='Dados experimentais', color='blue')
plt.plot(x_fit, y_fit, label='Ajuste mínimos quadrados', color='red')
plt.xlabel('Altura h (m)')
plt.ylabel('Tempo² (s²)')
plt.title('t_F² vs h (Equação 12)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()