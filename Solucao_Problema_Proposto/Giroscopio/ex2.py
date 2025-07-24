import numpy as np
import matplotlib.pyplot as plt

# Dados experimentais
dados = {
    'Massa (g)': [30, 30, 30, 30, 30],
    'Voltas': [39, 91, 44, 87, 25],
    'Tempo (s)': [10, 10, 10, 10, 10],
    'Volta_precessao': [0.5, 0.5, 0.5, 0.5, 0.5],
    'Tempo_precessao (s)': [8.75, 20, 9.56, 18, 6.09]
}

massa = np.array(dados['Massa (g)']) / 1000  # kg
voltas = np.array(dados['Voltas'])
tempo = np.array(dados['Tempo (s)'])
volta_precessao = np.array(dados['Volta_precessao'])
tempo_precessao = np.array(dados['Tempo_precessao (s)'])

# Parâmetros fixos
g = 9.81  # m/s²
r_estrela = 0.27  # m
m_estrela = 0.03  # kg

# Cálculos das velocidades angulares
omega_R = (2 * np.pi * voltas) / tempo  # rad/s
T_p = tempo_precessao / volta_precessao  # s
omega_p = (2 * np.pi) / T_p  # rad/s

# Linearização: ω_p vs 1/ω_R
x = 1 / omega_R  # 1/ω_R (variável independente)
y = omega_p      # ω_p (variável dependente)

# =============================================
# IMPLEMENTAÇÃO DOS MÍNIMOS QUADRADOS
# =============================================
n = len(x)
sum_x = np.sum(x)
sum_y = np.sum(y)
sum_xy = np.sum(x * y)
sum_x2 = np.sum(x**2)

# Cálculo dos coeficientes da reta y = a*x + b
a = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)  # Coef. angular
b = (sum_y - a * sum_x) / n                                 # Coef. linear

# Cálculo do momento de inércia
I_p = (m_estrela * g * r_estrela) / a

# Cálculo do coeficiente de determinação R²
y_medio = np.mean(y)
SST = np.sum((y - y_medio)**2)          # Soma total dos quadrados
SSR = np.sum((a * x + b - y_medio)**2)  # Soma dos quadrados da regressão
R2 = SSR / SST                          # Coeficiente de determinação

# Cálculo dos resíduos
y_pred = a * x + b
residuos = y - y_pred

# =========================
# SAÍDAS NO CONSOLE
# =========================
print("===== RESULTADOS DO AJUSTE MANUAL =====")
print(f"Coeficiente angular (a): {a:.6f} rad²/s²")
print(f"Coeficiente linear (b): {b:.6f} rad/s")
print(f"Coeficiente de determinação (R²): {R2:.6f}")
print(f"Momento de inércia (I_p): {I_p:.6f} kg·m²")

print("\n===== TABELA DE DADOS E AJUSTE =====")
print("| Exp | 1/ω_R (s/rad) | ω_p (rad/s) | ω_p_ajustado | Resíduo  |")
print("|-----|----------------|-------------|--------------|----------|")
for i in range(n):
    print(f"| {i+1:2d} | {x[i]:14.4f} | {y[i]:11.3f} | {y_pred[i]:11.3f} | {residuos[i]:8.3f} |")

# =============================================
# GRÁFICOS
# =============================================
plt.figure()
plt.scatter(x, y, color='blue', label='Dados experimentais')
plt.plot(x, a * x + b, 'r--', label=f'Ajuste: ω_p = {a:.3f}/ω_R + {b:.3f}')
plt.xlabel('1/ω_R (s/rad)')
plt.ylabel('ω_p (rad/s)')
plt.title('Ajuste Linear por Mínimos Quadrados')
plt.legend()
plt.grid(True)


plt.tight_layout()
plt.show()