import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path

# Cria a estrutura de diretórios
output_dir = Path("Solucao_Problema_Proposto") / "Pouso_na_Lua" / "figuras"
output_dir.mkdir(parents=True, exist_ok=True)

# Condições iniciais do V/E
theta = np.array([-30, 0, 30])  # [°]
g = 1.3                         # Gravidade da Lua [m/s^2]
m = 1                           # Massa do V/E [kg]
v0x = 2                         # Velocidade inicial horizontal [m/s]
v0y = -7                        # Velocidade inicial vertical [m/s]
x0 = 0                          # Posição inicial horizontal [m]
y0 = 70                      # Posição inicial vertical [m]
E = 2                           # Empuxo do motor foguete [N]
vy = np.arange(-15, 0, 0.001)   # Vetor de velocidades verticais para o envelope de voo

# Função para calcular as alturas dos envelopes de voo
def y_theta(theta):
    return (vy**2) / (2 * (E * np.cos(np.radians(theta)) / m - g))

# Cálculo das curvas de voo para os ângulos 0 e 30
y_theta0 = y_theta(theta[1])
y_theta30 = y_theta(theta[2])

# Plots dos Envelopes de Voo
plt.figure(figsize=(10, 6))
plt.plot(vy, y_theta0, 'r', linewidth=1.2)
plt.plot(vy, y_theta30, 'b', linewidth=1.2)
plt.plot(v0y, y0, 'k*', markersize=10)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('Velocidade vertical [m/s]')
plt.ylabel('Altitude [m]')
plt.legend(['θ = 0°', 'θ = 30°', 'Condição Inicial'])
plt.title('Envelopes de Voo')
plt.savefig(output_dir / "envelopes_de_voo.png")  # Salva a figura
plt.show()

# Função para cálculo das estratégias de pouso
def landing_strategy(E, g, m, theta, v0y, y0, strategy):
    a0 = 2 * (E * np.cos(np.radians(theta[1])) / m - g)
    a30 = 2 * (E * np.cos(np.radians(theta[2])) / m - g)
    
    if strategy == 1:
        vy1 = -np.sqrt((a0 * a30 / (a30 - a0)) * (v0y**2 / a0 - y0))
        y1 = vy1**2 / a30
        ty1 = (vy1 - v0y) / (E * np.cos(np.radians(theta[1])) / m - g)
        ty2 = -vy1 / (E * np.cos(np.radians(theta[2])) / m - g)
    elif strategy == 2:
        vy1 = -np.sqrt((a0 * a30 / (a0 - a30)) * (v0y**2 / a30 - y0))
        y1 = vy1**2 / a0
        ty1 = (vy1 - v0y) / (E * np.cos(np.radians(theta[1])) / m - g)
        ty2 = -vy1 / (E * np.cos(np.radians(theta[2])) / m - g)
    
    tytotal = ty1 + ty2
    vy1_str = np.linspace(v0y, vy1, 200)
    vy2_str = np.linspace(vy1, 0, 200)
    y1_str = (vy1_str**2 - v0y**2 + a0 * y0) / a0
    y2_str = (vy2_str**2) / (2 * (E * np.cos(np.radians(theta[2])) / m - g))
    
    return vy1, y1, ty1, ty2, tytotal, vy1_str, vy2_str, y1_str, y2_str

# Estratégia de pouso suave na vertical
print('================================')
print('Análise Movimento Eixo Y')

strategy = 1  # Defina a estratégia de pouso aqui
vy1, y1, ty1, ty2, tytotal, vy1_str, vy2_str, y1_str, y2_str = landing_strategy(E, g, m, theta, v0y, y0, strategy)

print(f'Vy_intersec = {vy1:.4f} m/s')
print(f'y_intersec = {y1:.4f} m')
print(f't_intersec = {ty1:.4f} s')
print(f't_curva = {ty2:.4f} s')
print(f't_voo = {tytotal:.4f} s')

# Plots das Estratégias de Pouso Vertical
plt.figure(figsize=(10, 6))
plt.plot(vy, y_theta30, 'b', linewidth=2)
plt.plot(vy, y_theta0, 'r', linewidth=2)
plt.plot(v0y, y0, 'g*', markersize=10)
plt.plot(vy1, y1, 'ko', markersize=10)
plt.plot(vy1_str, y1_str, 'k', linewidth=2)
plt.plot(vy2_str, y2_str, 'k', linewidth=2)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('Velocidade vertical [m/s]')
plt.ylabel('Altitude [m]')
plt.legend(['θ = 30°', 'θ = 0°', 'Condição Inicial', 'Interseção Parábolas'])
plt.title('Estratégias de Pouso Vertical')
plt.savefig(output_dir / "estrategias_pouso_vertical.png")  # Salva a figura
plt.show()

# Função para análise do movimento no eixo X
def analyze_horizontal_movement(E, m, theta, v0x, strategy, ty1, ty2, tmin):
    print('================================')
    print('Análise Movimento Eixo X')

    ax = -E * np.sin(np.radians(theta[0])) / m
    print(f'ax = {ax:.4f} m/s^2')
    print(f'tmin = {tmin:.4f} s')

    if strategy == 1:
        t1 = (ty2 - tmin) / 2
        t2 = t1 + tmin
        tstop = t1 + t2
        vx1 = v0x + ax * t1
        vx2 = vx1 - ax * t2
    elif strategy == 2:
        axa = -E * np.sin(np.radians(theta[2])) / m
        taa = (ty2 - tmin) / 2
        t1a = taa + tmin
        t2a = taa
        tstopa = t1a + t2a
        vx1a = v0x + axa * t1a
        vx2a = vx1a - axa * t2a

    print(f'Instante de desligamento motor = {tstop:.4f} s')
    print(f'Tempo total de motores ligados = {tstop:.4f} s')

# Análise do movimento no eixo X
tmin = abs(v0x / (-E * np.sin(np.radians(theta[0])) / m))
analyze_horizontal_movement(E, m, theta, v0x, strategy, ty1, ty2, tmin)