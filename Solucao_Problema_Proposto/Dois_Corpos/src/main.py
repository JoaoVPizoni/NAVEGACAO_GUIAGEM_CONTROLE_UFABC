from sgp4.api import Satrec, jday
from datetime import datetime
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os

# Constantes - MODELO WGS 84
mu = 3.986004418e14  # [m³/s²] - Terra
req = 6378.137       # [km]
rp = 6356.752        # [km]
f = (req - rp) / req # Achatamento
e = np.sqrt(1 - (rp/req)**2) #Excentricidade

# --------------------------------------
# Função para obter posição e velocidade de um TLE usando sgp4
def obter_estado_tle(linha1, linha2):
    sat = Satrec.twoline2rv(linha1, linha2)
    ano, mes, dia, hora, minuto, segundo = 2025, 5, 30, 0, 0, 0  # data de referência
    jd, fr = jday(ano, mes, dia, hora, minuto, segundo)
    e, r, v = sat.sgp4(jd, fr)
    if e != 0:
        raise RuntimeError("Erro ao processar TLE")
    return np.array(r) * 1000, np.array(v) * 1000  # posição e velocidade [m, m/s]

# --------------------------------------
# EDOs do problema de dois corpos
def equacoes_orbitais(t, x):
    r_vec = x[:3]
    v_vec = x[3:]
    r = np.linalg.norm(r_vec)

    # Termo de J2 (achatamento da Terra)
    J2 = 1.08263e-3
    R = req * 1000  # converte para metros
    x_, y_, z_ = r_vec

    fator_J2 = (3/2) * J2 * mu * R**2 / r**5
    ax_J2 = fator_J2 * x_ * (1 - 5 * (z_**2) / r**2)
    ay_J2 = fator_J2 * y_ * (1 - 5 * (z_**2) / r**2)
    az_J2 = fator_J2 * z_ * (3 - 5 * (z_**2) / r**2)

    a_vec = -mu * r_vec / r**3 + np.array([ax_J2, ay_J2, az_J2])

    return np.concatenate((v_vec, a_vec))


# --------------------------------------
# Resolução via solve_ivp
def simular_orbita(estado_inicial, T, passo, titulo, nome_sat):
    num_pontos = int(np.ceil(T / passo)) + 1
    t_eval = np.linspace(0, T, num_pontos)
    sol = solve_ivp(equacoes_orbitais, [0, T], estado_inicial, t_eval=t_eval,
                    method='RK45', atol=1e-10, rtol=1e-8) #RANGE KUTTA 45
    plotar_orbita(sol, titulo, nome_sat)


# --------------------------------------
# Plotagem da órbita e da Terra
def plotar_orbita(sol, titulo, nome_sat):
    # Criar diretório se não existir
    pasta = "Solucao_Problema_Proposto/Dois_Corpos/figuras"
    os.makedirs(pasta, exist_ok=True)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = sol.y[0] / 1000
    y = sol.y[1] / 1000
    z = sol.y[2] / 1000
    ax.plot(x, y, z, label='Órbita')

    # Esfera oblata (Terra)
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    xe = req * np.cos(u) * np.sin(v)
    ye = req * np.sin(u) * np.sin(v)
    ze = rp * np.cos(v)
    ax.plot_surface(xe, ye, ze, color='b', alpha=0.3)

    igualar_escalas(ax)
    ax.set_title(titulo)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_zlabel("z [km]")
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()

    # Salvar figura
    nome_arquivo = os.path.join(pasta, f"{nome_sat.replace(' ', '_')}.png")
    plt.savefig(nome_arquivo, dpi=300)
    plt.close()
    print(f"Figura salva em: {nome_arquivo}")

# --------------------------------------
# Igualar escalas para evitar distorção visual
def igualar_escalas(ax):
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    x_range = abs(xlim[1] - xlim[0])
    y_range = abs(ylim[1] - ylim[0])
    z_range = abs(zlim[1] - zlim[0])

    max_range = max(x_range, y_range, z_range)

    x_middle = np.mean(xlim)
    y_middle = np.mean(ylim)
    z_middle = np.mean(zlim)

    ax.set_xlim3d([x_middle - max_range/2, x_middle + max_range/2])
    ax.set_ylim3d([y_middle - max_range/2, y_middle + max_range/2])
    ax.set_zlim3d([z_middle - max_range/2, z_middle + max_range/2])


# --------------------------------------
# Simular para um satélite
def executar_para_tle(nome, linha1, linha2):
    print(f"Simulando para: {nome}")
    r, v = obter_estado_tle(linha1, linha2)
    estado_inicial = np.concatenate((r, v))
    mod_r = np.linalg.norm(r)
    mod_v = np.linalg.norm(v)
    epsilon = (mod_v**2) / 2 - mu / mod_r
    a = -mu / (2 * epsilon)
    T = 2 * np.pi * np.sqrt(a**3 / mu)  # período orbital
    print(f"Período orbital (s): {T:.2f}")
    simular_orbita(estado_inicial, T, T * 0.01, f"Órbita: {nome}", nome)


# --------------------------------------
# TLEs dos satélites
satelites = {
    "ISS": (
        "1 25544U 98067A   25150.54603503  .00012878  00000-0  23439-3 0  9999",
        "2 25544  51.6399  34.4830 0002197 166.0379 262.3840 15.49859072512427"
    ),
    "CBERS-4A": (
        "1 44883U 19093E   25150.83982066  .00001032  00000-0  13795-3 0  9998",
        "2 44883  97.7932 224.3455 0001774   6.2371 353.8864 14.81582034294467"
    ),
    "MOLNIYA 1-91": (
        "1 25485U 98054A   25150.46666369 -.00000126  00000-0  00000-0 0  9999",
        "2 25485  64.4919 341.5937 6788400 287.6338  12.7988  2.36440192204520"
    ),
    "STARONE D2": (
        "1 49055U 21069A   25150.42901189 -.00000270  00000-0  00000+0 0  9997",
        "2 49055   0.0110 322.8887 0001735 134.5547 235.1911  1.00271589 14066"
    )
}

# --------------------------------------
# Rodar simulações para todos
if __name__ == "__main__":
    for nome, (l1, l2) in satelites.items():
        executar_para_tle(nome, l1, l2)

#---------------------------------------
# Fluxograma de execução
# ┌────────────────────────────┐
# │ Início (main)              │
# └────────────┬───────────────┘
#              │
#              ▼
# ┌────────────────────────────┐
# │ Para cada satélite         │
# │ (nome, linha1, linha2)     │
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ obter_estado_tle()         │
# │ - Converte TLE para        │
# │   posição e velocidade     │
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ Calcular:                  │
# │ - módulo de r e v          │
# │ - energia específica       │
# │ - semi-eixo maior (a)      │
# │ - período orbital (T)      │
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ simular_orbita()           │
# │ - define t_eval            │
# │ - resolve EDO com solve_ivp│
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ equacoes_orbitais()        │
# │ - calcula aceleração       │
# │   gravitacional (-μr/r³)   │
# │ - (opcional: termo J2)     │
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ plotar_orbita()            │
# │ - Plota órbita 3D          │
# │ - Plota Terra oblata       │
# │ - Salva como imagem        │
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ Fim do loop                │
# │ (volta ao próximo satélite)│
# └────────────┬───────────────┘
#              ▼
# ┌────────────────────────────┐
# │ Fim do programa            │
# └────────────────────────────┘

