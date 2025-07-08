import matplotlib.pyplot as plt
import numpy as np
import os

# Constantes - MODELO WGS 84
mu = 3.986004418e14  # [m³/s²] - Terra
req = 6378.137       # [km]
rp = 6356.752        # [km]
f = (req - rp) / req # Achatamento
e = np.sqrt(1 - (rp/req)**2) #Excentricidade

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
xe = req * np.cos(u) * np.sin(v)
ye = req * np.sin(u) * np.sin(v)
ze = rp * np.cos(v)
ax.plot_surface(xe, ye, ze, color='b', alpha=0.65)
ax.set_title("Terra WGS 84")
ax.set_xlabel("x [km]")
ax.set_ylabel("y [km]")
ax.set_zlabel("z [km]")
ax.set_box_aspect([1, 1, 1])
plt.tight_layout()


# Salvar figura
pasta = "Solucao_Problema_Proposto/Dois_Corpos/figuras"
os.makedirs(pasta, exist_ok=True)
nome_arquivo = os.path.join(pasta, f"TERRA_WGS84.png")
plt.savefig(nome_arquivo, dpi=300)
print(f"Figura salva em: {nome_arquivo}")

plt.show()