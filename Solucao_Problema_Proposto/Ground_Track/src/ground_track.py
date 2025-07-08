import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, radians, degrees, atan2, asin, pi
import datetime
from sgp4.api import Satrec, jday
import os

class CalculadorRastroTerrestre:
    def __init__(self):
        # Constantes físicas
        self.G = 6.67430e-11  # constante gravitacional [m^3/kg/s^2]
        self.M_terra = 5.972e24  # massa da Terra [kg]
        self.mu = self.G * self.M_terra
        self.omega_terra = 7.2921150e-5  # velocidade angular da Terra [rad/s]
        
        # Parâmetros do elipsoide terrestre (GRS80)
        self.a = 6378137.0  # semi-eixo maior [m]
        self.b = 6356752.0  # semi-eixo menor [m]
    
    def obter_estado_do_tle(self, linha1, linha2, ano, mes, dia, hora, minuto, segundo):
        """Obtém posição e velocidade a partir de TLE usando SGP4"""
        sat = Satrec.twoline2rv(linha1, linha2)
        jd, fr = jday(ano, mes, dia, hora, minuto, segundo)
        e, r, v = sat.sgp4(jd, fr)
        if e != 0:
            raise RuntimeError("Erro ao processar TLE")
        return np.array(r) * 1000, np.array(v) * 1000  # converte km para m
    
    def dinamica_orbital(self, t, estado):
        """Equações diferenciais para movimento orbital"""
        x, y, z, vx, vy, vz = estado
        r = sqrt(x**2 + y**2 + z**2)
        
        dxdt = vx
        dydt = vy
        dzdt = vz
        dvxdt = -self.mu * x / r**3
        dvydt = -self.mu * y / r**3
        dvzdt = -self.mu * z / r**3
        
        return [dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt]
    
    def calcular_orbita(self, r0, v0, periodo, n_pontos=1000):
        """Calcula órbita por um período"""
        estado0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]]
        t_eval = np.linspace(0, periodo, n_pontos)
        sol = solve_ivp(self.dinamica_orbital, [0, periodo], estado0, 
                        t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)
        return sol.y.T, t_eval
    
    def dia_juliano(self, dt):
        """Calcula Dia Juliano a partir de datetime"""
        a = (14 - dt.month) // 12
        y = dt.year + 4800 - a
        m = dt.month + 12*a - 3
        
        jd = dt.day + (153*m + 2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
        jd += (dt.hour - 12)/24 + dt.minute/1440 + dt.second/86400
        
        return jd
    
    def tempo_sideral_greenwich(self, jd):
        """Calcula Tempo Sideral de Greenwich"""
        T = (jd - 2451545.0) / 36525
        theta_g = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T**2 - T**3/38710000
        return theta_g % 360  # graus
    
    def calcular_rastro_terrestre(self, orbita, tempo_inicio):
        """Calcula coordenadas do rastro terrestre"""
        latitudes = []
        longitudes = []
        
        jd_inicio = self.dia_juliano(tempo_inicio)
        
        for i in range(len(orbita)):
            x, y, z = orbita[i,0], orbita[i,1], orbita[i,2]
            r = sqrt(x**2 + y**2 + z**2)
            
            # Calcula Dia Juliano atual
            jd = jd_inicio + (i/len(orbita)) * (1/1440)  # assumindo 1 período = 1 dia
            
            # Calcula GST em graus
            theta_g = self.tempo_sideral_greenwich(jd)
            
            # Calcula latitude e longitude
            lat = degrees(asin(z / r))
            lon = degrees(atan2(y, x)) - theta_g
            
            # Normaliza longitude para [-180, 180]
            while lon < -180:
                lon += 360
            while lon > 180:
                lon -= 360
            
            latitudes.append(lat)
            longitudes.append(lon)
        
        return latitudes, longitudes
    
    def plotar_rastro_terrestre(self, latitudes, longitudes, nome_sat, titulo="Rastro Terrestre"):
        """Plota o rastro terrestre"""
        plt.figure(figsize=(12, 6))
        
        # Plota rastro terrestre
        plt.scatter(longitudes, latitudes, s=5, c='b', alpha=0.5)
        
        # Formatação
        plt.title(titulo, fontsize=14)
        plt.xlabel("Longitude [°]", fontsize=12)
        plt.ylabel("Latitude [°]", fontsize=12)
        plt.xlim(-180, 180)
        plt.ylim(-90, 90)
        
        # Adiciona grade e linhas de referência
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.axhline(0, color='k', linestyle='-', linewidth=0.5)  # Equador
        plt.axvline(0, color='k', linestyle='-', linewidth=0.5)  # Meridiano de Greenwich
        
        # Adiciona linhas de latitude de referência
        for lat in [-60, -30, 30, 60]:
            plt.axhline(lat, color='gray', linestyle=':', linewidth=0.5)
        
        plt.tight_layout()
        
        # Salvar figura
        pasta = "Solucao_Problema_Proposto/Ground_Track/figuras"
        os.makedirs(pasta, exist_ok=True)
        nome_arquivo = os.path.join(pasta, f"ground_track_{nome_sat.replace(' ', '_')}.png")
        plt.savefig(nome_arquivo, dpi=300)
        print(f"Figura do rastro terrestre salva em: {nome_arquivo}")
        
        plt.show()

    def plotar_orbita_3d(self, orbita, nome_sat, titulo="Órbita 3D"):
        """Plota órbita 3D com a Terra"""
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plota órbita
        ax.plot(orbita[:,0]/1000, orbita[:,1]/1000, orbita[:,2]/1000, 'b', linewidth=1)
        
        # Plota Terra
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = self.a * np.outer(np.cos(u), np.sin(v)) / 1000
        y = self.a * np.outer(np.sin(u), np.sin(v)) / 1000
        z = self.b * np.outer(np.ones(np.size(u)), np.cos(v)) / 1000
        ax.plot_surface(x, y, z, color='lightblue', alpha=0.3)
        
        ax.set_xlabel('X [km]')
        ax.set_ylabel('Y [km]')
        ax.set_zlabel('Z [km]')
        ax.set_title(titulo)
        
        # Equaliza eixos
        max_range = np.array([orbita[:,0].max()-orbita[:,0].min(), 
                             orbita[:,1].max()-orbita[:,1].min(), 
                             orbita[:,2].max()-orbita[:,2].min()]).max() / 1000 / 2.0
        mid_x = (orbita[:,0].max()+orbita[:,0].min()) * 0.5 / 1000
        mid_y = (orbita[:,1].max()+orbita[:,1].min()) * 0.5 / 1000
        mid_z = (orbita[:,2].max()+orbita[:,2].min()) * 0.5 / 1000
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        # Salvar figura
        pasta = "Ground_Track/figuras"
        os.makedirs(pasta, exist_ok=True)
        nome_arquivo = os.path.join(pasta, f"orbita_3d_{nome_sat.replace(' ', '_')}.png")
        plt.savefig(nome_arquivo, dpi=300)
        print(f"Figura da órbita 3D salva em: {nome_arquivo}")
        
        plt.tight_layout()
        plt.show()

    def calcular_periodo_orbital(self, r0, v0):
        """Calcula período orbital a partir de posição e velocidade"""
        epsilon = (np.linalg.norm(v0)**2)/2 - self.mu/np.linalg.norm(r0)
        a = -self.mu/(2*epsilon)
        return 2 * pi * sqrt(a**3 / self.mu)


if __name__ == "__main__":
    # Inicializa calculador
    rt = CalculadorRastroTerrestre()
    
    # TLEs de satélites
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
    
    # Data de referência
    data_ref = datetime.datetime(2006, 9, 11, 0, 0, 0)
    
    # Processa cada satélite
    for nome_sat, (linha1, linha2) in satelites.items():
        print(f"\nProcessando {nome_sat}...")
        
        # Obtém estado inicial do TLE
        r0, v0 = rt.obter_estado_do_tle(linha1, linha2, 2023, 5, 30, 0, 0, 0)
        
        # Calcula período orbital
        periodo = rt.calcular_periodo_orbital(r0, v0)
        print(f"Período orbital: {periodo/3600:.2f} horas")
        
        # Calcula órbita
        orbita, t = rt.calcular_orbita(r0, v0, periodo)
        
        # Calcula rastro terrestre
        lats, lons = rt.calcular_rastro_terrestre(orbita, data_ref)
        
        # Plota resultados
        rt.plotar_orbita_3d(orbita, nome_sat, f"Órbita 3D - {nome_sat}")
        rt.plotar_rastro_terrestre(lats, lons, nome_sat, f"Rastro Terrestre - {nome_sat}")