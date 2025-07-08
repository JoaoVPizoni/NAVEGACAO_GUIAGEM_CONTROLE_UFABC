import numpy as np
from math import sqrt, radians, cos, sin, acos, asin
import datetime
import math
import unittest

class LocalizaPonto: #classe que contém as funções para geolocalizar um ponto
    def __init__(self):
        # Parâmetros do elipsoide GSW84
        self.a = 6378137.0  # semi-eixo maior [m]
        self.b = 6356752.0  # semi-eixo menor [m]
        self.e = sqrt((self.a**2 - self.b**2) / self.a**2)  # excentricidade
    
    def gms_para_decimal(self, grau, minutos, segundos, direcao):
        """Converte coordenadas de grau, minutos, segundos para decimal"""
        decimal = grau + minutos/60 + segundos/3600
        if direcao in ['S', 'W']:
            decimal *= -1
        return decimal
    
    def geografica_para_terrestre(self, lat_grau, lat_min, lat_seg, lat_dir, 
                                lon_grau, lon_min, lon_seg, lon_dir, alt):
        """Converte coordenadas geográficas para terrestres"""
        # Converter DMS para decimal
        phi = self.gms_para_decimal(lat_grau, lat_min, lat_seg, lat_dir)
        lam = self.gms_para_decimal(lon_grau, lon_min, lon_seg, lon_dir)
        
        # Converter para radianos
        phi_rad = radians(phi)
        lam_rad = radians(lam)
        
        # Calcular raio de curvatura
        N = self.a / sqrt(1 - (self.e**2 * sin(phi_rad)**2))
        
        # Calcular coordenadas terrestres
        x = (N + alt) * cos(phi_rad) * cos(lam_rad)
        y = (N + alt) * cos(phi_rad) * sin(lam_rad)
        z = (N * (1 - self.e**2) + alt) * sin(phi_rad)
        
        return np.array([x, y, z])
    
    def terrestre_para_inercial(self, r_terrestre, data, longitude):
        """Converte coordenadas terrestres para inerciais"""
        # Calcular tempo sideral de Greenwich
        theta_g = self.tempo_sideral_greenwich(data)
        
        # Matriz de rotação
        R = np.array([
            [cos(theta_g), sin(theta_g), 0],
            [-sin(theta_g), cos(theta_g), 0],
            [0, 0, 1]
        ])
        
        # Aplicar rotação
        r_inercial = R.T @ r_terrestre
        
        return r_inercial
    
    import math

    def tempo_sideral_greenwich(self, data):
        """Calcula o tempo sideral de Greenwich em radianos para uma data/hora específica"""
        
        # Extrair componentes da data/hora
        ano = data.year
        mes = data.month
        dia = data.day
        hora = data.hour
        minuto = data.minute
        segundo = data.second
        
        # Ajustar mês e ano se necessário
        if mes <= 2:
            ano -= 1
            mes += 12
        
        # Calcular parte inteira do dia juliano
        parte_a = ano // 100
        parte_b = 2 - parte_a + (parte_a // 4)
        
        # Calcular dia juliano (formula completa)
        termo1 = 365.25 * (ano + 4716)
        termo2 = 30.6001 * (mes + 1)
        parte_inteira = int(termo1) + int(termo2) + dia + parte_b - 1524
        
        # Calcular fração do dia
        fração_dia = (hora + minuto/60.0 + segundo/3600.0) / 24.0
        
        # Dia juliano completo
        dia_juliano = parte_inteira + fração_dia - 0.5
        
        # Séculos julianos desde J2000.0
        T = (dia_juliano - 2451545.0) / 36525.0
        
        # Calcular tempo sideral de Greenwich em grau
        termo_constante = 280.46061837
        termo_linear = 360.98564736629 * (dia_juliano - 2451545.0)
        termo_quadratico = 0.000387933 * (T * T)
        termo_cubico = (T * T * T) / 38710000.0
        
        theta_g = termo_constante + termo_linear + termo_quadratico - termo_cubico
        
        # Normalizar para o intervalo [0, 360) grau
        theta_g = theta_g - 360.0 * (theta_g // 360)
        if theta_g < 0:
            theta_g += 360.0
        
        # Converter para radianos
        return theta_g * (math.pi / 180.0)
    

    def processar_pontos(self, pontos):
        """Processa uma lista de pontos e retorna os resultados"""
        resultados = []
        
        for ponto in pontos:
            # Extrair dados do ponto
            dados_lat = ponto['latitude']
            dados_lon = ponto['longitude']
            altura = ponto['altitude']
            data_hora = ponto['datetime']
            
            # Converter para coordenadas terrestres
            r_terrestre = self.geografica_para_terrestre(
                dados_lat['grau'], dados_lat['min'], dados_lat['seg'], dados_lat['dir'],
                dados_lon['grau'], dados_lon['min'], dados_lon['seg'], dados_lon['dir'],
                altura
            )
            
            # Converter longitude para decimal
            lon_decimal = self.gms_para_decimal(
                dados_lon['grau'], dados_lon['min'], dados_lon['seg'], dados_lon['dir']
            )
            
            # Converter para coordenadas inerciais
            r_inercial = self.terrestre_para_inercial(r_terrestre, data_hora, lon_decimal)
            
            # Calcular norma do vetor terrestre
            norma_terrestre = math.sqrt(
                r_terrestre[0]**2 + r_terrestre[1]**2 + r_terrestre[2]**2
            )
            r_terrestre_unitario = [
                r_terrestre[0] / norma_terrestre,
                r_terrestre[1] / norma_terrestre,
                r_terrestre[2] / norma_terrestre
            ]
            
            # Calcular norma do vetor inercial
            norma_inercial = math.sqrt(
                r_inercial[0]**2 + r_inercial[1]**2 + r_inercial[2]**2
            )
            r_inercial_unitario = [
                r_inercial[0] / norma_inercial,
                r_inercial[1] / norma_inercial,
                r_inercial[2] / norma_inercial
            ]
            
            resultados.append({
                'terrestre': r_terrestre,
                'terrestre_unitario': r_terrestre_unitario,
                'inercial': r_inercial,
                'inercial_unitario': r_inercial_unitario
            })
        
        return resultados
    
    def inercial_para_terrestre(self, r_inercial, data, longitude):
        """Converte coordenadas inerciais de volta para terrestres: teste"""
        theta_g = self.tempo_sideral_greenwich(data)
        R = np.array([
            [cos(theta_g), sin(theta_g), 0],
            [-sin(theta_g), cos(theta_g), 0],
            [0, 0, 1]
        ])
        return R @ r_inercial  # Multiplicação direta (sem transposta)




if __name__ == "__main__":
    # Definir pontos da atividade
    pontos = [
        {
            'latitude': {'grau': 5, 'min': 55, 'seg': 23, 'dir': 'S'},
            'longitude': {'grau': 35, 'min': 9, 'seg': 51, 'dir': 'W'},
            'altitude': 39,
            'datetime': datetime.datetime(2025, 6, 3, 15, 54, 10)
        },
        {
            'latitude': {'grau': 2, 'min': 20, 'seg': 20, 'dir': 'S'},
            'longitude': {'grau': 44, 'min': 24, 'seg': 18, 'dir': 'W'},
            'altitude': 44,
            'datetime': datetime.datetime(2024, 7, 10, 11, 23, 10)
        },
        {
            'latitude': {'grau': 23, 'min': 40, 'seg': 37, 'dir': 'S'},
            'longitude': {'grau': 46, 'min': 33, 'seg': 46, 'dir': 'W'},
            'altitude': 778,
            'datetime': datetime.datetime(2025, 6, 24, 21, 45, 25)
        }
    ]
    
    localizador = LocalizaPonto()
    resultados = localizador.processar_pontos(pontos)
    
    def format_array(arr, decimals=6):
        return [f"{float(x):.{decimals}f}" for x in arr]
    
    # Exibir resultados
    for i, resultado in enumerate(resultados, 1):
        print(f"\nPonto {i}:")
        print("Coordenadas terrestres [m]:", format_array(resultado['terrestre']))
        print("Versor terrestre:", format_array(resultado['terrestre_unitario']))
        print("Coordenadas inerciais [m]:", format_array(resultado['inercial']))
        print("Versor inercial:", format_array(resultado['inercial_unitario']))
        
        # Teste de conversão inversa para cada ponto
        r_terrestre_recuperado = localizador.inercial_para_terrestre(
            resultado['inercial'],
            pontos[i-1]['datetime'],
            localizador.gms_para_decimal(
                pontos[i-1]['longitude']['grau'],
                pontos[i-1]['longitude']['min'],
                pontos[i-1]['longitude']['seg'],
                pontos[i-1]['longitude']['dir']
            )
        )
        print(f"Ponto {i} - terrestre recuperado:", format_array(r_terrestre_recuperado))
        print(f"Diferença terrestre original vs recuperado:", 
              format_array(np.array(resultado['terrestre']) - np.array(r_terrestre_recuperado)))