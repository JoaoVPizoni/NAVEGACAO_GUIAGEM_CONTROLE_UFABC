# Simulador de Órbitas com TLEs e Integração Numérica

Este projeto realiza simulações da órbita de satélites em torno da Terra a partir de seus dados TLE (Two-Line Element), utilizando o modelo de dois corpos e integração numérica com `solve_ivp`. A órbita é visualizada em 3D juntamente com uma representação da Terra oblata.

## 📌 Funcionalidades

- Leitura e propagação de TLEs usando `sgp4`.
- Resolução da equação do problema dos dois corpos usando integração numérica (`solve_ivp`).
- Visualização em 3D da órbita e da Terra com escala realista.
- Cálculo do período orbital com base na energia orbital específica.

## 🚀 Requisitos

Crie um ambiente virtual e instale os pacotes com:

```bash
pip install -r requirements.txt