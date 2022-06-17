from turtle import color
import numpy as np
from matplotlib import pyplot as plt

def bezier(arrayDePontos, numeroPontosCurva):
    grauDaBezier = len(arrayDePontos) - 1

    arrayPontosControleEixoX = [ponto[0] for ponto in arrayDePontos]
    arrayPontosControleEixoY = [ponto[1] for ponto in arrayDePontos]

    # Retorna len(numeroPontosCurva) números espaçados entre 0 e 1
    arrayInicializacaoCurva = np.linspace(0, 1, numeroPontosCurva)

    arrayCurvaBezierEixoX = [0] * numeroPontosCurva
    arrayCurvaBezierEixoY = [0] * numeroPontosCurva

    for k in range(0, len(arrayInicializacaoCurva)):
        for j in range(0, len(arrayDePontos)):
            # Cálculo do Binômio de Newton para a resolução dos coeficientes
            auxiliar = binomioNewton(j, grauDaBezier) * ((1 - arrayInicializacaoCurva[k]) ** (grauDaBezier - j)) * (arrayInicializacaoCurva[k] ** j)
            arrayCurvaBezierEixoX[k] = arrayCurvaBezierEixoX[k] + arrayPontosControleEixoX[j] * auxiliar
            arrayCurvaBezierEixoY[k] = arrayCurvaBezierEixoY[k] + arrayPontosControleEixoY[j] * auxiliar

    return arrayCurvaBezierEixoX, arrayCurvaBezierEixoY

def binomioNewton(i,n):
    return fatorial(n)/(fatorial(i)*fatorial(n-i))

def fatorial(n):
    nFat = 1
    for i in range(2, n+1):
        nFat = nFat * i
    return nFat

if __name__ == "__main__":
    # Definição do número de pontos de controle para a Bezier. Como minha curva é de grau 3, defini 4 pontos.
    numeroPontosDeControle = 4

    # Seleção do usuário por click dos pontos para formação da curva
    plt.title("Clique em 4 (quatro) pontos")

    # Obtém a array de pontos selecionados pelo usuário
    arrayDePontos = np.array(plt.ginput(numeroPontosDeControle))

    # Divide a array em eixo X e eixo Y
    arrayEixoX = arrayDePontos[:,0]
    arrayEixoY = arrayDePontos[:,1]

    # Chama a função para retornar a curva de Bezier
    # Obs.: o número de pontos na curva a torna proporcionalmente sinuosa (+ pontos + contínua; - pontos + quebrada)
    arrayCurvaBezEixoX, arrayCurvaBezEixoY = bezier(arrayDePontos, 1000)

    # Fecho o plot de seleção dos pontos para abrir outro com a curva
    plt.close()

    # Plot da curva na tela
    plt.plot(arrayCurvaBezEixoX, arrayCurvaBezEixoY, color = 'blue')
    plt.plot(arrayEixoX, arrayEixoY, 'o', color = 'red')
    plt.title("Bezier Grau 3 - Fernanda Maria de Souza")
    plt.show()