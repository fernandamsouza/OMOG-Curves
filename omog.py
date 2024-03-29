import math
import sys
from turtle import down, left
import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
import scipy.special

v = sp.symbols('v')

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

def getKnotsVetor(control_points, spline_deg):
    vetorKnot = []
    for i in range(spline_deg):
        vetorKnot.append(0)

    m = spline_deg + len(control_points) + 1 - 2 * pNURBS

    for i in range(m):
        vetorKnot.append(i / (m - 1))

    for i in range(spline_deg):
        vetorKnot.append(1)

    return vetorKnot

def projetarPontosControle(pontosControle, funcaoLambdaMul):
    novosPontosControle = []

    for pt in pontosControle:
        novoPonto = (funcaoLambdaMul(pt[0], 1), funcaoLambdaMul(pt[1], 1), 1)
        novosPontosControle.append(novoPonto)

    return novosPontosControle

def de_boor(indexVetorKnot: int, x: int, vetorKnotPos, vetorNovosPontosDeControle, grauNURBS: int):

    vetorNovosPontosDeControle = list(map(list, vetorNovosPontosDeControle))
    arrayCurvaNurbs = [vetorNovosPontosDeControle[j + indexVetorKnot - grauNURBS] for j in range(0, grauNURBS+1)]

    funcaoSoma = lambda x,y: x + y

    linha = []
    for r in range(1, grauNURBS+1):
        if len(linha) > 1:
            linha = []
        for j in range(grauNURBS, r-1, -1):
            alpha = (x - vetorKnotPos[j+indexVetorKnot-grauNURBS]) / (vetorKnotPos[j+1+indexVetorKnot-r] - vetorKnotPos[j+indexVetorKnot-grauNURBS])
            a = list(map(lambda x: x * (1.0 - alpha), arrayCurvaNurbs[j-1]))
            b = list(map(lambda x: x * alpha, arrayCurvaNurbs[j]))
            arrayCurvaNurbs[j] = list(map(funcaoSoma, a, b))
            linha.append(tuple(arrayCurvaNurbs[j]))
    assert (len(linha) == 1)

    return arrayCurvaNurbs[grauNURBS]

def nurbs(vetorKnotCortado, numeroPontosCurva):

    xComeco = 0
    xFinal = 1
    tamanhoPasso = (xFinal - xComeco) / numeroPontosCurva

    for step_num in range(numeroPontosCurva):

        xCurva = xComeco + step_num * tamanhoPasso
        indexKnot = pNURBS - 1

        for pontoKnot in vetorKnotCortado[1:len(vetorKnotCortado)]:
            indexKnot += 1
            if pontoKnot >= xCurva:
                break

        funcOpMul = lambda x, y: x * y
        novosPontosDeControle = projetarPontosControle(arrayDePontosNURBS, funcOpMul)
        retorno = de_boor(indexKnot, xCurva, vetorKnot, novosPontosDeControle, pNURBS)

        arrayEixoXNURBS.append(retorno[0])
        arrayEixoYNURBS.append(retorno[1])
        
    return arrayEixoXNURBS, arrayEixoYNURBS
 
def rotacionarBezierNurbs(arrayDePontosNURBS, arrayDePontos):

    # C0
    ultimoPontoNURBS = arrayDePontosNURBS[4]
    primeiroPontoBezier = arrayDePontos[0]
    deltaX = ultimoPontoNURBS[0] - primeiroPontoBezier[0]
    deltaY = ultimoPontoNURBS[1] - primeiroPontoBezier[1]

    for ponto in arrayDePontos:
        ponto[0] = ponto[0] + deltaX
        ponto[1] = ponto[1] + deltaY

    # G1
    penultimoNURBS = arrayDePontosNURBS[3]
    primeiroBEZIER = arrayDePontos[0]

    deltaXC1 = penultimoNURBS[0] - primeiroBEZIER[0]
    deltaYC1 = penultimoNURBS[1] - primeiroBEZIER[1]
    
    arrayDePontos[1][0] = primeiroBEZIER[0] - deltaXC1
    arrayDePontos[1][1] = primeiroBEZIER[1] - deltaYC1

    # G2
    antiPenultimoNURBS = arrayDePontosNURBS[2]
    primeiroBEZIER = arrayDePontos[0]

    deltaXC1 = antiPenultimoNURBS[0] - primeiroBEZIER[0]
    deltaYC1 = antiPenultimoNURBS[1] - primeiroBEZIER[1]
    
    arrayDePontos[2][0] = primeiroBEZIER[0] - deltaXC1
    arrayDePontos[2][1] = primeiroBEZIER[1] - deltaYC1

    # deltaX = arrayDePontosNURBS[2][0] - arrayDePontosNURBS[1][0]
    # deltaY = arrayDePontosNURBS[2][1] - arrayDePontosNURBS[1][1]
    # seclast_pos = arrayDePontos[-2]

    # arrayDePontos[-3] = seclast_pos[0] - deltaX, seclast_pos[1] - deltaY

    return arrayDePontos

def b (i, numeroPontos):
    return scipy.special.binom(numeroPontos - 1, i) * (v ** i) * ((1 - v) ** (numeroPontos - i - 1))

def c1 (arrayDePontos):
    numeroPontos = len(arrayDePontos)
    binomio = [b(i, numeroPontos) for i in range(numeroPontos)]
    inicializacao = np.linspace(0, 1, 1000)

    bezierXC1 = []
    bezierYC1 = []
    tuplaFinal = []

    xf = 0
    yf = 0

    for i in range(numeroPontos):
        xf += arrayDePontos[i][0] * binomio[i]
        yf += arrayDePontos[i][1] * binomio[i]
    
    for j in inicializacao:
        x = xf.subs(v, j)
        y = yf.subs(v, j)

        bezierXC1.append(x)
        bezierYC1.append(y)
         
        tuplaFinal.append([x,y])

    return tuplaFinal, bezierXC1, bezierYC1, xf

if __name__ == "__main__":

    # NURBS Grau 4
    arrayEixoXNURBS = []
    arrayEixoYNURBS = []

    # Definição do número de pontos de controle para NURBS
    numeroPontosDeControleNURBS = 5

    # Seleção do usuário por click dos pontos para formação da curva NURBS
    plt.title("NURBS - Clique em 5 (cinco) pontos")

    # Obtém a array de pontos selecionados pelo usuário
    arrayDePontosNURBS = np.array(plt.ginput(numeroPontosDeControleNURBS))
    
    # Divide a array em eixo X e eixo Y
    arrayEixoXNURBSPontos = arrayDePontosNURBS[:,0]
    arrayEixoYNURBSPontos = arrayDePontosNURBS[:,1]

    # Grau da NURBS
    pNURBS = 4

    vetorKnot = getKnotsVetor(arrayDePontosNURBS, pNURBS)
    vetorKnotCortado = [vetorKnot[j] for j in range(pNURBS, len(vetorKnot) - pNURBS)]

    arrayEixoXNURBS, arrayEixoYNURBS = nurbs(vetorKnotCortado, 1000)

    # Fecho o plot de seleção dos pontos para abrir outro com a curva
    plt.close()

    # Plot da curva na tela - NURBS
    plt.figure(figsize=(8, 6), dpi=80)
    plt.xlim(right=1)
    plt.ylim(top=1)
    plt.plot(arrayEixoXNURBS, arrayEixoYNURBS, color = 'green')
    plt.plot(arrayEixoXNURBSPontos, arrayEixoYNURBSPontos, 'o', color = 'yellow')

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

    # Obtendo continuidade C0 e C1
    arrayDePontos = rotacionarBezierNurbs(arrayDePontosNURBS, arrayDePontos)

    arrayCurvaBezEixoX, arrayCurvaBezEixoY = bezier(arrayDePontos, 1000)

    # Divide a array em eixo X e eixo Y
    arrayEixoX = arrayDePontos[:,0]
    arrayEixoY = arrayDePontos[:,1]

    # Fecho o plot para gerar as curvas com continuidade C0
    plt.close()

    plt.plot(arrayEixoXNURBS, arrayEixoYNURBS, color = 'green')
    plt.plot(arrayEixoXNURBSPontos, arrayEixoYNURBSPontos, 'o', color = 'yellow')
    plt.plot(arrayCurvaBezEixoX, arrayCurvaBezEixoY, color = 'blue')
    plt.plot(arrayEixoX, arrayEixoY, 'o', color = 'red')
    plt.title("Bezier Grau 3 + Nurbs Grau 4 - G2 - Fernanda Maria de Souza")
    plt.show()