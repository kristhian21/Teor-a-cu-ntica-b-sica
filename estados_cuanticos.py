from Complex_Vector_Spaces import *
import numpy as np


def transpuesta(matriz):
    """
    :param matriz:
    :return matriz transpuesta:
    """
    if isinstance(matriz[0], tuple):
        resultado = [[None for j in range(1)] for i in range(len(matriz))]

        for i in range(len(resultado)):
            for j in range(1):
                resultado[i][j] = matriz[i]

    return resultado


def convertirMatrizNP(matriz):
    """
    :param matriz:
    :return matriz:
    """
    dim = len(matriz)
    nMatriz = [[None for j in range(dim)] for i in range(dim)]

    for i in range(dim):
        for j in range(dim):
            nMatriz[i][j] = complex(matriz[i][j][0], matriz[i][j][1])

    return nMatriz


def probabilidad_pos(posicion, ket):
    """
    :param posicion:
    :param ket:
    :return numero real:
    """
    numerador = moduloComplejos(ket[posicion][0]) ** 2
    suma = 0

    for i in range(len(ket)):
        for j in range(1):
            suma += moduloComplejos(ket[i][j]) ** 2

    return numerador / math.sqrt(suma) ** 2


def bra(ket):
    """
    :param ket:
    :return vector bra:
    """
    resultado = [0 for i in range(len(ket))]

    for i in range(len(resultado)):
        resultado[i] = conjugadoComplejos(ket[i][0])

    return resultado


def amplitudTransicion(vector_bra, ket):
    """
    :param vector_bra:
    :param ket:
    :return numero complejo:
    """
    resultado = (0, 0)
    for i in range(len(ket)):
        resultado = sumaComplejos(resultado, productoComplejos(vector_bra[i], ket[i][0]))

    return resultado


def valor_esperado(omega, ket):
    """
    :param omega:
    :param ket:
    :return numero real:
    """
    producto = productoMatricesComp(omega, ket)
    bra_producto = bra(producto)

    return redondearComplejo(productoEscalar(bra_producto, ket))


def operadorDelta(ket, omega):
    """
    :param ket:
    :param omega:
    :return: matriz:
    """
    identidad = [[(0, 0) for j in range(len(omega[0]))] for i in range(len(omega))]
    for i in range(len(omega)):
        for j in range(len(omega[0])):
            if i == j:
                identidad[i][j] = (1, 0)

    valor_e = valor_esperado(omega, ket)
    m_identidadOmega = multiplicacionEscalar_M(valor_e, identidad)
    delta = adicionMatrices(omega, inversoAditivo_M(m_identidadOmega))

    return delta


def varianza(ket, omega):
    """
    :param ket:
    :param omega:
    :return numero real:
    """
    delta = operadorDelta(ket, omega)
    producto = productoMatricesComp(delta, delta)
    return valor_esperado(producto, ket)


def valoresPropios(matriz):
    """
    :param matriz:
    :return lista de valores propios:
    :return lista de vectores propios:
    """
    valoresNumpy, vectoresNumpy = np.linalg.eig(convertirMatrizNP(matriz))
    valores, vectores = [], []

    for i in range(len(valoresNumpy)):
        valores.append((valoresNumpy[i].real, valoresNumpy[i].imag))

    for i in range(len(vectoresNumpy)):
        vectores.append([])
        for j in range(len(vectoresNumpy[0])):
            vectores[i].append((vectoresNumpy[i][j].real, vectoresNumpy[i][j].imag))

    return valores, vectores


def normalizarKet(ket):
    """
    :param ket:
    :return ket normalizado:
    """
    vKet = []
    ketNormalizado = [[None for j in range(1)] for i in range(len(ket))]
    for i in range(len(ket)):
        vKet.append(ket[i][0])

    normaKet = normaVector(vKet)
    for i in range(len(ketNormalizado)):
        valores = []
        for k in range(2):
            valores.append(ket[i][0][k] / normaKet)
        ketNormalizado[i][0] = valores[0], valores[1]

    return ketNormalizado


def probabilidadVectoresPropios(ket, vectorP):
    """
    :param ket:
    :param vectorP:
    :return probabilidad:
    """
    braKet = bra(ket)
    return moduloComplejos(productoEscalar(braKet, transpuesta(vectorP))) ** 2


def imprimir(matriz):
    """
    :param matriz:
    :return None:
    """
    for i in matriz:
        print(i)


def main():
    """
    Omega: matriz
    Ket: vector columna
    """
    """
    omega = [[(1, 0), (0, -1)], [(0, 1), (2, 0)]]
    ket = [[(math.sqrt(2)/2, 0)], [(0, math.sqrt(2)/2)]]
    print(varianza(ket, omega))
    --------------------------
    matriz = [[(-1, 0), (0, -1)],
              [(0, 1), (1, 0)]]

    print(valoresPropios(matriz))
    --------------------------
    ket = [[(2, -3)],
           [(1, 2)]]

    print(normalizarKet(ket))
    --------------------------
    matriz = [[(-1, 0), (0, -1)],
              [(0, 1), (1, 0)]]

    ket = [[(-1, 0)],
           [(-1, -1)]]

    valoresP, vectoresP = valoresPropios(matriz)
    print(probabilidadVectoresPropios(ket, vectoresP[1]))
    """


if __name__ == '__main__':
    main()
