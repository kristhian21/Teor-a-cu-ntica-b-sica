import math


def sumaComplejos(c1, c2):
    return c1[0] + c2[0], c1[1] + c2[1]


def restaComplejos(c1, c2):
    return c1[0] - c2[0], c1[1] - c2[1]


def divisionComplejos(c1, c2):
    return c1[0] * c2[0] + c1[1] * c2[1] + c2[0] * c1[1] - c1[0] * c2[1]


def productoComplejos(c1, c2):
    a = (c1[0] * c2[0]) + (-(c1[1] * c2[1]))
    b = (c1[0] * c2[1]) + (c1[1] * c2[0])

    return a, b


def moduloComplejos(c):
    m = math.sqrt((c[0]) ** 2 + (c[1]) ** 2)

    return round(m, 3)


def conjugadoComplejos(c):
    return c[0], -c[1]


def cartesiano_polar(c):
    p = moduloComplejos(c)
    g = faseComplejos(c)

    return round(p, 3), g


def polar_cartesiano(polar):
    a = polar[0] * math.cos(math.radians(polar[1]))
    b = polar[0] * math.sin(math.radians(polar[1]))

    return round(a, 3), round(b, 3)


def faseComplejos(c):
    fase = math.degrees(math.atan(c[1] / c[0]))

    if c[0] < 0 and c[1] < 0:
        fase = 180 + fase
        return round(fase, 3)

    elif c[1] < 0 < c[0]:
        fase = 270 + fase
        return round(fase, 3)

    elif c[0] < 0 < c[1]:
        fase = 90 + fase
        return round(fase, 3)

    elif c[0] > 0 and c[1] > 0:
        return round(fase, 3)


def sumaVectores(v1, v2):
    result = []
    longitud = len(v1)

    for c in range(longitud):
        c1 = v1[c]
        c2 = v2[c]
        result.append((sumaComplejos(c1, c2)))

    return result


def inversoAditivo(v):
    result = []
    longitud = len(v)

    for i in range(longitud):
        c1 = -1 * v[i][0]
        c2 = -1 * v[i][1]

        result.append((c1, c2))

    return result


def multiplicacionEscalar(escalar, v):
    result = []
    longitud = len(v)

    for i in range(longitud):
        producto = productoComplejos(escalar, v[i])
        result.append(producto)

    return result


def adicionMatrices(m1, m2):
    result = []
    m = len(m1)

    for i in range(m):
        suma = sumaVectores(m1[i], m2[i])
        result.append(suma)

    return result


def inversoAditivo_M(m):
    result = []
    filas = len(m)

    for i in range(filas):
        result.append(inversoAditivo(m[i]))

    return result


def multiplicacionEscalar_M(escalar, matriz):
    filas = len(matriz)
    result = []

    for i in range(filas):
        producto = multiplicacionEscalar(escalar, matriz[i])
        result.append(producto)

    return result

def redondearComplejo(complejo):
    return round(complejo[0], 3), round(complejo[1], 3)


def transpuesta(matriz):
    filas = len(matriz)
    columnas = len(matriz[0])

    if filas == columnas:
        result = [[0 for j in range(columnas)] for i in range(filas)]

    else:
        result = [[0 for j in range(filas)] for i in range(columnas)]

    for i in range(filas):
        for j in range(columnas):
            result[j][i] = matriz[i][j]

    return result


def conjugada(matriz):
    filas = len(matriz)
    columnas = len(matriz[0])

    if type(matriz[0]) is tuple:
        matriz2 = []

        for i in range(filas):
            matriz2.append(conjugadoComplejos(matriz[i]))

    else:
        matriz2 = [[0 for j in range(columnas)] for i in range(filas)]

        for i in range(filas):
            for j in range(columnas):
                matriz2[i][j] = conjugadoComplejos(matriz[i][j])

    return matriz2


def adjunta(matriz):
    result = conjugada(matriz)
    result = transpuesta(result)

    return result


def productoMatricesComp(m1, m2):
    filasM1 = len(m1)
    filasM2 = len(m2)
    columnasM1 = len(m1[0])
    columnasM2 = len(m2[0])

    assert columnasM1 == filasM2, "Las matrices no son compatibles para hacer un producto entre ellas"

    resultado = [[0 for j in range(columnasM2)] for i in range(filasM1)]

    for i in range(filasM1):
        for j in range(columnasM2):
            suma = (0, 0)
            for k in range(filasM2):
                suma = sumaComplejos(suma, productoComplejos(m1[i][k], m2[k][j]))
            resultado[i][j] = suma

    return resultado


def accionMatrizVector(matriz, vector):
    filasM = len(matriz)
    result = [[0 for j in range(1)] for i in range(filasM)]

    for i in range(filasM):
        for j in range(1):
            suma = (0, 0)
            for k in range(len(vector)):
                suma = sumaComplejos(suma, productoComplejos(matriz[i][k], vector[k]))
            result[i][j] = suma

    return result


def productoInternoVectores(v1, v2):
    v3 = conjugada(v1)

    suma = (0, 0)
    for i in range(len(v2)):
        producto = productoComplejos(v3[i], v2[i])
        suma = sumaComplejos(suma, producto)

    return suma


def productoEscalar(vector, vectorCol):
    suma = (0, 0)
    for i in range(len(vector)):
        suma = sumaComplejos(suma, productoComplejos(vector[i], vectorCol[i][0]))

    return suma


def normaVector(vector):
    result = 0
    for i in vector:
        result += i[0] ** 2 + i[1] ** 2
    return math.sqrt(result)


def distanciaVectores(v1, v2):
    result = 0
    for i in range(len(v1)):
        result += (v1[i][0] - v2[i][0]) ** 2 + (v1[i][1] - v2[i][1]) ** 2
    return math.sqrt(result)


def matrizIdentidad(long):
    mUnitaria = [[(0, 0) for j in range(long)] for i in range(long)]
    for i in range(long):
        mUnitaria[i][i] = (1, 0)
    return mUnitaria


def esUnitaria(a):
    if productoMatrices(a, adjunta(a)) == matrizIdentidad(len(a)):
        ans = True
    else:
        ans = False
    return ans


def esHermitiana(matriz):
    return matriz == adjunta(matriz)


def productoTensor(a, b):
    result = [[[[]] for j in range(len(a[0]) * len(b[0]))] for i in range(len(a) * len(b))]

    for i in range(len(a) * len(b)):
        for j in range(len(a[0]) * len(b[0])):
            x, y = i // len(b), j // len(b[0])
            res = multiplicacionEscalar_M(a[x][y], b)
            x1, y1 = i % len(b), j % len(b[0])
            result[i][j] = res[x1][y1]
    return result


def imprimir(matriz):
    for i in matriz:
        print(i)


def main():

    omega = [[(-1, 0), (0, -1)], [(0, 1), (1, 0)]]
    ket = [[(-1, 0)], [(-1, -1)]]
    imprimir(productoMatricesComp(omega, ket))



if __name__ == '__main__':
    main()
