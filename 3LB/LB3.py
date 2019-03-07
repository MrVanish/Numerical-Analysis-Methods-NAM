from func import print_m
import numpy as np
import math
from Constants import *


def symmetrize(a,b):
    tr = np.array(list(zip(*a[:])))
    matixA = np.dot(a, tr)
    matrixB = np.dot(b, tr)
    return matixA, matrixB


def U(a, b):
    u = np.array([[0.0] * len(row) for row in a])
    for i in range(len(u)):
        u[i, i] = math.sqrt(a[i, i] - sum([u[j, i] ** 2 for j in range(len(u))]))
        print("U[%d][%d] = %f" % (i, i, u[i, i]))
        for j in range(i + 1, len(u[i])):
            u[i, j] = (a[i, j] - sum([u[k, i] * u[k, j] for k in range(len(u))])) / u[i, i]
            print("U[%d][%d] = %f" % (i, j, u[i, j]))
    u_t = np.array(list((zip(*u))))
    return u, u_t


def lower(m, b):
    matrix = [[item for item in row[:]] for row in m[:]]
    for i in range(len(matrix)):
        matrix[i].append(b[i])
    answer = []
    for i in range(len(matrix)):
        var = matrix[i][-1]
        templist = [x for x in matrix[i] if x]
        if var == 0:
            templist.append(0)
        templist.reverse()
        sub = templist[0]
        for j in range(0, len(answer)):
            sub -= templist[j + 2] * answer[len(answer)-j-1]
        answer.append(sub/templist[1])
    return answer


def upper(m, b):
    matrix = [[item for item in row[:]] for row in m[:]]
    for i in range(len(matrix)):
        matrix[i].append(b[i])
    answer = []
    matrix.reverse()
    for i in range(len(matrix)):
        var = matrix[i][-1]
        templist = [x for x in matrix[i] if x]
        if var == 0:
            templist.append(0)
        templist.reverse()
        sub = templist[0]
        for j in range(0, len(answer)):
            sub -= templist[j + 1] * answer[j]
        answer.append(sub/templist[len(templist)-1])
    answer.reverse()
    return answer


def invert(u, u_t):
    res = []
    for n in range(len(u)):
        ort = [0 for i in range(len(u) - 1)]
        ort.insert(n, 1)
        y = lower(u_t, ort)
        x = upper(u, y)
        res.append(x)
    res = np.array(list(zip(*res[:])))
    return np.array(res)


if np.linalg.det(A) == 0:
    raise ValueError(" det A == 0 ")
ASymmetrized,BSymmetrized = symmetrize(A,B)
print("Symmetrized matrix :")
print_m(ASymmetrized)
print()
print(BSymmetrized)
print("------------------")
U, U_T = U(ASymmetrized, BSymmetrized)
print()
print_m(U)
print()
print_m(U_T)
print("------------------")
print("Y-vect :")
AnswerY = lower(U_T, BSymmetrized)
print(AnswerY)
print("------------------")
print("X-vect :")
AnswerX = upper(U, AnswerY)
print(AnswerX)
print("------------------")
print("Inverted matrix :")
print_m(invert(U, U_T))

