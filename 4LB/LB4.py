import numpy as np
from Constants import *
from func import *
from math import *

def symmetrize(a,b):
    tr = np.array(list(zip(*a[:])))
    print_m(tr)
    matixA = np.dot(a, tr)
    matrixB = np.dot(b, tr)
    return matixA, matrixB


def MaxElement(Matrix):
    result = None
    for i in range(len(Matrix)):
        for j in range(i+1, len(Matrix)):
                if result is None:
                    result = (Matrix[i][j], i, j)
                else:
                    if fabs(result[0]) < fabs(Matrix[i][j]):
                        result = (Matrix[i][j], i, j)
    return result


def GetAngleRotation(Matrix):
    i = MaxElement(Matrix)[1]
    j = MaxElement(Matrix)[2]
    return atan(2*Matrix[i][j]/(Matrix[i][i]-Matrix[j][j]))/2

def GetMatixRotation(maxElement, angleRotation):
    size = len(A)
    result = [[1 if j == i else 0 for j in range(size)] for i in range(size)]
    for i in range(size):
        for j in range(size):
            if ((maxElement[1] == i and maxElement[1] == j) or (maxElement[2] == i and maxElement[2] == j)):
                result[i][j] = cos(angleRotation)
            if (maxElement[1] == i and maxElement[2] == j):
                result[i][j] = 0-sin(angleRotation)
            if (maxElement[2] == i and maxElement[1] == j):
                result[i][j] = sin(angleRotation)
    return result


def MultMatrix(matrix1, matrix2):
    temp1 = np.array(matrix1)
    temp2 = np.array(matrix2)
    return list(np.dot(temp1, temp2))


def GetEigenvectors(allUMatrix, iteration, size):
    result = allUMatrix.pop(0)
    for i in range(iteration-1):
        temp = allUMatrix.pop(0)
        result = MultMatrix(result, temp)
    return [[result[j][i] for j in range(size)] for i in range(size)]



ASymmetrized,BSymmetrized = symmetrize(A,B)
print()
print("Symmetrized matrix :")
print_m(ASymmetrized)
print()
iteration = 0
UN = []
while fabs(MaxElement(ASymmetrized)[0]) > 0.01:
    print("                ======%d======" % iteration)
    print("---------------------------------------------")
    angle = GetAngleRotation(ASymmetrized)
    print("angle on iteration %d = %f" % (iteration, angle))
    print("MaxElemen A[%d][%d]= %f" % (MaxElement(ASymmetrized)[1],MaxElement(ASymmetrized)[2],MaxElement(ASymmetrized)[0]))
    U = GetMatixRotation(MaxElement(ASymmetrized), angle)
    UN.append(U)
    print("---------------------------------------------")
    print("U")
    print_m(U)
    print("---------------------------------------------")
    U_T = np.transpose(U)
    print("U_T")
    print_m(U_T)
    print("---------------------------------------------")
    ASymmetrized = MultMatrix(MultMatrix(U_T, ASymmetrized), U)
    print_m(ASymmetrized)
    print("---------------------------------------------")
    iteration += 1
for i in range(len(A)):
    print(ASymmetrized[i][i])
print("---------------------------------------------")
UT = UN.copy()
print_m(np.dot(ASymmetrized,GetEigenvectors(UT,iteration,len(A))))
print_m(GetEigenvectors(UN,iteration,len(A)))