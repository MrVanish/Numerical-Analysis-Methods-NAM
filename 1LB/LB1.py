from Constants import *
from numpy import linalg
from func import print_m

def solution(A):
    n = len(A)

    for k in range(1, n):
        for j in range(k, n):
            m = A[j][k-1] / A[k-1][k-1]
            for i in range(n+1):
                A[j][i] = A[j][i] - m*A[k-1][i]
        print(f"----------{k:} step ---------------")
        print_m(A)
        print("-------------------------------")
        print()
    x = [0]*n
    for i in reversed(range(n)):
        x[i] = A[i][n] / A[i][i]
        for c in reversed(range(i + 1, n)):
            x[i] -= A[i][c]*x[c]/A[i][i]
    return x


A = testMatrix
a = [row[:-1] for row in A]
b = [row[-1] for row in A]
A_1 = linalg.inv(a)
print("Your matrix:")
print_m(testMatrix)
print()
print(solution(testMatrix))
print()
print("pow(Matrix,-1)")
print_m(A_1)
print()
A_max = max([sum(row)] for row in A)
Inv_max = max([sum(row)] for row in A_1)
B_max = max(b)
d_x = Inv_max[0] * 0.001
_d_x = A_max[0] * Inv_max[0] * 0.001 / B_max
print("absolute and relative error")
print(d_x)
print(_d_x)
