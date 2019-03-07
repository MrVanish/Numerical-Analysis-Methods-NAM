from Constants import *
import numpy as np
from func import print_m


def transform(m, e=0.001):
    m = [[e / row[i] for e in row]for i, row in enumerate(m)]
    a = [[-e for e in row[:-1]] for row in m]
    b = [row[-1] for row in m]
    a = [row[:i] + [0] + row[i+1:]for i, row in enumerate(a)]
    print_m(a)
    return a, b


def iterarion(a, b, acc=0.001):
    px = cx = b[:]
    print("Смотри сюда сука0")
    print(px)
    diff = acc
    while diff >= acc:
        cx = [0.0] * len(b)
        for i, row in enumerate(a):
            for j, e in enumerate(row):
                cx[i] += e * px[j]
            cx[i] += b[i]
        diff = max([abs(e - px[i]) for i, e in enumerate(cx)])
        px = cx
        print(px, diff)
    return diff, cx


def seidel(a, b, acc=0.01):
    px = b[:]
    diff = acc
    while diff >= acc:
        temp = px[:]
        cx = [0.0] * len(b)
        for i,row in enumerate(a):
            for j,e in enumerate(row):
                cx[i] += e * px[j]
            cx[i] += b[i]
            px[i] = cx[i]
        diff = max([abs(e-temp[i]) for i, e in enumerate(cx)])
        px = cx
        print(px,diff)
    return diff, px


print("Your matrix")
print_m(testMatrix)
print("------------")

a, b = transform(testMatrix)
print("Transformed matrix :")
print_m(a)
print(np.matrix(b))
print("------------")

print("Seidel method:")
max_val, result = seidel(a, b)
print(f"Max value is {max_val}")
print("Result:")
print(result)
print("------------")

print("Iteration method:")
max_val, result = iterarion(a, b)
print(f"Max value is {max_val}")
print("Result:")
print(result)

#[0.18527050565463515, 0.3452181748974222, 0.11436415253599609, 0.04931875027699819]