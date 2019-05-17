from Constants import *
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import matplotlib.patches as mpatches

print(X)
print(Y)

def printLP(x, y):
    for x_n in x:
        lStrUp = ""
        lStrDown = ""
        _x = x.copy()
        _i = _x.pop(x.index(x_n))
        for i in _x:
            lStrDown += f"({_i:.3f} - {i:.3f})"
        for i in _x:
            lStrUp += f"( x - {i:.3f})"
        print("     "+lStrUp)
        print("p[{0}]=".format(x.index(x_n)) + ('-' * max([len(lStrUp), len(lStrDown)])))
        print("     " + lStrDown + "\n")

#########################################################################################################################
print("Базисные многочлены \n============================================================================================"
      +"=================")
printLP(X, Y)
#########################################################################################################################
print('Многочлен Лагранжа')
def getL_N_X(n,x):
    s = 0
    for i in range(n):
        p_i = Y[i]
        _x = X.copy()
        _i = _x.pop(i)
        for j in _x:
            p_i *= (x - j)/(_i - j)
        s += p_i
    return s
#########################################################################################################################
print(f"L[4]({X[0]+X[1]})={getL_N_X(len(X),X[0]+X[1])}")
print("=============================================================================================================")
#########################################################################################################################
#plt.figure('lagrange')
lagrangeGR = (np.linspace(X[0], X[-1], 100), [getL_N_X(len(X),x) for x in np.linspace(X[0], X[-1], 100)])
plt.plot(lagrangeGR[0], lagrangeGR[1])
plt.plot(X, [getL_N_X(len(X), x) for x in X], 'ro')
plt.axhline(0, color='black')
plt.axvline(0,color='black')
plt.grid(True, which='both')
plt.savefig('lagrange.png')
#########################################################################################################################


def table_1():
    head = ['xk', 'yk']
    for i in range(len(Y)-1):
        head.append(f"∆{i+1}yk")
    cl = []
    for i in range(len(X)):
        cl.append([X[i], Y[i]])
    y = Y.copy()
    for i in range(len(Y)-1):
        y_new = []
        for j in range(len(y)-1):
            y_new.append(y[j+1]-y[j])
        y = y_new.copy()
        for k in range(len(y_new)):
            cl[k].append(y_new.pop(0))
    print(tabulate(cl, head, tablefmt="grid", stralign='center'))
    return cl

print("Таблица конечных разностей")
t_1 = table_1()
print("=============================================================================================================")
#########################################################################################################################
def table_2():
    head = ['xk', 'yk']
    cl = []
    for i in range(len(Y)-1):
        head.append(f"{i+1}-го порядка")
    for i in range(len(X)):
        cl.append([X[i], Y[i]])
    y = Y.copy()
    for i in range(1, len(Y)):
        y_new = []
        for j in range(len(y)-1):
            y_new.append((y[j+1]-y[j])/(X[j+i]-X[j]))
        y = y_new.copy()
        for k in range(len(y_new)):
            cl[k].append(y_new.pop(0))
    print(tabulate(cl, head, tablefmt="grid", stralign='center'))
    return cl


print("Таблица разделенных разностей")
t_2 = table_2()
print("=============================================================================================================")
#########################################################################################################################
print("Полином Ньютона")


def NP():
    print(f"N[0] = {Y[0]}")
    for i in range(1, len(X)):
        str_ = f"N[{i}] = {t_2[0][i+1]} *"
        for j in range(i):
            str_ +=f"(x - {X[j]})"
        print(str_)


NP()


def N_N_X(n, x):
    s = Y[0]
    for i in range(1, n):
        p_i = t_2[0][i + 1]
        for j in range(i):
            p_i *= (x-X[j])
        s += p_i
    return s


print(f"N[4]({X[0]+X[1]})={N_N_X(len(X), X[0]+X[1])}")

########################################################################################################################
#plt.figure('newton')
newtonGR = (np.linspace(X[0], X[-1], 100), [N_N_X(len(X), x) for x in np.linspace(X[0], X[-1], 100)])
plt.plot(newtonGR[0], newtonGR[1])
plt.plot(X, [N_N_X(len(X), x) for x in X], 'ro')
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.grid(True, which='both')
plt.savefig('newton.png')
print("=============================================================================================================")

########################################################################################################################
print("Кусочно-линейный сплайн.")


def solve_liner(n, x_values, y_values):
    a = np.array([[x_values[n-1], 1],
                  [x_values[n], 1]])
    b = np.array([y_values[n-1], y_values[n]])
    an, bn = np.linalg.solve(a, b)
    return an, bn


def getsolve_liner(x, y):
    coefs = [solve_liner(k_i, x, y) for k_i in range(1, len(X))]
    for i, c in enumerate(coefs):
        print(f"a{i+1} = {c[0]} b{i+1} = {c[1]}")
    return coefs


coef_liner = getsolve_liner(X, Y)


def get_linear_spline(x_values, coef):
    def inner_fun(x):
        i = 0
        while i < len(x_values):
            if x_values[i] <= x and x <= x_values[i+1]:
                return coef[i][0]*x + coef[i][1]
            i += 1
    return inner_fun


linear_spline = get_linear_spline(X, coef_liner)
print(f"liner({X[0]+X[1]}) = {linear_spline(X[0]+X[1])}")

print("=============================================================================================================")
########################################################################################################################
#plt.figure('liner')
linerGR = (np.linspace(X[0], X[-1], 100), [linear_spline(x) for x in np.linspace(X[0], X[-1], 100)])
plt.plot(linerGR[0], linerGR[1])
plt.plot(X, [linear_spline(x) for x in X], 'ro')
plt.axhline(0, color='black')
plt.axvline(0,color='black')
plt.grid(True, which='both')
plt.savefig('liner.png')
########################################################################################################################
print("Кусочно-квадратичный сплайн.")

def solve_quadratic(k, x_values, y_values):
    try:
        a = np.array([[x_values[2*k - 2] ** 2, x_values[2*k-2], 1],
                      [x_values[2*k - 1] ** 2, x_values[2*k-1], 1],
                      [x_values[2*k] ** 2, x_values[2*k], 1]])
        b = np.array([y_values[2*k-2],
                      y_values[2*k-1],
                      y_values[2*k]])
    except:
        a = np.array([[x_values[2*k - 3] ** 2, x_values[2*k-3], 1],
                      [x_values[2*k - 2] ** 2, x_values[2*k-2], 1],
                      [x_values[2*k - 1] ** 2, x_values[2*k-1], 1]])
        b = np.array([y_values[2*k-3],
                      y_values[2*k-2],
                      y_values[2*k-1]])
    ak, bk, ck = np.linalg.solve(a, b)
    return ak, bk, ck, k


def getsolve_quadratic(x, y):
    coefs = [solve_quadratic(k_i+1, x, y) for k_i in range(len(X)//2)]
    for i, c in enumerate(coefs):
        print(f"a{i+1} = {c[0]} b{i+1} = {c[1]} c{i+1} = {c[2]}")
    return coefs


def get_quadratic_spline(x_values, coefs):
    def inner_fun(x):
        i = 0
        try:
            while i < len(x_values):
                if x < x_values[i] or (x_values[i] <= x and x <= x_values[i+2]):
                    return coefs[(i+1)//2][0] * x ** 2 + coefs[(i+1)//2][1] * x + coefs[(i+1)//2][2]
                i += 1
        except:
            print()
    return inner_fun

coef_quadratic = getsolve_quadratic(X, Y)
quadratic_spline = get_quadratic_spline(X, coef_quadratic)
print(f"quadratic({X[0]+X[1]}) = {quadratic_spline(X[0]+X[1])}")
print("=============================================================================================================")
########################################################################################################################
#plt.figure('quadratic')
quadraticGR = (np.linspace(X[0], X[-1], 100), [quadratic_spline(x) for x in np.linspace(X[0], X[-1], 100)])
plt.plot(quadraticGR[0], quadraticGR[1])
plt.plot(X, [quadratic_spline(x) for x in X], 'ro')
plt.axhline(0, color='black')
plt.axvline(0,color='black')
plt.grid(True, which='both')
plt.savefig('quadratic.png')
########################################################################################################################
def cub_int(n, x, y):
    a = y.copy()
    b = np.zeros(n)
    c = np.zeros(n)
    d = np.zeros(n)
    h = np.zeros(n)
    l = np.zeros(n)
    de = np.zeros(n)
    la = np.zeros(n)
    for i in range(1, n):
        h[i] = x[i] - x[i - 1]
        l[i] = (y[i] - y[i - 1]) / h[i]
    c[0] = 0
    c[n - 1] = 0
    for i in range(2, n):
        de[i - 1] = h[i] / (2 * h[i - 1] + 2 * h[i] + h[i - 1] * de[i - 2])
        la[i - 1] = (2 * l[i] - 3 * l[i - 1] - h[i - 1] * la[i - 2]) / (2 * h[i - 1] + 2 * h[i] + h[i - 1] * de[i - 2])
    for i in range(n - 1, 1, -1):
        c[i - 1] = de[i - 1] * c[i] + la[i - 1]
    for i in range(0, n):
        b[i] = l[i] + 2 * c[i] * h[i] / 3 + h[i] * c[i - 1] / 3
        d[i] = (c[i] - c[i - 1]) / 3 * h[i]
    return a, b, c, d


def get_cub_int(x, a, b, c, d, x0):
    for i in range(1, 5):
        if x[i - 1] <= x0 <= x[i]:
            return a[i] + b[i] * (x0 - x[i]) + c[i] * (x0 - x[i]) ** 2 + d[i] * (x0 - x[i])
    return 0

print()

a3, b3, c3, d3 = cub_int(len(X), X, Y)
print(a3,b3,c3,d3)
########################################################################################################################
#plt.figure('cubic')
x_ = np.linspace(X[0], X[len(X) - 1])
y_ = np.array([get_cub_int(X, a3, b3, c3, d3, x_[i]) for i in range(len(x_))])
plt.plot(x_, y_)
plt.plot(X, Y, 'ro')
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.grid(True, which='both')
#plt.savefig('cubic.png')
plt.plot(X, Y, 'ro')
dots_legend = mpatches.Patch(color="red", label='Dots')
lagrange_legend = mpatches.Patch(color="blue", label='Lagrange polynom')
newton_legend = mpatches.Patch(color="orange", label='Newton polynom')
linear_legend = mpatches.Patch(color="green", label='Linear splines')
qudratic_legend = mpatches.Patch(color="red", label='Quadratic splines')
cubic_legend = mpatches.Patch(color="purple", label='Cubic splines')
plt.legend(handles=[lagrange_legend, newton_legend, linear_legend, qudratic_legend, cubic_legend, dots_legend], loc='upper left');
plt.savefig('all.png')
#plt.show()
########################################################################################################################