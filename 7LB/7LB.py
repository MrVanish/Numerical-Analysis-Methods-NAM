from Constants import *
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import matplotlib.patches as mpatches
import math


Answ = integrate.quad(f, A, B)[0]
print(f'Ответ модуля scipy = {Answ}')
print('=============================================================')
########################################################################################################################
e = 0.001
M = abs(max([f2(x) for x in range(A, B+1)]))
H = math.sqrt((e * 12)/(M * abs(B-A)))
N = abs(B - A)//H
while N % 4 != 0:
    N += 1
H = abs(B - A)/N
N = int(N)
########################################################################################################################


def trapezium(f, x0, n, dx):
    return sum([f(A)/2, f(B)/2] + [f(X) for X in np.linspace(A+dx, B-dx, n-1)])*dx


print('Метод Трапеций h')
print(f"x0 = {A}, n = {N}, h = {H}")
integral1 = trapezium(f, A, N, H)
print(f"Полученное значение = {integral1}")
print(f"∆ = {abs(integral1 - Answ)}")
print()
print('=============================================================')
print('Метод Трапеций 2h')
print(f"x0 = {B}, n = {N // 2}, h = {H * 2}")
integral2 = trapezium(f, A, N // 2, H * 2)
print(f"Полученное значение = {integral2}")
print(f"∆ = {abs(integral2 - Answ)}")
print()
print('=============================================================')
########################################################################################################################


def simpson(a, b, n, dx):
    return (1.0 / 3.0) * dx * (f(a) + f(b) + 4 * sum([f(a + k * dx) for k in range(1, n, 2)]) + 2 * sum([f(a + k * dx) for k in range(2, n, 2)]))


print('Метод Симпсона h')
print(f"a = {A}, b = {B}, n = {N}, h = {H}")
res1 = simpson(A, B, N, H)
print(f'Полученое значение = {res1}')
print(f"∆ = {abs(res1 - Answ)}")
print()
print('=============================================================')
print('Метод Симпсона 2h')
print(f"a = {A}, b = {B}, n = {N // 2}, h = {H * 2}")
res2 = simpson(A, B, N // 2, H * 2)
print(f'Полученое значение = {res2}')
print(f"∆ = {abs(res2 - Answ)}")
print()
print('=============================================================')
########################################################################################################################


def dy(x, y):
    return ((y**2 * math.log(x))-y)/x


eps = 0.0001
A = 1
B = 5
H = eps**(1/4)
N = (B-A)/H
X = 1
Y = 0.5

print("Грубое значение шага = ", H)
print("Интервалов = ", N)


def runge_kut_op(x, y, h):
    f1 = dy(x, y)
    f2 = dy(x + h / 2, y + f1 * h / 2)
    f3 = dy(x + h / 2, y + f2 * h / 2)
    f4 = dy(x + h, y + h * f3)
    return y + (h / 6) * (f1 + 2 * f2 + 2 * f3 + f4)


def H_true(x, y, h):
    while True:
        y1 = runge_kut_op(x, y, h)
        y2 = runge_kut_op(x + h, y1, h)
        y2_ = runge_kut_op(x, y, 2 * h)
        if abs(y2 - y2_) / 15 < eps:
            break
        else:
            h /= 2
    return h


def runge_kut(a, b, h, n):
    x = np.linspace(a, b, n + 1)
    y = [Y]
    yk = Y
    for i in range(1, len(x)):
        yk = runge_kut_op(x[i - 1], yk, h)
        y.append(yk)
    return y


H = H_true(X, Y, H)
N = int(N)
RK_X_h = np.linspace(A, B, N + 1)
RK_X_2h = np.linspace(A, B, N//2 + 1)
RK_h = runge_kut(A, B, H, N)
RK_2h = runge_kut(A, B, H * 2, N // 2)
########################################################################################################################
plt.figure('Метод Рунге-Кутта) , метод Адамса ,метод Эйлера ')
plt.axhline(0, color='black')
plt.axvline(0,color='black')
plt.grid(True, which='both')
plt.plot(RK_X_h, RK_h, 'black')
#plt.plot(RK_X_2h, RK_2h, 'black')
########################################################################################################################


def adams(a, b, h, n):
    x = np.linspace(a, b, n + 1)
    y = [Y]
    y.append(runge_kut_op(x[0], y[0], h))

    for i in range(2, len(x)):
        yk = y[i - 1] + (h / 2) * (3 * dy(x[i - 1], y[i - 1]) - dy(x[i - 2], y[i - 2]))
        y.append(yk)
    return y


A_h = adams(A, B, H, N)
A_2h = adams(A, B, H * 2, N // 2)
########################################################################################################################
plt.plot(RK_X_h, A_h, 'blue')
#plt.plot(RK_X_2h, A_2h, 'yellow')
########################################################################################################################
print("Метод Рунге-Кутта")
head = ['x', 'yh','y2h','∆']
table = [[RK_X_2h[i], RK_h[i*2], RK_2h[i], abs(RK_h[i*2]-RK_2h[i])] for i in range(1, len(RK_X_2h))]
print(tabulate(table, head, tablefmt="grid", stralign='center'))
########################################################################################################################
print("Метод Адамса")
head = ['x', 'yh','y2h','∆']
table = [[RK_X_2h[i], A_h[i*2], A_2h[i], abs(A_h[i*2]-A_2h[i])] for i in range(1,len(RK_X_2h))]
print(tabulate(table, head, tablefmt="grid", stralign='center'))
########################################################################################################################


def eiler(a, b, h, n):
    x = np.linspace(a, b, n + 1)
    y = [Y]
    for i in range(1, len(x)):
        y.append(y[i - 1] + h * dy(x[i - 1], y[i - 1]))
    return y


E_h = eiler(A, B, H, N)
E_2h = eiler(A, B, H * 2, N // 2)
########################################################################################################################
plt.plot(RK_X_h, E_h, 'brown')
#plt.plot(RK_X_2h, E_2h, 'green')
########################################################################################################################
print("Метод Эйлера")
head = ['x', 'yh', 'y2h', '∆']
table = [[RK_X_2h[i], E_h[i*2], E_2h[i], abs(E_h[i*2]-E_2h[i])] for i in range(1, len(RK_X_2h))]
print(tabulate(table, head, tablefmt="grid", stralign='center'))
########################################################################################################################
h =[
    mpatches.Patch(color="black", label='Точное решение'),
    mpatches.Patch(color="brown", label='Эйлер'),
    #mpatches.Patch(color="green", label='Эйлер 2h'),
    mpatches.Patch(color="blue", label='Адамс'),
    #mpatches.Patch(color="yellow", label='Адамс 2h'),
    mpatches.Patch(color="red", label='Рунг-Кутт')
    #mpatches.Patch(color="black", label='Рунг-Кутт 2h')
    ]
plt.legend(handles=h, loc='best')
########################################################################################################################
d_true_y = lambda x: 1/(x + math.log(x) + 1)
true_y = [d_true_y(i) for i in RK_X_h]
plt.plot(RK_X_h, true_y, 'red')
########################################################################################################################
print("Конечная таблица")
head = ['x', 'yh-Рунге','yh','∆','yh-Адамс','yh','∆','yh-Эйлер','yh','∆']
table = []
for i in range(1, len(RK_X_h)):
    table.append([RK_X_h[i], RK_h[i], true_y[i], abs(true_y[i]-RK_h[i]), A_h[i], true_y[i], abs(true_y[i]-A_h[i]), E_h[i], true_y[i], abs(true_y[i]-E_h[i])])
print(tabulate(table, head, tablefmt="grid", stralign='center'))
########################################################################################################################
plt.show()