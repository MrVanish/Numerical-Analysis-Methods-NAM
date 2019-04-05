import matplotlib as mpl
import math
import matplotlib.pyplot as plt
import numpy as np
#-----------------------------------------------------------------
def func_val(x):
    return 5*math.log(x, 2.7182818284)-math.sqrt(x)

def doPlot():
    x = np.arange(0.1, 10, 0.1)
    y = [func_val(i) for i in x]
    plt.figure()
    plt.plot(x, y)
    plt.savefig('{}.{}'.format('plot', 'png'), fmt='png')

#doPlot()
#-----------------------------------------------------------------
def chords_method(func, xk, xk_inc, e=0.001):
    print("Chords method:")
    print(f'a = {xk}, b = {xk_inc}')
    f = func
    index = 1
    while abs(xk_inc - xk) > e:
        result = (xk - f(xk)*(xk_inc-xk)/(f(xk_inc)-f(xk)))
        print(f" x{index} = {xk:.3f} - ( {f(xk)*(xk_inc - xk):.3f}/{f(xk_inc) - f(xk):.3f} ) = {result:.3f}")
        xk, xk_inc = xk_inc, result
        print(f"|ε| = {abs(xk-xk_inc)}")
        index += 1
    print(f"f({result:.3f}) = {f(result):.7f}")
    print(f"x := {result}")
    print()

chords_method(func_val,1,1.5)
#-----------------------------------------------------------------

def func_val_df(x):
    return (10-math.sqrt(x))/(2*x)

def tangent_method_calc(x):
    return x - func_val(x)/func_val_df(x)

def tangent_method(func, x, e=0.001):
    print(" Tangent method:")
    f = func
    xk = x
    xk_inc = tangent_method_calc(xk)
    index = 1
    print(f" x0 = {xk}")
    print(f" x{index} = {xk:.3f} - ( {f(xk):.2f}/{func_val_df(xk):.3f} ) = {xk_inc:.3f}")
    index += 1
    print(f"|ε| = {abs(xk_inc-xk)}")
    while abs(xk_inc-xk) > e:
        xk = xk_inc
        xk_inc = tangent_method_calc(xk)
        print(f" x{index} = {xk:.3f} - ( {f(xk):.2f}/{func_val_df(xk):.3f} ) = {xk_inc:.3f}")
        print(f"|ε| = {abs(xk_inc-xk)}")
    print(" x = ", xk_inc)

tangent_method(func_val,1)
#-----------------------------------------------------------------
def func_val_1(x):
    return math.sqrt((1-0.6*(x**2))/2)

def func_val_2(x):
    return (np.arctan(x**2)-0.2)/x

def doPlot1():
    plt.figure()
    x = [-math.sqrt(1/0.6)]
    for i in np.arange(-1.2, 1.2, 0.1):
        x.append(i)
    x.append(math.sqrt(1/0.6))
    y = [func_val_1(i) for i in x]
    plt.plot(x, y)
    y =[-i for i in y]
    plt.plot(x, y)
    plt.savefig('{}.{}'.format('plot2', 'png'), fmt='png')

#doPlot1()

#-----------------------------------------------------------------
def f1_x(x, y):
    return math.tan(x*y+0.2)/x


def f2_y(x, y=None):
    return math.sqrt(1-0.8*x**2)/2


def partial_derivative(fun, i, delta=0.0000001):
    def inner_derivative(*args):
        d_args = list(args)
        d_args[i] += delta
        return (fun(*d_args)-fun(*args))/delta
    return inner_derivative


def simple_iterations_method(x, y):
    x_inc = math.tan(x*y+0.2)/x
    y_inc = math.sqrt(1-0.8*x**2)/2
    return x_inc, y_inc


def iteration_method(x_0_p, y_0_p, e=0.001):
    d_dx = abs(partial_derivative(f1_x, 0)(x_0_p, y_0_p)) + abs(partial_derivative(f2_y, 0)(x_0_p, y_0_p))
    d_dy = abs(partial_derivative(f1_x, 1)(x_0_p, y_0_p)) + abs(partial_derivative(f2_y, 1)(x_0_p, y_0_p))
    print(f"Converges check(d/dx): {d_dx} < 1 : {d_dx < 1}")
    print(f"Converges check(d/dy): {d_dy} < 1 : {d_dy < 1 }")
    print()
    if d_dx > 1 or d_dy > 1:
        return
    x_0 = x_0_p
    y_0 = y_0_p
    index = 0
    print(f"x{index} = {x_0:.4f} y{index} = {y_0:.4f}")
    while True:
        x_1, y_1 = simple_iterations_method(x_0, y_0)
        index += 1
        print(f"x{index} = {x_1:.4f} y{index} = {y_1:.4f}")
        print(f"|x{index} - x{index-1}| = {abs(x_1-x_0):.4f}")
        print(f"|y{index} - y{index-1}| = {abs(y_1-y_0):.4f}")
        if abs(x_1-x_0) < e and abs(y_1-y_0) < e:
            break
        x_0, y_0 = x_1, y_1
    print(f"x = x{index} = {x_1:.4f} \ny = y{index} = {y_1:.4f}")


def show_matrix(message, matrix):
    print(message)
    for row in matrix:
        print("".join([(" " if element >= 0 else "") + f" {element:.4f}" for element in row]))

print()


iteration_method(0.9, 0.5)
iteration_method(0.2, -0.7)


#-----------------------------------------------------------------

def f1(x, y):
    return math.tan(x*y+0.2)-x**2


def f2(x, y):
    return 0.6*x**2 + 2*y**2-1


def f_m(xv):
    x, y = xv[0, 0], xv[1, 0]
    return np.array([f1(x, y), f2(x, y)]).reshape(-1, 1)


def jacobian(xv):
    x, y = xv[0, 0], xv[1, 0]
    return np.array([[partial_derivative(f1, 0)(x, y), partial_derivative(f1, 1)(x, y)],
                     [partial_derivative(f2, 0)(x, y), partial_derivative(f2, 1)(x, y)]])


def newton_calc(xv, index):
    show_matrix(f"F(x{index}):", f_m(xv))
    show_matrix(f"J(x{index}):", jacobian(xv))
    show_matrix(f"J(x{index})^-1:", np.linalg.inv(jacobian(xv)))
    print()
    return xv - np.dot(np.linalg.inv(jacobian(xv)), f_m(xv))


def dif_greater_than_e(vect1, vect2, epsi):
    return np.any(np.absolute(np.subtract(vect1[:, 0], vect2[:, 0])) > epsi)


def newton_method(x, y, e=0.0001):
    print("Newton method:")
    vect_0 = np.array([x, y]).reshape(-1, 1)
    index = 0
    show_matrix(f"x{index}", vect_0)
    while True:
        vect_1 = newton_calc(vect_0, index)
        index += 1
        show_matrix(f"x{index} ", vect_1)
        if not dif_greater_than_e(vect_0, vect_1, e):
            break
        vect_0 = vect_1
    print(f"x = {vect_1[0, 0]:.5f}, y = {vect_1[1, 0]:.5f}")
    return vect_1

newton_method(-0.9,-0.5)
newton_method(0.2, -0.7)

def modified_newton_method(xv, j_inv, index):
    show_matrix(f"F(x{index}):", f_m(xv))
    return xv - np.dot(j_inv, f_m(xv))


def extended_newton_method(x, y, e=0.0001):
    print("Extended Newton method:")
    index = 0
    vect_0 = np.array([x, y]).reshape(-1, 1)
    show_matrix(f"x{index}", vect_0)
    show_matrix(f"J(x{index})", jacobian(vect_0))
    J_inv = np.linalg.inv(jacobian(vect_0))
    show_matrix(f"J(x{index}^-1)", J_inv)
    while True:
        vect_1 = modified_newton_method(vect_0, J_inv, index)
        index += 1
        show_matrix(f"x{index} ", vect_1)
        if not dif_greater_than_e(vect_0, vect_1, e):
            break
        vect_0 = vect_1
    print(f"x = {vect_1[0,0]:.5f}, y = {vect_1[1,0]:.5f}")
    return vect_1

print()
extended_newton_method(-0.9,-0.5)
extended_newton_method(0.2, -0.7)