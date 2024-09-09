#!/usr/bin/python3
# -*- encoding: utf-8 -*-

import numpy as np # библиотека для работы с матрицами
import matplotlib.pyplot as plt # библиотека для постоения графиков

# Константы
m = 1 # масса груза, [кг]
L = 5 # длмна стержня, [м]
g0 = 9.81 # модуль ускорения свободного падения, [м/с^2]
T = 1 # модуль силы натяжения стержня, [Н]
v = 1 # модуль скорости груза маятника, [м/с]

# Система алгебро-дифференциальных уравнений
def system(t, z):
    '''
    :param t: время
    :param z: начальные условия
    :return: система
    Дифференциальные уравнения второго порядка решаются через сведения их
        к дифференциальным уравнениям первого порядка.
    '''
    x, vx, y, vy = z
    g = g0 + 0.05*np.sin(2*np.pi*t)
    T = m*g
    dxdt = vx
    dvxdt = -(x / (m * L)) * T
    dydt = vy
    dvydt = -(y / (m * L)) * T - g
    return np.array([dxdt, dvxdt, dydt, dvydt])

# Начальные условия
z0 = [3, -v/np.sqrt(2), -4, -v/np.sqrt(2)]

# Промежуток времени
t_span = (0, 100) # Собствено промежуток времени
t_eval = np.linspace(*t_span, 1000) # Разбиение
h = t_eval[1] - t_eval[0]  # Размер шага

# Реализация метода Эйлера
def euler_method(system, z0, t_eval):
    '''
    :param system: система алгебро-дифференциальных уравнений
    :param z0: начальные условия
    :param t_eval: интервал с разбиением
    :return: численное решение дифференциальных уравнений
    '''
    z = np.zeros((len(t_eval), len(z0)))
    z[0, :] = z0
    for i in range(1, len(t_eval)):
        t = t_eval[i-1]
        z[i, :] = z[i-1, :] + h * system(t, z[i-1, :])
    return z

# Решение системы
sol = euler_method(system, z0, t_eval)

# Построение графиков
plt.figure(figsize=(12, 6))

'''
Построим четыре графика на одной фигуре
'''

# x(t)
plt.subplot(2, 2, 1)
plt.plot(t_eval, sol[:, 0], label='$x(t)$')
plt.plot(t_eval, sol[:, 2], label='$y(t)$')
plt.xlabel('Время $t$')
plt.ylabel('$x(t)$ и $y(t)$')
plt.grid(True)
plt.title('График 1. Зависимость $x$ и $y$ состовляющих радиус-вектора $\\vec{r}\\left(x, y\\right)$ от времени $t$')
plt.legend()

# vx(t) and vy(t)
plt.subplot(2, 2, 2)
plt.plot(t_eval, sol[:, 1], label='$v_x(t)$')
plt.plot(t_eval, sol[:, 3], label='$v_y(t)$')
plt.xlabel('Время $t$')
plt.ylabel('Компоненты скорости $v_x$ и $v_y$')
plt.grid(True)
plt.title('График 2. Завивисомсть составляющих скорости $v_x$ и $v_y$ от времени $t$')
plt.legend()

# y(x)
plt.subplot(2, 2, 3)
plt.plot(sol[:, 0], sol[:, 2], label='$y(x)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.grid(True)
plt.title('График 3. Зависимость $y$ от $x$')
plt.legend()

# Построение $\sqrt{x^2(t) + y^2(t)}$
plt.subplot(2, 2, 4)
r = np.sqrt(sol[:, 0]**2 + sol[:, 2]**2)
plt.plot(t_eval, r, label='$\\left|\\vec{r}\\left(x,y\\right)\\right|$')
plt.xlabel('Время $t$')
plt.ylabel('$\\left|\\vec{r}\\left(x,y\\right)\\right|$')
plt.grid(True)
plt.title('График 4. Зависимость $\\left|\\vec{r}\\right|=\\sqrt{x^2(t)+y^2(t)}$ от времени $t$')
plt.legend()

plt.tight_layout()
plt.show()
