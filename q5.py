# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 20:14:04 2025

@author: 17680
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

epsilon = 1e-8

def transformed_func(u):
    x = (u + epsilon)**(5/3)
    dx_du = (5/3) * (u + epsilon)**(2/3)
    return (5/3) * (u + epsilon)**(-1/3) * np.cos(x**2 + 5*x)

exact_integral, _ = quad(transformed_func, epsilon, 1)

def trapezoid_rule(f, A, B, N):
    h = (B - A) / (N - 1)
    integral_sum = (f(A + epsilon) + f(B)) / 2
    for i in range(1, N - 1):
        integral_sum += f(A + i * h + epsilon)
    return h * integral_sum


def simpson_rule(f, A, B, N):
    if (N - 1) % 2 == 1:
        N += 1
    h = (B - A) / (N - 1)
    integral_sum = (f(A + epsilon) + f(B)) / 3
    for i in range(1, N - 1, 2):
        integral_sum += 4 / 3 * f(A + i * h + epsilon)
    for i in range(2, N - 1, 2):
        integral_sum += 2 / 3 * f(A + i * h + epsilon)
    return h * integral_sum


N_values = []
trap_errors = []
simp_errors = []

N = 3
maxpoints = 10000
while N < maxpoints:
    print(f"Computing for N = {N}")
    N_values.append(N)


    I_trap = trapezoid_rule(transformed_func, epsilon, 1, N)
    trap_errors.append(abs(I_trap - exact_integral))


    I_simp = simpson_rule(transformed_func, epsilon, 1, N)
    simp_errors.append(abs(I_simp - exact_integral))


    N = int(N * 1.2) + 1
    if N % 2 == 0:
        N += 1


plt.loglog(N_values, trap_errors, label="Trapezoid Error", linestyle='dashed', marker='o')
plt.loglog(N_values, simp_errors, label="Simpson Error", linestyle='solid', marker='s')
plt.loglog(N_values, 0.1 * np.array(N_values) ** (-2.0), label="O(N^-2)", linestyle='dotted')
plt.loglog(N_values, 0.1 * np.array(N_values) ** (-4.0), label="O(N^-4)", linestyle='dotted')


plt.xlabel("N")
plt.ylabel("Error")
plt.title("Error Scaling for Transformed Integral")
plt.legend()
plt.grid(True)
plt.show()

print(f"Best estimate of the integral: {exact_integral:.8f}")
