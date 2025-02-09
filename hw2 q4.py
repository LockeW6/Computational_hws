# -*- coding: utf-8 -*-

""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Adapted by Lev Kaplan 2019"""

# Lagrange.py: Langrange interpolation tabulated data
import numpy as np
import matplotlib.pyplot as plt

def legendrepol(x, xin, yin, beg, finish):
    y = 0.0
    for i in range(beg, finish + 1):
        term = yin[i]
        for j in range(beg, finish + 1):
            if i != j:
                term *= (x - xin[j]) / (xin[i] - xin[j])
        y += term
    return y

n_input = 21
xin = np.linspace(-1, 1, n_input)
yin = 1.0 / (1.0 + 25 * xin**2)


x_eval = np.linspace(-1, 1, 1000)
y_true = 1.0 / (1.0 + 25 * x_eval**2)

y_global = []
for x in x_eval:
    y = legendrepol(x, xin, yin, 0, len(xin)-1) 
    y_global.append(y)
error_global = np.array(y_global) - y_true

y_cubic = []
for x in x_eval:
    i = np.searchsorted(xin, x) - 1
    i = max(0, min(i, len(xin)-2))
    if i == 0:
        beg = 0
    elif i >= len(xin) - 3:
        beg = len(xin) - 4
    else:
        beg = i - 1
    end = min(beg + 3, len(xin)-1)

    if end - beg < 3:
        beg = max(0, end - 3)
    
    y = legendrepol(x, xin, yin, beg, end)
    y_cubic.append(y)
error_cubic = np.array(y_cubic) - y_true

y_linear = []
for x in x_eval:
    i = np.searchsorted(xin, x) - 1
    i = max(0, min(i, len(xin)-2))
    y = legendrepol(x, xin, yin, i, i+1)
    y_linear.append(y)
error_linear = np.array(y_linear) - y_true

plt.figure(figsize=(12, 6))

plt.plot(x_eval, error_global, label="Global Lagrange Error", linewidth=1)
plt.plot(x_eval, error_cubic, label="Local Cubic Error", linewidth=1)
plt.plot(x_eval, error_linear, label="Linear Error", linewidth=1)

plt.title("Interpolation Error Comparison (vs 1/(1+25x²))")
plt.xlabel("x")
plt.ylabel("Error (Interpolated - True)")
plt.legend()
plt.grid(True)
plt.show()