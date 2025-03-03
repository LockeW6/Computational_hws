# -*- coding: utf-8 -*-


from pylab import *
from scipy.optimize import curve_fit


def func(x, A, B, C):
    return A / (1 + B*(x-C)**2)


data = loadtxt("lagrange.dat")
xdata = data[:, 0]
ydata = data[:, 1]

sig = sqrt(ydata)


A0 = max(ydata)
C0 = xdata[argmax(ydata)]
B0 = 1.0
p0 = [A0, B0, C0]

popt, pcov = curve_fit(func, xdata, ydata, p0=p0, sigma=sig, absolute_sigma=True)


print("Best fit parameters A, B, C =", popt)
print("Covariance matrix:\n", pcov)
print("Parameter uncertainties =", sqrt(diag(pcov)))


figure()
errorbar(xdata, ydata, sig, fmt='o', label='Data')
x_fit = linspace(min(xdata), max(xdata), 100)
y_fit = func(x_fit, *popt)
plot(x_fit, y_fit, 'r-', label='Best Fit: A/(1+B*(x-C)^2)')
xlabel('x (units)')
ylabel('y (units)')
legend(loc='best')
title('Nonlinear Fit to lagrange.dat')
grid(True)
show()
