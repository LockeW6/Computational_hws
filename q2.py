# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 15:24:53 2025

@author: 17680
"""

from pylab  import*
x = range(5,120,10)    # time from 5 to 115 in steps of 10 (12 points)
Nd = len(x)   # number of data points
#y = log([32,17,21,7,8,6,5,2,2,0.1,4,1])   # log of number of counts
#sig = [1] * 12   # error bars all set to 1
n_counts = [32, 17, 21, 7, 8, 6, 5, 2, 2, 0.1, 4, 1]
y = log(n_counts)
sig = [1.0 / sqrt(n) if n > 0 else 1.0 for n in n_counts]

plot(x, y, 'bo' )                                   # Plot data in blue

errorbar(x,y,sig)                                     # Plot error bars
title('Linear Least Squares Fit with Realistic Uncertainties')                        # Plot figure
xlabel( 't(ns)' )                                            # Label axes
ylabel( 'log(dN/dt )' )
grid(True)                                               # plot grid
xlim(0,120)                                              # x range for plot

ss = sx = sxx = sy = sxy = 0   # initialize various sums

for i in range(0, Nd):         # compute various sums over data points                              
    sig2 = sig[i] * sig[i]
    ss += 1. / sig2;    sx += x[i]/sig2;        sy += y[i]/sig2
    sxx += x[i] * x[i]/sig2;    sxy += x[i]*y[i]/sig2;
         
delta = ss*sxx-sx*sx;
slope = (ss*sxy-sx*sy) / delta    #calculate best fit slope
inter = (sxx*sy-sx*sxy) / delta   # calculate best fit intercept

slope_err = sqrt(ss / delta)
inter_err = sqrt(sxx / delta)
tau = - 1/slope
tau_err = slope_err / (slope ** 2)

chi2 = 0.0
for i in range(Nd):
    fit_y = inter + slope * x[i]
    chi2 += ((y[i] - fit_y)**2) / (sig[i]**2)
dof = Nd - 2
      
print('Linear Fit Final Results\n') 
print('y(x) = a + b x')                          # Desired fit
print('a = ', inter, '+/-', sqrt(sxx/delta))                  
print('b = ', slope, '+/-', sqrt(ss/delta))
print('correlation =',-sx/sqrt(sxx*ss))
print('Lifetime tau =', tau, '+/-', tau_err, 'ns')
print('chi^2 =', chi2, 'for', dof, 'degrees of freedom')
# red line is the fit, red dots the fits at y[i]
t = range(0,120,1)
curve  = inter + slope*t
points = inter + slope*x
plot(t, curve,'r', x, points, 'ro')
plot(x, points, 'ro', label='Fit at Data Points')
legend(('Original Data', 'Best Fit', 'Fit at Data Points'), loc='best')
show()