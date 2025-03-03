""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    adapted by Lev Kaplan in 2019 for linear decay fit y(x) = a + b x """
 
from pylab  import*

x = range(5,120,10)    # time from 5 to 115 in steps of 10 (12 points)
Nd = len(x)   # number of data points
y = log([32,17,21,7,8,6,5,2,2,0.1,4,1])   # log of number of counts
sig = [1] * 12   # error bars all set to 1

plot(x, y, 'bo' )                                   # Plot data in blue

errorbar(x,y,sig)                                     # Plot error bars
title('Linear Least Squares Fit')                        # Plot figure
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

tau = - 1/slope
      
print('Linear Fit Final Results\n') 
print('y(x) = a + b x')                          # Desired fit
print('a = ', inter, '+/-', sqrt(sxx/delta))                  
print('b = ', slope, '+/-', sqrt(ss/delta))
print('correlation =',-sx/sqrt(sxx*ss))
print('tau =',tau, 'ns')
# red line is the fit, red dots the fits at y[i]
t = range(0,120,1)
curve  = inter + slope*t
points = inter + slope*x
plot(t, curve,'r', x, points, 'ro')
legend(('Original Data', 'Least Squares Fit', 'Fit at Data Points'), loc='best')
show()
