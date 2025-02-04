# -*- coding: utf-8 -*-
""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Original code: TrapMethods.py
    
    2019: Extended by Lev Kaplan to include Simpson integration and loop over number of points"""

# integ.py: trapezoid and Simpson integration, a<x<b, N pts, N-1 intervals 

import numpy as np # import numpy and matplotlib
from pylab import * 
import matplotlib.pyplot as plt

def func(x):          # function to be integrated
    return np.exp(-x)

#ef rapezoid(A,B,N):   #integrate from A to B using N points
#   h = (B - A)/(N - 1)                     # step size 
#   sum = (func(A)+func(B))/2               # (1st + last)/2
#   for i in range(1, N-1):        # i goes from 1 to (N-1)-1
#      sum += func(A+i*h)
#      #sum=float32(sum)            # to simulate single-precision (32 bit) calculation
#   return h*sum  

def simpson(A,B,N):
    if ((N-1)%2==1):     #  if number of intervals odd
        print("Simpson's rule requires even number of intervals")
        return 0
    h = (B - A)/(N - 1)                     # step size 
    sum = (func(A)+func(B))/3               # (1st + last)/3
    for i in range(1, N-1,2):        # i loops over odd integers  from 1 to (N-1)-1
       sum += 4/3*func(A+i*h)
       #sum=float32(sum)
    for i in range(2, N-1,2):        # i loops over even integers starting with 2
       sum += 2/3*func(A+i*h)
       #sum=float32(sum)
    return h*sum  

def romberg_extrapolation(A, B, N):
    I_N = simpson(A, B, N)
    I_2N = simpson(A, B, 2 * N - 1)
    I_R = (4 * I_2N - I_N) / 3
    return I_R
              
A = 0.0
B = 1.0

maxpoints = 10000

Nvalues = []
simpson_error = []
romberg_error = []

exact=1-np.exp(-1)

N=3
while N<maxpoints:    # loop over number of points
   print(f"Computing for N = {N}...")
   Nvalues.append(N)

   I_simpson = simpson(A, B, N)
   simpson_error.append(abs(I_simpson - exact))

   I_romberg = romberg_extrapolation(A, B, N)
   romberg_error.append(abs(I_romberg - exact))
   
   N=int(N*1.2)+1    # N grows roughly by 1.1 factor each time
   if N%2 == 0:
       N=N+1    # make sure N is odd
        
plt.loglog(Nvalues, simpson_error, label="Simpson Error", linestyle='dashed', marker='o')
plt.loglog(Nvalues, romberg_error, label="Romberg Extrapolation Error", linestyle='solid', marker='s')
#plt.loglog(Nvalues, 0.1 * np.array(Nvalues) ** (-2.0), label="O(N^-2) (Expected Scaling")
plt.loglog(Nvalues, 0.1 * np.array(Nvalues) ** (-4.0), label="O(N^-4)", linestyle='dotted')
plt.loglog(Nvalues, 0.1 * np.array(Nvalues) ** (-6.0), label="O(N^-6)", linestyle='dotted')
        
#loglog(Nvalues,traperror,label="Trapezoid error")   # log log plot of error in trapezoid method
#loglog(Nvalues,0.1*array(Nvalues)**(-2.0),label="0.1/N^2")   # plot 0.1/N^2 for comparison 
#ylim([1e-8,1])   # set range of y values
plt.xlabel("N")
plt.ylabel("Error")
plt.title("Simpson vs. Romberg Extrapolation Error Scaling")
legend(loc="upper right")
show()    #show the graph
