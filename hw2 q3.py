""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Adapted by Lev Kaplan 2019"""

# Lagrange.py: Langrange interpolation tabulated data
    
from pylab import *

def legendrepol (x,beg,finish):          # poly interpolation at x 
    y = 0.                               # using input points from beg to finish
    for  i in range(beg,finish+1): 
       lambd = 1.0;
       for j in range(beg,finish+1):
           if i != j:                       #Lagrange polynom formed here
              lambd = lambd * ((x - xin[j])/(xin[i] - xin[j]))
       y += yin[i] * lambd
    return y

NMAX = 100  # max number of input points

xin = zeros(NMAX)  # each is array of length NMAX, all elements set to zero
yin = zeros(NMAX)

inputfile = open("lagrange.dat","r")  # read in the input x,y values
r = inputfile.readlines()  # read the whole file into list (one item per line)
inputfile.close()
        # input has the form: x0 y0
        #                     x1 y1
        #                     ...

m = 0
for line in r:
    #print(line)
    s = line.split() # split line and split into list of items(assume items separated by spaces)
    xin[m] = s[0] # first number in each line is the x value
    yin[m] = s[1]
    print(xin[m],yin[m])
    m+=1         # m is total number of input data points
                 # will be stored in xin[0]..xin[m-1],yin[0].yin[m-1]

xvalues=range(0,201,5)

yvalues_cubic = []
for x in xvalues:
    i = 0
    while i < m - 1 and x > xin[i + 1]:
        i += 1

    if i == 0:
        beg = 0 
    elif i >= m - 3:
        beg = m - 4
    else:
        beg = max(0, i - 1)
    
    end = min(beg + 3, m - 1)
    if end - beg < 3:
        beg = max(0, end - 3)

    y = legendrepol(x, beg, end)
    yvalues_cubic.append(y)

#numpoints = m # use all points for interpolation
#firstpoint = 0 # first point to use for interpolation
yvalues_linear = []
for x in xvalues:         # now interpolate  
    for i in range(m-1):
        if xin[i] <=x<=xin[i+1]:
            firstpoint=i
            break
    y=legendrepol(x, firstpoint, firstpoint+1)
    yvalues_linear.append(y)

plot(xin[:m], yin[:m], "o", markersize=8, label="Input Data")
plot(xvalues, yvalues_cubic, "-", linewidth=2, label="Cubic Interpolation (Local)")
plot(xvalues, yvalues_linear, "--", linewidth=2, label="Linear Interpolation")
xlabel("x")
ylabel("y")
title("Interpolation Comparison: Linear vs Cubic")
legend(loc="upper right")
grid(True)
show()
