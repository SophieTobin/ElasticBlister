#This code was used to create figure 5.

import csv
import matplotlib.pyplot as plt
import math

directory = './'
files = [['de-2.txt',2],['de-3.txt',3],['de-4.txt',4],['de-4pt5.txt',4.5],['de-5.txt',5],['de-5pt5.txt',5.5]]

plt.rcParams['text.usetex'] = True

def triangle(x0,y0,xchange,ychange,scale):
    x = [x0, x0, x0, x0]
    y = [y0, y0, y0, y0]
    x[1] = x[0]*10**(xchange/scale)
    x[2] = x[1]
    y[2] = y[1]*10**(ychange/scale)
    plt.plot(x,y,color='k')
    return x, y

first = True
for fname in files:
    f = open(directory+fname[0])
    exponent = fname[1]
    csvfile = csv.reader(f,quoting=csv.QUOTE_NONNUMERIC)
    t = []
    x = []
    for row in csvfile:
        t.append(math.exp(row[0]))
        x.append(row[3])
    if first:
        tgravity = t
        xgravity = [(27*ti/4)**(1.0/3.0) for ti in t]
        first = False
    plt.loglog(t,x,label=r'$\delta = 10^{-'+str(exponent)+'}$')
plt.xlabel(r'Time $t$')
plt.xlim(left=1)
plt.ylabel(r'Front Position $x$')
plt.legend()
x, y = triangle(14.1426*10**(-0.25),2.08689*10**0.25,-17,-2,22)
plt.annotate(r'$17$',((x[0]*x[1])**0.5/10**0.15,y[0]*10**0.05))
plt.annotate(r'$2$',(x[1]*10**(-0.1),(y[1]*y[2])**0.5/10**0.025))
x, y = triangle(3126.79*10**0.1,5.043,1,1,2)
plt.annotate(r'$1$',((x[0]*x[1])**0.5/10**0.05,y[0]*10**(-0.1)))
plt.annotate(r'$1$',(x[1]*10**(0.1),(y[1]*y[2])**0.5/10**0.025))
x, y = triangle(177828,88.0125,3,1,5)
plt.annotate(r'$3$',((x[0]*x[1])**0.5/10**0.05,y[0]*10**(-0.1)))
plt.annotate(r'$1$',(x[1]*10**(0.05),(y[1]*y[2])**0.5/10**0.025))
plt.xlim([1.0, 10**6.0])
plt.show()
