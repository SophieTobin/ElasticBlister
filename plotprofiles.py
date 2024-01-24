#Code used to produce figure 2

import csv
import matplotlib.pyplot as plt
import math
import matplotlib.cm as cm
import matplotlib.colors as colors

plt.rcParams['text.usetex'] = True

directory = './'
#Plots in paper used delta = 10^{-3}
#fnames = ['hde-4.txt','xde-4.txt']
#exponent = 4.0
#To plot output of blister.py without changing delta from the default 10^{-3}, use these:
fnames = ['hde-3.txt','xde-3.txt']
exponent = 3.0

f = open(directory+fnames[0])
csvfile = csv.reader(f,quoting=csv.QUOTE_NONNUMERIC)
xf = open(directory+fnames[1])
xcsvfile = csv.reader(xf,quoting=csv.QUOTE_NONNUMERIC)
cmapearly = cm.ScalarMappable(norm=colors.LogNorm(vmin=1,vmax=10.0**(exponent/2.0)),cmap='viridis')
cmapmid = cm.ScalarMappable(norm=colors.LogNorm(vmin=10.0**(exponent/2.0),vmax=10**5.0),cmap='viridis')
cmaplate = cm.ScalarMappable(norm=colors.LogNorm(vmin=10.0**(5),vmax=10.0**6),cmap='viridis')
cycle = 1
spcycle = 1
trcycle = 3
for row in csvfile:
    t = math.exp(row[0])
    h = [value + 10.0**(-exponent) for value in row[1:]]
    xlong = next(xcsvfile)
    if cycle >= 1:
        if t < 10.0**(exponent/2.0):
            if spcycle >= 2:
                plt.figure(1)
                plt.plot(xlong[1:],h,color=cmapearly.to_rgba(t))
                print(str(t)+' spreading')
                spcycle = 0
            spcycle += 1
        elif t < 10.0**5.0:
            if trcycle >= 3:
                plt.figure(2)
                plt.plot(xlong[1:],h,color=cmapmid.to_rgba(t))
                print(str(t)+' translating')
                trcycle = 0
            trcycle += 1
        else:
            plt.figure(3)
            plt.plot(xlong[1:],h,color=cmaplate.to_rgba(t))
            print(str(t)+' gravity')
        cycle = 0
    cycle += 1
plt.figure(1)
plt.colorbar(cmapearly,orientation='horizontal',label=r'Time $t$')
plt.xlabel(r'Downslope position $x$')
plt.ylabel(r'Height')
plt.ylim(bottom=0)
ax = plt.gca()
ax.set_aspect(1.0/(4.0*ax.get_data_ratio()))
plt.figure(2)
plt.colorbar(cmapmid,orientation='horizontal',label=r'Time $t$')
plt.xlabel(r'Downslope position $x$')
plt.ylabel(r'Height')
plt.ylim(bottom=0)
ax = plt.gca()
ax.set_aspect(1.0/(4.0*ax.get_data_ratio()))
plt.figure(3)
plt.colorbar(cmaplate,orientation='horizontal',label=r'Time $t$')
plt.xlabel(r'Downslope position $x$')
plt.ylabel(r'Height')
plt.ylim(bottom=0)
ax = plt.gca()
ax.set_aspect(1.0/(4.0*ax.get_data_ratio()))
plt.show()
