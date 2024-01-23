import numpy as np
import scipy.linalg as la
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import sys
import csv

itmax = 10000000
#epsilon = 0.00001
epsilon = 0.005

def main(delta,itmax,epsilon,readfiles=['Gaussian','Gaussian']):
    pointsinwave = 10.0 #number of grid points per wavelength
    ttest = 1.0 #length of initial run to find best grid spacing
    largedx, smalldx, wavelength, x, h, dx, zeroloc, N, t, plotexponent, tplot = initialise(readfiles,delta)
    dt = 0.0000001
    xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
    h,x,dx,N,largegridindex = regrid(h,x,dx,N,wavelength,frontindex,zeroloc,largedx,smalldx)
    xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
    checkinterval = 1
    lastchecked = 0
    profilefreq = 4
    fullprofile =  profilefreq #cycles from 1 to profilefreq repeatedly, saving a full profile only when = profilefreq
    negative = False
    for it in range(0,itmax):
        wavelength = checkwavelength(h,frontindex,dx,N,delta)
        #increase size of domain if necessary
        if xback < 0.9*x[0]:
            h, x, dx, N, zeroloc, frontindex, largegridindex = doubleupslope(h,x,dx,N,zeroloc,frontindex,largegridindex)
            print('Doubled upslope extent')
        if xfront > x[N]-10.0*wavelength: #I don't think this actually ever happened. Regridding due to round off error typically increased domain size first. Kept in as in theory could run off the edge of the domain otherwise.
            h,x,dx,N,largegridindex = regrid(h,x,dx,N,wavelength,frontindex,largegridindex,largedx,smalldx)
            print('Increased downslope extent')
        #avoid changing gridsize early
        if t > ttest:
            if wavelength < pointsinwave*smalldx:
                print('resolution may be too low')
                print(t)
            elif wavelength > 2.0*pointsinwave*smalldx:
                print('wavelength is:')
                print(wavelength)
                print('smalldx is:')
                print(smalldx)
                xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
                smalldx = 2.0*smalldx
                h,x,dx,N,largegridindex = regrid(h,x,dx,N,wavelength,frontindex,largegridindex,largedx,smalldx)
                print('Decreased front resolution')
        if checkgrid(h,dx,largegridindex,frontindex):
            print('grid too small')
            xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
            h,x,dx,N,largegridindex = regrid(h,x,dx,N,wavelength,frontindex,largegridindex,largedx,smalldx)
        if (it-lastchecked) == checkinterval:
            t += dt
            dt, checkinterval, h, P = checktimestep(h,x,dt,dx,N,delta,epsilon,checkinterval)
            lastchecked = it
            if dt < 10.0**(-12.0):
                print('stepsize too small')
                print(dt)
                sys.exit()
        else:
            isnegative, h, P = timestep(h,x,dt,dx,N,delta)
            if isnegative:
                dt = dt/1.5
                print('decreased timestep due to negative h')
            else:
                t += dt
        xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
        if t > tplot:
            print(t)
            f = open(fileloc,'a')
            f.write(str(plotexponent))
            f.write(",")
            f.write(str(xback))
            f.write(",")
            f.write(str(xbackhead))
            f.write(",")
            f.write(str(xfront))
            f.write(",")
            f.write(str(np.amax(h)))
            f.write('\n')
            f.close()
            if fullprofile >= profilefreq:
                f = open(xfileloc,'a')
                g = open(hfileloc,'a')
                g2 = open(pfileloc,'a')
                f.write(str(plotexponent))
                g.write(str(plotexponent))
                g2.write(str(plotexponent))
                for i in range(0,N-3):
                    f.write(',')
                    f.write(str(x[i]))
                    g.write(',')
                    g.write(str(h[i]))
                    g2.write(',')
                    g2.write(str(P[i]))
                for i in range(N-3,N+1):
                    f.write(',')
                    f.write(str(x[i]))
                    g.write(',')
                    g.write(str(h[i]))
                f.write('\n')
                g.write('\n')
                g2.write('\n')
                f.close()
                g.close()
                g2.close()
                fullprofile = 1
            else:
                fullprofile += 1
            plotexponent += 0.1
            tplot = np.exp(plotexponent)
            V = np.dot(h,(dx[:N+1]+dx[1:])/2.0)
            if V > 1.1 or V < 0.9:
                print('conservation of volume error')
                sys.exit()
        if t > 1000000.0:
            break
    return

def regrid(h,x,dx,N,wavelength,frontindex,largegridindex,largedx,smalldx):
    #Create a grid spacing which varies linearly from largedx to smalldx between x[largegridindex] and
    #x[smallgridindex], assuming the previous grid was already in this form. smallgridindex will be chosen
    #to be the front of the current, while x[largegridindex] will be the max of 0 and x[frontindex]-0.5.
    #Interpolate h to get new values
    hinterp = interpolate.interp1d(x,h)
    newextent = x[frontindex] + 50.0*wavelength
    oldextent = x[N]
    xnew = x[0:largegridindex+1]
    hnew = h[0:largegridindex+1]
    dxnew = dx[0:largegridindex+1]
    largegridx = max(x[frontindex]-0.5,0)
    decaylength = x[frontindex]-largegridx
    i = largegridindex
    while xnew[i] < newextent:
        scale = (xnew[i]-largegridx)/decaylength
        if scale < 0.0:
            dxnew = np.append(dxnew,[largedx])
            largegridindex = i+1
        elif scale < 1.0:
            dxnew = np.append(dxnew,[largedx+scale*(smalldx-largedx)])
        else:
            dxnew = np.append(dxnew,[smalldx])
        xnew = np.append(xnew,[xnew[i]+dxnew[i+1]])
        if xnew[i+1] < oldextent:
            hnew = np.append(hnew,[hinterp(xnew[i+1])])
        else:
            hnew = np.append(hnew,[0])
        i += 1
    dxnew = np.append(dxnew,[smalldx])
    N = i
    return hnew,xnew,dxnew,N,largegridindex

def checkgrid(h,dx,largegridindex,frontindex):
    #checks whether grid is large enough to avoid error due to precision of floating point
    i = largegridindex
    problem = False
    while i <= frontindex:
        volumeerror = h[i]**3.0/dx[i]**5.0
        if volumeerror > 10.0**12.0: #0.1% of 10^15
            problem = True
            i = frontindex
        i += 1
    return problem
    
def checkwavelength(h,frontindex,dx,N,delta):
    #wavelength based on h < -10^(-15)*delta to avoid h<0 being down to round off error
    startwave = N
    eps = delta*10.0**(-15.0)
    for i in range(frontindex,N+1):
        if h[i] < -eps:
            startwave = i
            break
    endwave = N
    for i in range(startwave,N+1):
        if h[i] > -eps:
            endwave = i
            break
    wavelength = np.sum(dx[startwave+1:endwave+1])
    if wavelength == 0 or wavelength > 1.0:
        wavelength = delta**0.5
        print('using wavelength estimate')
    return wavelength

def doubleupslope(h,x,dx,N,zeroloc,frontindex,largegridindex):
    h = np.append(np.zeros((zeroloc)),h)
    adddx = dx[0]
    dx = np.append(np.zeros((zeroloc))+adddx,dx)
    N = N+zeroloc
    frontindex += zeroloc
    largegridindex += zeroloc
    x = np.append(np.zeros((zeroloc)),x)
    for i in range(zeroloc-1,-1,-1):
        x[i] = x[i+1]-dx[i+1]
    zeroloc = zeroloc*2
    return h,x,dx,N,zeroloc,frontindex,largegridindex

def findedges(h,x,N,delta):
    imax = np.argmax(h)
    for i in range(imax,0,-1):
        if h[i] < delta/10.0:
            xback = x[i]
            backindex = i
            break
    for i in range(imax,N+1):
        if h[i] < delta/10.0:
            xfront = x[i]
            frontindex = i
            break
    xbackhead = xback
    beginning = True
    for i in range(imax,backindex,-1):
        if beginning:
            if h[i-2] < h[imax] - (x[imax]-x[i-2])*delta/10.0:
                beginning = False
        if beginning:
            continue
        if h[i] < h[i-2] + (x[i]-x[i-2])*delta/10.0:
            xbackhead = x[i-1]
            break
    return xback, backindex, xfront, frontindex, xbackhead

def checktimestep(h,x,dt,dx,N,delta,epsilon,checkinterval):
    #THIS DOES NOT GIVE ACCURATE RESULTS IF THE LARGE TIMESTEP ERRONEOUSLY GIVES NEGATIVE H
    #If h after the timestep happens to be negative when the timestep is checked, the check will be inaccurate as the function 'timestep' will not update h.
    #This is unlikely to happen, and as far as I know did not happen in any of my original runs.
    #I will fix this in an update after uploading this to github, but have left it for now so that there is a record of the version with the bug.
    isnegative, hneworig, Porig = timestep(h,x,dt,dx,N,delta)
    isnegative, hnewsmall, Psmall = timestep(h,x,dt/2.0,dx,N,delta)
    isnegative, hnewsmall, Psmall = timestep(hnewsmall,x,dt/2.0,dx,N,delta)
    error = errorsize(hneworig,hnewsmall,N,delta)
    if error > 2.0*epsilon:
        hnew = hnewsmall
        P = Psmall
        dt = dt/1.5
        checkinterval = (1+checkinterval)//2
    elif error < epsilon/10.0:
        hnew = hneworig
        P = Porig
        dt = dt*1.5
        checkinterval = (1+checkinterval)//2
    else:
        hnew = hneworig
        P = Porig
        checkinterval += 1
    return dt, checkinterval, hnew, P

def errorsize(h1, h2, N,delta):
    diff = h1-h2
    #error = (np.dot(diff,diff)/N)**0.5 #You might like to calculate the error this way instead.
    frac = abs(diff/(h1+delta))
    error = np.amax(frac)
    if np.amin(h1+delta) < 0:
        print(error)
    return error

def timestep(h,x,dt,dx,N,delta):
    hnew = naivetimestep(h,h+0.0,dt,dx,N,delta)
    hmid = (hnew+h)/2.0
    hnew = naivetimestep(h,hmid,dt,dx,N,delta)
    if min(hnew+delta) < 0.0:
        isnegative = True
        hnew = h
    else:
        isnegative = False
    P = makePressure(x,(h+hnew)/2.0,dx,N)
    return isnegative, hnew, P

def naivetimestep(h,hmid,dt,dx,N,delta):
    hmid += delta
    #Create Ad, the advection coefficient at i
    Ad = hmid**2.0
    #Create D, the diffusion coefficient at i-1/2
    D = np.zeros((N+2))
    D[0] = (0.5*delta + 0.5*hmid[0])**3.0
    D[N+1] = (0.5*delta + 0.5*hmid[N])**3.0
    D[1:N+1] = ((hmid[1:N+1]+hmid[0:N])/2.0)**3.0
    A, b = makeAb(h,D,Ad,dt,dx,N,delta)
    hnew = la.solve_banded((3,3),A,b,overwrite_ab=True,overwrite_b=True)
    return hnew

def makeAb(h,D,Ad,dt,dx,N,delta):
    Flux = makeFlux(D,Ad,dt,dx,N)
    ident = np.zeros((7,N+1))
    for i in range(0,N+1):
        ident[3,i] = 1.0
    A = ident+0.5*Flux
    b = multiplyheptadiagonal(ident-0.5*Flux,h,N)
    longFlux = makeFlux(makelong(D,3.0,delta,delta),makelong(Ad,2.0,delta,delta),dt,makelong(dx+0.0,1.0,dx[0],dx[N+1]),N+6)
    #effect of prewetted layer on b
    blayer = multiplyheptadiagonal(-longFlux,np.zeros((N+7))+delta,N+6)
    b += blayer[3:N+4]
    return A, b

def makelong(v,power,delta1,delta2):
    #make a vector longer by appending three lots of delta1**power to the start and delta2**power to the end
    v = np.append(np.zeros((3))+delta1**power,v)
    v = np.append(v,np.zeros((3))+delta2**power)
    return v

def multiplyheptadiagonal(bmat,h,N):
    #performs matrix multiplication using a heptadiagonalmatrix in the condensed form required by la.solve_banded
    b = np.zeros((N+1))
    for i in range(3,N-2):
        b[i] = bmat[6,i-3]*h[i-3]+bmat[5,i-2]*h[i-2]+bmat[4,i-1]*h[i-1]+bmat[3,i]*h[i]+bmat[2,i+1]*h[i+1]+bmat[1,i+2]*h[i+2]+bmat[0,i+3]*h[i+3]
    b[0] = bmat[3,0]*h[0]+bmat[2,1]*h[1]+bmat[1,2]*h[2]+bmat[0,3]*h[3]
    b[1] = bmat[4,0]*h[0]+bmat[3,1]*h[1]+bmat[2,2]*h[2]+bmat[1,3]*h[3]+bmat[0,4]*h[4]
    b[2] = bmat[5,0]*h[0]+bmat[4,1]*h[1]+bmat[3,2]*h[2]+bmat[2,3]*h[3]+bmat[1,4]*h[4]+bmat[0,5]*h[5]
    b[N-2] = bmat[6,N-5]*h[N-5]+bmat[5,N-4]*h[N-4]+bmat[4,N-3]*h[N-3]+bmat[3,N-2]*h[N-2]+bmat[2,N-1]*h[N-1]+bmat[1,N]*h[N]
    b[N-1] = bmat[6,N-4]*h[N-4]+bmat[5,N-3]*h[N-3]+bmat[4,N-2]*h[N-2]+bmat[3,N-1]*h[N-1]+bmat[2,N]*h[N]
    b[N] = bmat[6,N-3]*h[N-3]+bmat[5,N-2]*h[N-2]+bmat[4,N-1]*h[N-1]+bmat[3,N]*h[N]
    return b

def makeFlux(D,Ad,dt,dx,N):
    Flux = np.zeros((7,N+1))
    Flux += makeBending(D,dt,dx,N)
    Flux += makeAdvection(Ad,dt,dx,N)
    return Flux

def makeBending(D,dt,dx,N):
    dxmid = np.zeros((N+1))
    dxmid[:] = (dx[:N+1]+dx[1:])/2.0 # dx at i
    Bending = np.zeros((7,N+1))
    Bendingtemp = np.zeros((7,N+1))
    #-ve 2nd derivative
    Bending[2,1:N+1] = -dt/(dx[1:N+1]*dxmid[0:N])
    Bending[3,0:N+1] = 2.0*dt/(dx[0:N+1]*dx[1:N+2])
    Bending[4,0:N] = -dt/(dx[1:N+1]*dxmid[1:N+1])
    #4th derivative
    Bendingtemp[1,2:N+1] = Bending[2,2:N+1]/(dx[1:N]*dxmid[0:N-1])
    Bendingtemp[2,1:N+1] = -2.0*Bending[2,1:N+1]/(dx[0:N]*dx[1:N+1])+Bending[3,1:N+1]/(dx[1:N+1]*dxmid[0:N])
    Bendingtemp[3,0:N+1] = Bending[2,0:N+1]/(dx[0:N+1]*dxmid[0:N+1])-2.0*Bending[3,0:N+1]/(dx[0:N+1]*dx[1:N+2])+Bending[4,0:N+1]/(dx[1:N+2]*dxmid[0:N+1])
    Bendingtemp[4,0:N] = Bending[3,0:N]/(dx[1:N+1]*dxmid[1:N+1])-2.0*Bending[4,0:N]/(dx[1:N+1]*dx[2:N+2])
    Bendingtemp[5,0:N-1] = Bending[4,0:N-1]/(dx[2:N+1]*dxmid[2:N+1])
    #6th derivative with diffusion coefficient
    Bending[0,3:N+1] = (D[1:N-1]*Bendingtemp[1,3:N+1]/dx[1:N-1])/dxmid[0:N-2]
    Bending[1,2:N+1] = (-D[0:N-1]*Bendingtemp[1,2:N+1]/dx[0:N-1]+D[1:N]*(-Bendingtemp[1,2:N+1]+Bendingtemp[2,2:N+1])/dx[1:N])/dxmid[0:N-1]
    Bending[2,1:N+1] = (D[0:N]*(Bendingtemp[1,1:N+1]-Bendingtemp[2,1:N+1])/dx[0:N]+D[1:N+1]*(-Bendingtemp[2,1:N+1]+Bendingtemp[3,1:N+1])/dx[1:N+1])/dxmid[0:N]
    Bending[3,0:N+1] = (D[0:N+1]*(Bendingtemp[2,0:N+1]-Bendingtemp[3,0:N+1])/dx[0:N+1]+D[1:N+2]*(-Bendingtemp[3,0:N+1]+Bendingtemp[4,0:N+1])/dx[1:N+2])/dxmid[0:N+1]
    Bending[4,0:N] = (D[1:N+1]*(Bendingtemp[3,0:N]-Bendingtemp[4,0:N])/dx[1:N+1]+D[2:N+2]*(-Bendingtemp[4,0:N]+Bendingtemp[5,0:N])/dx[2:N+2])/dxmid[1:N+1]
    Bending[5,0:N-1] = (D[2:N+1]*(Bendingtemp[4,0:N-1]-Bendingtemp[5,0:N-1])/dx[2:N+1]-D[3:N+2]*Bendingtemp[5,0:N-1]/dx[3:N+2])/dxmid[2:N+1]
    Bending[6,0:N-2] = (D[3:N+1]*Bendingtemp[5,0:N-2]/dx[3:N+1])/dxmid[3:N+1]
    return Bending

def makeAdvection(Ad,dt,dx,N):
    Advection = np.zeros((7,N+1))
    twodxmid = np.zeros((N+1))
    twodxmid[:] = dx[:N+1]+dx[1:] # 2*dx at i
    Advection[2,1:N+1] = dt*Ad[1:N+1]/twodxmid[0:N]
    Advection[4,0:N] = -dt*Ad[0:N]/twodxmid[1:N+1]
    return Advection

def makePressure(x,h,dx,N):
    #NOTE: Pressure excluding last 2 points at each end
    derivatives = np.zeros((N+1))
    #2nd derivative
    derivatives[1:N] = 2*((h[2:N+1]-h[1:N])/dx[2:N+1] - (h[1:N]-h[0:N-1])/dx[1:N])/(dx[1:N]+dx[2:N+1])
    #4th derivative
    derivatives[2:N-1] = 2*((derivatives[3:N]-derivatives[2:N-1])/dx[3:N] - (derivatives[2:N-1]-derivatives[1:N-2])/dx[2:N-1])/(dx[2:N-1]+dx[3:N])
    #pressure
    P = derivatives[2:N-1] - x[2:N-1]
    return P

def initialise(filelocs,delta):
    if filelocs[0] == 'Gaussian':
        #If no files specified to read in, initialise with a Gaussian.
        largedx = 1.0/128.0
        smalldx = 1.0/128.0
        N = 512
        wavelength = delta**0.5
        while smalldx > wavelength/pointsinwave:
            smalldx = smalldx/2.0
        print(smalldx)
        zeroloc = N//2
        x = np.zeros((N+1))
        h = np.zeros((N+1))
        dx = np.zeros((N+2))
        dx = dx+largedx
        for i in range(zeroloc+1,N+1):
            x[i] = x[i-1]+dx[i]
            h[i] = np.exp(-(4.0*x[i])**2)*4.0/(np.pi**0.5)
        for i in range(zeroloc-1,-1,-1):
            x[i] = x[i+1]-dx[i+1]
            h[i] = np.exp(-(4.0*x[i])**2)*4.0/(np.pi**0.5)
        h[zeroloc] = 4.0/(np.pi**0.5)
        xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
        h,x,dx,N,largegridindex = regrid(h,x,dx,N,wavelength,frontindex,zeroloc,largedx,smalldx)
        t = 0.0
        plotexponent = 0.1
        tplot = np.exp(plotexponent)
    else:
        #Otherwise, read in some old files and start wherever they left off
        xf = open(filelocs[0])
        xcsvfile = csv.reader(xf,quoting=csv.QUOTE_NONNUMERIC)
        hf = open(filelocs[1])
        hcsvfile = csv.reader(hf,quoting=csv.QUOTE_NONNUMERIC)
        for row in xcsvfile:
            plotexponent = row[0]
            x = np.array(row[1:])
            row2 = next(hcsvfile)
            h = np.array(row2[1:])
        t = np.exp(plotexponent)
        plotexponent += 0.1
        tplot = np.exp(plotexponent)
        largedx = x[1]-x[0]
        N = len(x)-1
        smalldx = x[N]-x[N-1]
        dx = x[1:]-x[:N]
        dx = np.append([largedx],dx)
        dx = np.append(dx,[smalldx])
        zeroloc = np.abs(x).argmin()
        xback, backindex, xfront, frontindex, xbackhead = findedges(h,x,N,delta)
        wavelength = checkwavelength(h,frontindex,dx,N,delta)
    return largedx, smalldx, wavelength, x, h, dx, zeroloc, N, t, plotexponent, tplot
    

fileloc = 'de-5.txt'
xfileloc = 'xde-5.txt'
#oldxfileloc = 'xde-5-old.txt'
hfileloc = 'hde-5.txt'
#oldhfileloc = 'hde-5-old.txt'
pfileloc = 'Pde-5.txt'
f = open(fileloc,'x')
f.close()
f = open(xfileloc,'x')
f.close()
f = open(hfileloc,'x')
f.close()
f = open(pfileloc,'x')
f.close()
main(0.1**5.0,itmax,epsilon)
#main(0.1**5.0,itmax,epsilon,readfiles=[oldxfileloc,oldhfileloc]) #Use this if you want to read in an old file for initial conditions.
