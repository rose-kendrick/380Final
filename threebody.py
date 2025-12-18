import numpy as np 
from matplotlib import pyplot as pp 


def F(m1,m2,r,dx):
    F = (-G*m1*m2*(dx))/((r**3))
    return F

def dvdt(F1,F2,m):
    dvxdt = (F1+F2)/m
    return dvxdt

def r(dx,dy):
    r = np.sqrt((dx)**2+(dy)**2)
    return r


'''TIME'''
dt = 0.0001 # Time step (years)
tFinal = 40 # End time (years)
N = int((tFinal+dt)/dt) # Number of time steps
time = np.linspace(0,tFinal+dt,N)


'''EARTH initial conditions'''
mE = 3.0027e-6 # Mass of Earth (Solar Mass)

xE = np.zeros(N)
yE = np.zeros(N)
vEx = np.zeros(N)
vEy = np.zeros(N)

xE[0] = 0.9950378913690774 # Earth initial x-position (AU)
yE[0] = 0 # Earth initial y-position (AU)
vEx[0] = 0 # Earth initial x-velocity (AU/year)
vEy[0] = 2*np.pi # Earth initial y-velocity (AU/year)
pEy = 1.886652e-5 #Earth initial y-momentum (Solar Mass*AU/year) - Used to find v0y-sun

'''JUPITER initial conditions'''
mJ = 0.000954588*1000 # Mass of jupiter (Solar Mass)
xJ = np.zeros(N)
yJ = np.zeros(N)
vJx = np.zeros(N)
vJy = np.zeros(N)

xJ[0] = 5.195037891369077 # Jupiter initial x-position (AU)
yJ[0] = 0 # Jupiter initial y-position (AU)
vJx[0] = 0 # Jupiter initial x-velocity (AU/year)
vJy[0] = 2*np.pi*xJ[0]/12 # Jupiter initial y-velocity (AU/year)
pJy = 0.0025990698 #Jupiter initial y-momentum (Solar Mass*AU/year) - Used to find v0y-sun


'''SUN initial conditions'''
mS = 1 # Mass of Sun (Solar Mass)

xS = np.zeros(N)
yS = np.zeros(N)
vSx = np.zeros(N)
vSy = np.zeros(N)

xS[0] = -0.004962108630922639 # Sun initial x-position (AU)
yS[0] = 0 # Sun initial y-position (AU)
vSx[0] = 0 # Sun initial x-velocity (AU/year)
vSy[0] = -0.0026179363 # Sun initial y-velocity (AU/year)
pSy = -0.0026179363 #Sun initial y-momentum (Solar Mass*AU/year) - Used to find v0y-sun


Xcm = (xS[0]*mS+xE[0]*mE+xJ[0]*mJ)/(mS+mE+mJ)
G = 4*np.pi**2 # Gravitational constant (AU^3/(mS*year^2))

for i in range(N-1):
    '''Earth'''
    dxES = xE[i]-xS[i]
    dyES = yE[i]-yS[i]
    dxEJ = xE[i]-xJ[i]
    dyEJ = yE[i]-yJ[i]

    FSunEx = F(mE,mS,r(dxES,dyES),dxES)
    FJupiterEx = F(mE,mJ,r(dxEJ,dyEJ),dxEJ)
    aEx = dvdt(FSunEx,FJupiterEx,mE)

    FSunEy = F(mE,mS,r(dxES,dyES),dyES)
    FJupiterEy = F(mE,mJ,r(dxEJ,dyEJ),dyEJ)
    aEy = dvdt(FSunEy,FJupiterEy,mE)

    '''Jupiter'''
    dxJS = xJ[i]-xS[i]
    dyJS = yJ[i]-yS[i]
    dxJE = xJ[i]-xE[i]
    dyJE = yJ[i]-yE[i]

    FSunJx = F(mJ,mS,r(dxJS,dyJS),dxJS)
    FEarthJx =F(mJ,mE,r(dxJE,dyJE),dxJE)
    aJx = dvdt(FSunJx,FEarthJx,mJ)

    FSunJy = F(mJ,mS,r(dxJS,dyJS),dyJS)
    FEarthJy =F(mJ,mE,r(dxJE,dyJE),dyJE)
    aJy = dvdt(FSunJy,FEarthJy,mJ)

    '''Sun'''
    dxSJ = xS[i]-xJ[i]
    dySJ = yS[i]-yJ[i]
    dxSE = xS[i]-xE[i]
    dySE = yS[i]-yE[i]

    FEarthSx = F(mS,mE,r(dxSE,dySE),dxSJ)
    FJupiterSunx =F(mS,mJ,r(dxSJ,dySJ),dySE)
    aSx = dvdt(FEarthSx,FJupiterSunx,mS)

    FEarthSy = F(mS,mE,r(dxSE,dySE),dySJ)
    FJupiterSuny = F(mS,mJ,r(dxSJ,dySJ),dySE)
    aSy = dvdt(FEarthSy,FJupiterSuny,mS)

    vEx[i+1] = vEx[i]+aEx*dt
    vEy[i+1] = vEy[i]+aEy*dt

    vJx[i+1] = vJx[i]+aJx*dt
    vJy[i+1] = vJy[i]+aJy*dt

    vSx[i+1] = vSx[i]+aSx*dt
    vSy[i+1] = vSy[i]+aSy*dt

    xE[i+1] = xE[i]+vEx[i+1]*dt
    yE[i+1] = yE[i]+vEy[i+1]*dt

    xJ[i+1] = xJ[i]+vJx[i+1]*dt
    yJ[i+1] = yJ[i]+vJy[i+1]*dt

    xS[i+1] = xS[i]+vSx[i+1]*dt
    yS[i+1] = yS[i]+vSy[i+1]*dt
    

pp.figure()
pp.plot(xE,yE,label='Earth',color='blue')
pp.plot(xJ,yJ,label='Jupiter',color='red')
pp.plot(xS,yS,label='Sun',color='orange')
pp.xlabel("x (AU)")
pp.ylabel("y (AU)")
pp.title("Position")
pp.legend(loc='lower right')
pp.axis("equal")
# pp.show()
pp.savefig('FULLpos3bdy1000MJdt0.0001.png')

pp.figure()
pp.plot(time,xE,label="xEarth(t)",color='green')
pp.plot(time,yE,label="yEarth(t)",color='blue')
pp.plot(time,xJ,label='xJupiter(t)',color='red')
pp.plot(time,yJ,label='yJupiter(t)',color='orange')
pp.plot(time,xS,label='xSun(t)',color='purple')
pp.plot(time,yS,label='ySun(t)',color='yellow')
pp.xlabel('Time (years)')
pp.ylabel('Position (AU)')
pp.title('x and y-Position vs. Time')
pp.legend(loc="lower left")
# pp.show()
pp.savefig('FULLposvt3bdy1000MJdt0.0001.png')

