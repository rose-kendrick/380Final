import numpy as np 
from matplotlib import pyplot as pp 

def FgS(m1,m2,r_pos,r):
    F = (-G*m1*m2*r_pos)/(r**3)
    return F

def Fej(mj,me,r,x1,x2):
    F = (-G*mj*me*(x2-x1))/(np.abs(r**3))
    return F

def dvdt(F1,F2,m2):
    dvxdt = (F1+F2)/m2
    return dvxdt

def r(dx,dy):
    r = np.sqrt((dx)**2+(dy)**2)
    return r

G = 4*np.pi**2 # Gravitational constant (AU^3/(mS*year^2))
mE = 3.0027e-6 # Mass of Earth (Solar Mass)
mS = 1 # Mass of Sun (Solar Mass)
mJ = 0.000954588*1000 # Mass of jupiter (Solar Mass)
dt = 0.001 # Time step (years)
tFinal = 12 # End time (years)
N = int((tFinal+dt)/dt) # Number of time steps
sunx = 0 # Sun x-position (AU)
suny = 0 # Sun y-position (AU)

time = np.linspace(0,tFinal+dt,N)

xE = np.zeros(N)
yE = np.zeros(N)
vEx = np.zeros(N)
vEy = np.zeros(N)

xE[0] = 1 # Earth initial x-position (AU)
yE[0] = 0 # Earth initial y-position (AU)
vEx[0] = 0 # Earth initial x-velocity (AU/year)
vEy[0] = 2*np.pi # Earth initial y-velocity (AU/year)

xJ = np.zeros(N)
yJ = np.zeros(N)
vJx = np.zeros(N)
vJy = np.zeros(N)

xJ[0] = 5.2 # Jupiter initial x-position (AU)
yJ[0] = 0 # Jupiter initial y-position (AU)
vJx[0] = 0 # Jupiter initial x-velocity (AU/year)
vJy[0] = 2*np.pi*xJ[0]/12 # Jupiter initial y-velocity (AU/year)


for i in range(N-1):

    dxES = xE[i]-sunx
    dyES = yE[i]-suny
    dxEJ = xE[i]-xJ[i]
    dyEJ = yE[i]-yJ[i]

    FgSunEx = FgS(mS,mE,dxES,r(dxES,dyES))
    FEarthJupiterx = Fej(mJ,mE,r(dxEJ,dyEJ),xE[i],xJ[i])
    aEx = dvdt(FgSunEx,FEarthJupiterx,mE)

    FgSunEy = FgS(mS,mE,dyES,r(dxES,dyES))
    FEarthJupitery = Fej(mJ,mE,r(dxEJ,dyEJ),yE[i],yJ[i])
    aEy = dvdt(FgSunEy,FEarthJupitery,mE)

    dxJS = xJ[i]-sunx
    dyJS = yJ[i]-suny
    dxJE = xJ[i]-xE[i]
    dyJE = yJ[i]-yE[i]

    FgSunJx = FgS(mS,mJ,dxJS,r(dxJS,dyJS))
    FJupiterEarthx =Fej(mJ,mE,r(dxJE,dyJE),xJ[i],xE[i])
    aJx = dvdt(FgSunJx,FJupiterEarthx,mJ)

    FgSunJy = FgS(mS,mJ,dyJS,r(dxJS,dyJS))
    FJupiterEarthy = Fej(mJ,mE,r(dxJE,dyJE),yJ[i],yE[i])
    aJy = dvdt(FgSunJy,FJupiterEarthy,mJ)

    vEx[i+1] = vEx[i]+aEx*dt
    vEy[i+1] = vEy[i]+aEy*dt

    vJx[i+1] = vJx[i]+aJx*dt
    vJy[i+1] = vJy[i]+aJy*dt

    xE[i+1] = xE[i]+vEx[i+1]*dt
    yE[i+1] = yE[i]+vEy[i+1]*dt

    xJ[i+1] = xJ[i]+vJx[i+1]*dt
    yJ[i+1] = yJ[i]+vJy[i+1]*dt
    

pp.figure()
pp.plot(xE,yE,label='Earth',color='blue')
pp.plot(xJ,yJ,label='Jupiter',color='red')
pp.scatter(sunx,suny,color='orange',label='Sun')
pp.xlabel("x (AU)")
pp.ylabel("y (AU)")
pp.title("Position")
pp.legend(loc='upper right')
pp.axis("equal")
pp.show()
# pp.savefig('pos3bdyR1000MJ.png')

pp.figure()
pp.plot(time,xE,label="xEarth(t)",color='green')
pp.plot(time,yE,label="yEarth(t)",color='blue')
pp.plot(time,xJ,label='xJupiter(t)',color='red')
pp.plot(time,yJ,label='yJupiter(t)',color='orange')
pp.xlabel('Time (years)')
pp.ylabel('Position (AU)')
pp.title('x and y-Position vs. Time')
pp.legend(loc="upper right")
pp.show()
# pp.savefig('posvt3bdyR1000MJ.png')

