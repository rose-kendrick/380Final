# just the earth and the sun
import numpy as np 
from matplotlib import pyplot as pp 

def Fg(m1,m2,r_pos,r):
    F = (-G*m1*m2*r_pos)/(r**3)
    return F

def dvdt(Fgx,m2):
    dvxdt = Fgx/m2
    return dvxdt

def r(x_i,y_i):
    r = np.sqrt(x_i**2+y_i**2)
    return r

G = 4*np.pi**2 # Gravitational constant (AU^3/(mS*year^2))
mE = 3.0027e-6 #m2 Mass of Earth (Solar Mass)
mS = 1 #m1 Mass of Sun (Solar Mass)
dt = 0.002 #Time step (years)
tFinal = 5 #End time (years)
N = int((tFinal+dt)/dt) #Number of time steps
sunx = 0 #Sun x-position (AU)
suny = 0 #Sun y-position (AU)

time = np.linspace(0,tFinal+dt,N)
xE = np.zeros(N)
yE = np.zeros(N)
vEx = np.zeros(N)
vEy = np.zeros(N)
xE[0] = 1 #Earth initial x-position (AU)
yE[0] = 0 #Earth initial y-position (AU)
vEx[0] = 0 #Earth initial x-velocity (AU/year)
vEy[0] = 2*np.pi #Earth initial y-velocity (AU/year)


for i in range(N-1):

    ax = dvdt(Fg(mS,mE,xE[i],r(xE[i],yE[i])),mE)
    ay = dvdt(Fg(mS,mE,yE[i],r(xE[i],yE[i])),mE)

    vEx[i+1] = vEx[i]+ax*dt
    vEy[i+1] = vEy[i]+ay*dt

    xE[i+1] = xE[i]+vEx[i+1]*dt
    yE[i+1] = yE[i]+vEy[i+1]*dt
    

pp.figure()
pp.plot(xE, yE,label='Earth')
pp.scatter(sunx,suny,color='orange',label='Sun')
pp.xlabel("x (AU)")
pp.ylabel("y (AU)")
pp.title("Position")
pp.legend(loc='upper right')
pp.axis("equal")
pp.savefig('position2body.png')

pp.figure()
pp.plot(time,xE,color='green',label="x(t)")
pp.plot(time,yE,label="y(t)")
pp.xlabel('Time (years)')
pp.ylabel('Position (AU)')
pp.title('x and y-Position vs. Time')
pp.legend(loc="upper right")
pp.savefig('posVsTime2body.png')

