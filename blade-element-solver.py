import numpy as np
import math

# propeller geometry variables
chord_length = 0.10
pitch = 0.0
propeller_diameter = 1.6
radius = propeller_diameter / 2
RPM = 2200
n = RPM/60.0  # revs per second
tip_angle = 25.0  # pitch angle at tip
hub_angle = 65.0  # pitch angle at hub

xt = radius  # upper bound
xs = 0.1*radius  # lower bound

rho = 1.225  # density of air (kg/m^3)

omega = n*2.0*math.pi
tonc = 0.12*chord_length
coef1 = (tip_angle-hub_angle)/(xt-xs)
coef2 = hub_angle-coef1*xs

rstep = (xt-xs)/10
r1 = np.arange(xs, xt, rstep)

k = 0
eff = 0
V = 60
thrust = 0.0
torque = 0.0

for index, j in enumerate(r1):
    rad = r1[index]
    theta = coef1*rad+coef2+pitch
    # t2(j) = theta
    th = theta/180.0*math.pi
    sigma = 2.0*chord_length/2.0/math.pi/rad
    a = 0.1
    b = 0.01
    finished = 0
    sum = 1

    while finished == 0:
        V0 = V*(1+a)
        V2 = omega*rad*(1-b)
        phi = math.atan2(V0, V2)
        alpha = th-phi
        cl = 6.2*alpha
        cd = 0.008-0.003*cl+0.01*cl*cl
        Vlocal = math.sqrt(V0*V0+V2*V2)
        DtDr = 0.5*rho*Vlocal*Vlocal*2.0*chord_length * \
            (cl*math.cos(phi)-cd*math.sin(phi))
        DqDr = 0.5*rho*Vlocal*Vlocal*2.0*chord_length * \
            rad*(cd*math.cos(phi)+cl*math.sin(phi))
        tem1 = DtDr/(4.0*math.pi*rad*rho*V*V*(1+a))
        tem2 = DqDr/(4.0*math.pi*rad*rad*rad*rho*V*(1+a)*omega)
        anew = 0.5*(a+tem1)
        bnew = 0.5*(b+tem2)
        if (abs(anew-a) < 1.0e-5):
            if (abs(bnew-b) < 1.0e-5):
                finished = 1
        a = anew
        b = bnew
        sum = sum+1

        if (sum > 500):
            finished = 1

    # a2(j) = a
    # b2(j) = b
    thrust = thrust+DtDr*rstep
    torque = torque+DqDr*rstep

Ct = thrust/(rho*n*n*pow(propeller_diameter, 4))  # thrust coeff
Cq = torque/(rho*n*n*pow(propeller_diameter, 5))  # torque coeff
J = V/(n*propeller_diameter)


if (Ct < 0):
    eff = 0
else:
    eff = J/(2.0*math.pi)*(Ct/Cq)

print (thrust, " N")
print(torque, " Nm")
