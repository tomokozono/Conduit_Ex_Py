import numpy as np
import matplotlib.pyplot as plt
import sys
import math

# Parameters
pch  = 1.961E+8 # Chamber pressure (Pa)
z0   = -8000.0  # Chamber depth (m)
n0   = 0.04     # Initial H2O content
T    = 1173.0   # Temperature (K)
beta = 0.4      # Crystallinity
rc   = 10.0     # Conduit radius (m)

# Constants
R     = 462.0   # Gas constant (J/kg/K)
g     = 9.8     # Gravity acceleration (m/s2)
rhol  = 2500.0  # Liquid density (kg/m3)
rhoc  = 2500.0  # Crystal density (kg/m3)
s     = 4.11E-6 # Solubility constant (Pa^{-0.5})
pa    = 1.0E+5  # Atmospheric pressure (Pa)
phicr = 0.8     # Critical porosity for fragmentation

# Constants for q-loop and p-grid
qint  = 1.0E+2  # Initial value for loop of discharge rate (kg/m2/s)
dq    = 1.0E+1  # dq for loop of discharge rate (kg/m2/s)
dpmul = 10.0**(-0.01)   # For p-grid (p = p_old*dpmul)

# Constants for Hess and Dingwell (1996)
etaa1 = -3.545
etaa2 = 0.833
etab1 = 9601.0
etab2 = -2368.0
etac1 = 195.7
etac2 = 32.25

# Constants for Costa (2005)
betast  = 0.673
alphaCo = 0.999916
deltaCo = 16.9386
gammaCo = 3.98937

# Calculated constants
rholc = (1.0-beta)*rhol + beta*rhoc #Liquid+crystals density
nl = (1.0-beta)*rhol/rholc  # Mass fraction of liquid in liquid+crystals
# Crystal effect on vicsocity by Costa (2005)
funcb = (1.0 + (beta/betast)**deltaCo) \
    /(1.0 - alphaCo *math.erf(math.sqrt(math.pi)/(2.0*alphaCo)*beta/betast \
    *(1.0 + (beta/betast)**gammaCo)))**(2.5*betast)

# Reset of list
zpoints = []    # for z
ppoints = []    # for p (pressure)
vpoints = []    # for v (velocity)
phipoints = []  # for phi (porosity)

# Gas density (kg/m3)
def rhog(p):
    return p/(R*T)

# H2O concentration
def c(p):
    if n0 > (nl*s*p**0.5):
        return s*p**0.5
    else:
        return n0/nl

# Mass fraction of gas
def n(p):
    if n0 > (nl*c(p)):
        return (n0-nl*c(p))/(1.0-nl*c(p))
    else:
        return 0.0

# Magma density (kg/m3)
def rho(p):
    return (n(p)/rhog(p)+(1.0-n(p))/rholc)**(-1.0)

# Velocity (m/s)
def v(p):
    return q/rho(p)

# Porosity
def phi(p):
    return n(p)*rho(p)/rhog(p)

# Liquid vicsocity by Hess and Dingwell (1996)
def etal(p):
    return 10.0** \
        (etaa1+etaa2*math.log(100.0*c(p)) \
        +(etab1+etab2*math.log(100.0*c(p))) \
        /(T-(etac1+etac2*math.log(100.0*c(p)))))

# Porosity effect on vicsocity by Llewellin et al. (2002)
def funcphi(p):
    return (1.0-phi(p))**(5.0/3.0)

# Magma viscosity
def eta(p):
    return etal(p)*funcb*funcphi(p)

# Friction force
def Fric(p):
    if phi(p) < phicr:
        return 8.0*eta(p)*v(p)/(rc**2.0)
#        return 8.0*eta*v(p)/(rc**2.0)
    else:
        return 0.0

# Mach number squared (i.e., M^2)
def M2(p):
    if n0 > (nl*c(p)):
        return n(p)*R*T*(q**2.0)/(p**2.0) \
               * (1.0 + 0.5*nl*c(p)*(1.0-n(p)) \
                       /(n(p)*(1.0-nl*c(p)))*(1.0-p/(rholc*R*T)))
    else:
        return 0.0

# Function for calculation (dp/dz)
def f(p):
    return -(rho(p)*g + Fric(p))/(1.0-M2(p))

### Determine discharge rate  (q) by shooting method ###
q = qint
for i1 in range(1, 1000000, 1):
    q += dq

    # Initial conditions
    z = z0  #z
    p = pch #p (pressure)
    p_old = p
    z_old = z

    # Main loop
    for i2 in range(1, 1000000, 1):

    # Get p
        p = p_old * dpmul
        dp = p - p_old

    # Get z
        z = z_old + dp/f(p)

    # Set p_old and z_old
        p_old = p
        z_old = z

    # Stop if Mach number (M) is equal to 1 (i.e., sonic flow)
        if M2(p) > 1.0:
            break

    # Stop if M = 1 at the surface (z = 0)
    if z < 0.0:
        break

if i1 == 1:
    print("Calculation Failed!! Please set qint lower than",qint,"(kg/m2/s)")
else:
    print("Discharge rate =",q,"(kg/m2/s)")

    ### Output of steady solution ###
    file = open('output.txt', 'w')

    # Initial conditions
    z = z0  #z
    p = pch #p (pressure)
    p_old = p
    z_old = z

    # Main loop
    for i3 in range(1, 1000000, 1):

        # Add variables to list
        zpoints.append(z)
        ppoints.append(p)
        vpoints.append(v(p))
        phipoints.append(phi(p))

        # Get p
        p = p_old * dpmul
        dp = p - p_old

        # Get z
        z = z_old + dp/f(p)

        # Set p_old and z_old
        p_old = p
        z_old = z

        # Stop if Mach number is equal to 1 (i.e., sonic flow)
        if M2(p) > 1.0:
            break

        # Output of results
        # File name: "output.txt"
        # Column 1: Heights (m)     Column 2: Pressure (Pa)
        # Column 3: Velocity (m/s)  Column 4: Porosity
        output = "{:.7f} {:.7f} {:.7f} {:.7f}\n".format(z, p, v(p), phi(p))
        file.write(output)

    file.close()

    # For plotting
    plt.subplot(1,3,1)
    plt.plot(ppoints, zpoints)
    plt.xlabel("Pressure (Pa)")
    plt.ylabel("Heights (m)")

    plt.subplot(1,3,2)
    plt.plot(phipoints, zpoints)
    plt.xlabel("Porosity")
    plt.yticks(color="None")

    plt.subplot(1,3,3)
    plt.plot(vpoints, zpoints)
    plt.xlabel("Velocity (m/s)")
    plt.xscale('log')
    plt.yticks(color="None")

    plt.show()
