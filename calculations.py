import math as m
import numpy as np
import matplotlib.pyplot as pp

gravity = 9.8
#permeability of free space
mu_0 = 4*m.pi*10**(-7)

# Just some examples of realative permeability
# mu_r means realitive permeability to free space
# mu_r_ElectricalSteel = mu_0*4000
# mu_r_PureIron = mu_0*200000
# mu_r_Iron = mu_0*5000
# mu_r_CarbonSteel = mu_0*100

#NOTE: I tuned this number to get the calculation to match
mu_r_unknown = mu_0*700

# Measurements from my example solenoid
coreMu = mu_r_unknown
coreMass = 0.000785 #Kilograms
length = .035 #Meters
numLoops = 120 #Rough estimation
coreRadius = .0008 #Meters
crossSectionalArea = m.pi*(coreRadius**2)
resistance = 1.01
n = numLoops/length

# Since only a portion of the solenoid contains the core,
# we calculate a new mu = mu_0 * (% air) + mu_core * (% core)
def normalizedMu(lenCore):
    return (lenCore*coreMu + (length-lenCore)*coreMu)/length

# Find the inductance using the normalized mu ^^ given the
# length of the core that is enclosed in the solenoid
# http://hyperphysics.phy-astr.gsu.edu/hbase/magnetic/indcur.html#c1
def inductance(lenCore):
    return crossSectionalArea*(n**2)*normalizedMu(lenCore)

# Becasue of self inductance, the current in a solenoid
# takes time to reach the applied current. 
# http://hyperphysics.phy-astr.gsu.edu/hbase/electric/indtra.html#c1
def newCurrent(appliedCurrent, time, inductance):
    return appliedCurrent*(1-m.exp(-(resistance*time/inductance)))


# Using my calculations from previous logs
def force(appliedCurrent, time, lenCore, opposingForce):
    i = newCurrent(appliedCurrent, time, inductance(lenCore))
    f = (crossSectionalArea*(coreMu-mu_0)/2)*((n*i)**2)
    return f - opposingForce

# Solves for current needed to give a force equal to the opposing force
# based on the same force calculation above.
def getHoldingCurrent(opposingForce):
    return m.sqrt((2*opposingForce)/(crossSectionalArea * (coreMu - mu_0)))/n

# F = ma, so a = (calculated force)/mass
def acceleration(appliedCurrent, time, lenCore, opposingForce):
    return force(appliedCurrent, time, lenCore, opposingForce)/coreMass

# Velocity is just acceleration at a given time: v = a*t, but
# in our case, a changes because of the moving core
def velocity(appliedCurrent, time, lenCore, opposingForce):
    return acceleration(appliedCurrent, time, lenCore, opposingForce)*time

# Change in length = velocity * time, 
def newLenCore(lenCore, speed, time):
    return lenCore + (speed*time)

# Position is vertical, so opposing force is mass*gravity
opposingForce = coreMass * gravity
# Set the current to what is needed to oppose gravity
# should cause acceleration to 0 out
I = getHoldingCurrent(opposingForce)
# Other test current should cause continuous acceleration
I2 = I+.01
# Starting at time = 0
t = 0
# Time step
step = .00001
# Initial velocity
v = 0
# lists to store points for the plot
time = []
y = []

# We're going to go for (100,000 * .00001) seconds (1 second)
for i in range(100000):
    v = velocity(I, t, v, opposingForce)*t
    # Comment out above, and uncomment brlow to show graph 
    # of current slightly above holding force, allowing it 
    # to keep accelerating upwards.
    # v = velocity(I2, t, v, opposingForce)*t

    #this means only plot every thousanth step
    if i%1000 == 0:
        time.append(t)
        y.append(v)
    t += step
    a = np.array(time)
    b = np.array(y)
print("calculated holding current:",I)
# Plot
pp.plot(a, b)
pp.xlabel("time")
pp.ylabel("velocity")
pp.show()
