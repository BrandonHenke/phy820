from copy import copy
from time import sleep

import numpy as np
from p02_simulator import AstroObject, Simulator


from vpython import *


#################################################################################
# Main program
#################################################################################

# set up parameters
G    = 6.67*10**(-11)   # Newton's gravitational constant in m**3 kg**(-1) s**(-2)
M    = 1.99*10**30      # mass of the Sun in kg
m    = 5.97*10**24      # mass of the Earth in kg
Rmin = 147.1*10**9      # perihelion distance (initial point) in m

# express everything in natural units - use years for time
m0=5.97*10**24     # express all masses in terms of Earth's mass
R0=149.6*10**8     # 1/10 AU (experiment with this)
t0=24*3600*365.24

G=G/(R0**3) *m0 * t0**2  # G in R0^3 m_E**(-1) years**(-2) 
M    = M/m0
m    = m/m0
Rmin = Rmin/R0


# set up the VPython scene
scene = canvas(title='Solar System',
            width=600, height=400,
            center=vector(0,0,0), background=color.black)

# For some reason, the creators thought it would be a good idea to have y be the
# upward direction. We'll change that to the z direction.
scene.forward = vector(1,0,0)
scene.up = vector(0,0,1)

# Define and initiate the simulated objects - remember to translate to natural units!
sun     = AstroObject(G, 
                      mass = 1.99*10**30/m0, 
                      pos=vector(0,0,0),
                      velocity=vector(0,0,0), 
                      color=color.orange, radius=1)

earth   = AstroObject(G, 
                      mass = 5.97*10**24/m0, 
                      pos=vector(147.1*10**9/R0,0,0), 
                      velocity=vector(0,29800*t0/R0,0), 
                      color=color.blue, radius=0.2)

# Create the list of objects and initiate the simulator.
objects=[sun, earth]
sim = Simulator(objects, G, 0.001)


# Choose the time span of the simulation (in years).
tmax = 2


# Create a VPython graph object for the potential energy.
vgraph = graph(x=800, y=0,width=600,height=600,\
              title = 'Potential Energy', \
              xtitle = 't [yr]', ytitle = 'V [m_E  R0^2 yr^-2]', \
              foreground = color.black, background =color.white, \
              xmax = tmax, xmin = 0)

# All subsequently defined VPython curve objects are children of the
# same graph until the next graph object is created.
vcurves=[ ]
for obj in objects:
	vcurves.append(gcurve(color=obj.color))

# Same graph for the kinetic energy...
tgraph = graph(x=800, y=0,width=600,height=600,\
                  title = 'Kinetic Energy (radial + angular)', \
                  xtitle = 't [yr]', ytitle = 'T [m_E  R0^2 yr^-2]', \
              foreground = color.black, background =color.white, \
              xmax = tmax, xmin = 0)

tcurves=[ ]
for obj in objects:
	tcurves.append(gcurve(color=obj.color))


# ... and the total energy.
egraph = graph(x=800, y=0,width=600,height=600,\
                  title = 'Total Energy', \
                  xtitle = 't [yr]', ytitle = 'E [m_E  R0^2 yr^-2]', \
              foreground = color.black, background =color.white, \
              xmax = tmax, xmin = 0)

ecurves=[ ]
for obj in objects:
	ecurves.append(gcurve(color=obj.color))

# We add one curve that contains the total energy of the entire system.
ecurves.append(gcurve(color=vector(31,158,137)/255.))


# Initialize step counter...
steps = 0

# ... and start the simulation.
while steps * sim.dt < tmax:

	# VPython animation rate.
	rate(100)
	
	# Take a time step.
	sim.update_euler()

	# Update energy graphs.
	totE = 0
	for i, obj in enumerate(objects):
		# if obj == earth:
		tcurves[i].plot(steps*sim.dt, obj.T)
		vcurves[i].plot(steps*sim.dt, obj.V)
		ecurves[i].plot(steps*sim.dt, obj.V + obj.T)
		totE += obj.T + 0.5*obj.V

	ecurves[-1].plot(steps*sim.dt, totE)
        
	steps+=1


