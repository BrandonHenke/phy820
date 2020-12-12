from copy import copy
from time import sleep

import numpy as np

from vpython import *


class AstroObject(sphere):
	"""
		AstroObject holds parameters (mass, G) and trajectory variables (position, 
		velocity, acceleration), routines for computing forces and energies, as well 
		as the visualization functionality inherited from its parent class, the 
		VPython 7 sphere class.

		parameters:
		-----------
		G: 	float
			Newton's constant in appropriate units
		
		mass: float
		
		pos, velocity, accel: VPython vector
			position, velocity and acceleration of the object 

		color: VPython symbol or vector
			pre-defined colors (color.{white, black, ...}) or vector(r,g,b) (normalized 
			to [0,1])

		radius: float
			controls the size of the rendered object - use appropriate units for visibility

		make_trail: Boolean
			show trajectory

		functions:
		----------
		gravity: vector
			computes the gravitational force between the current object and another instance
			of the class
		
		total_force: vector
			computes the total gravitational force exerted on the current object by a list
			of AstroObjects

		kinetic_energy, potential_energy: float
			as the name suggests

		update:
			updates the position, velocity and acceleration vectors
	"""

	def __init__(self, G=1, mass=1, pos=vector(0,0,0), velocity=vector(0,0,0),
		color=color.white, radius=50):
		print("Initiate 1")
		# This calls the constructor of the parent class - VPython sphere - with the shared
		# arguments.
		print(super())
		super().__init__(pos=pos, color=color, radius=radius, make_trail=True)
		print("Initiate 2")
		self.G        = G
		self.velocity = velocity
		self.mass     = mass
		self.T        = 0
		self.V        = 0
    
	
	def kinetic_energy(self):
		"""
			Computes the object's kinetic energy.
		"""
		T = 0.5 * self.mass * mag(self.velocity) ** 2
		return T
    
	def potential_energy(self, objects):
		"""
			Computes the object's potential energy ( = interaction energy).
		"""
		V = 0
		for other in objects:
			if other is not self:
				r_vec    = self.pos - other.pos
				distance = mag(r_vec)
				V -= self.G * self.mass * other.mass / distance

		return V

	def update(self, pos, velocity):
		"""
			Updates the object's position, velocity and/or acceleration vectors.
		"""
		self.pos = copy(pos)
		self.velocity = copy(velocity)
	


class Simulator():
	"""
		The Simulator class controls the simulation of an N-body system interacting
		via gravitational forces.

		parameters:
		-----------
		objects: 	AstroObjects[]
			The list of objects that make up the N-body system.
		
		dt: float
			Size of the time step for the evolution (in appropriate units).

		r_vec, v_vec: ndarrays[]
			Numpy arrays holding the objects' positions and vectors for interfacing 
			with the ODE solvers.

		functions:
		----------
		get_state, set_state: 
			copies the present state of the objects' variables from AstroObjects
			to the work space and back
		
		update_energies: 
			Updates the energies of all objects.

		update_euler, update_verlet, update_forrest_ruth: 
			Integration routines for the time evolution: Forward Euler, velocity
			Verlet, and fourth-order Forrest-Ruth algorithms.

		transform_to_com:
			Transforms positions and velocities into the center-of-mass frame
			of the system.

		total_mass: float
			Computes the total mass of the system.

		center_of_mass, total_momentum, total_angular_momentum: vector
			Compute the center of mass, total momentum, or total angular momentum
			of the system.
	"""

	def __init__(self, objects, G=1, dt=0.001):
		self.objects = objects

		self.G     = G
		self.dt    = dt 

		self.r_vec = np.zeros(3*len(self.objects))
		self.v_vec = np.zeros(3*len(self.objects))
		


	###########################################
	# utility functions
	###########################################

	def get_state(self):
		"""
			Copies the present state of the objects' variables from the list of
			AstroObjects to the work space.
		"""

		for i, obj in enumerate(self.objects):
		    self.r_vec[3*i:3*i+3] = np.array([obj.pos.x, obj.pos.y, obj.pos.z]).copy()
		    self.v_vec[3*i:3*i+3] = np.array([obj.velocity.x, obj.velocity.y, obj.velocity.z]).copy()

			# self.r_vec.append(copy(obj.pos))
			# self.v_vec.append(copy(obj.velocity))
			# self.a_vec.append(copy(obj.accel))

	def set_state(self):
		"""
			Copies the present state of the objects' variables from the work space
			to the AstroObjects after the time evolution step has been evaluated.
		"""
		for i, obj in enumerate(self.objects):
			obj.update(pos=vector(self.r_vec[3*i],self.r_vec[3*i+1],self.r_vec[3*i+2]), 
					   velocity=vector(self.v_vec[3*i],self.v_vec[3*i+1],self.v_vec[3*i+2]))

			# obj.update(pos=self.work_r_vec[i], velocity=self.work_v_vec[i], accel=self.work_a_vec[i])

	def update_energies(self):
		"""
			Causes AstroObjects managed by the Simulator to recompute their energies.
		"""
		for obj in self.objects:
			obj.T = obj.kinetic_energy()
			obj.V = obj.potential_energy(self.objects)


	###########################################
	# integrators
	###########################################


	def dvdt(self, r_vec):
		"""
			Computes the accelerations a_vec = dv_vec/dt of all objects as a 
			3*N dimensional NumPy array.
		"""

		dvdt_vec = np.zeros(3*len(self.objects))

		for i, obj in enumerate(self.objects):

			a_vec = np.zeros(3)
			
			for j, other in enumerate(self.objects):
				if obj is not other:
			
					dr_vec = r_vec[3*i:3*i+3] - r_vec[3*j:3*j+3]
					dr     = np.linalg.norm(dr_vec)
					a_vec  -= self.G * other.mass * dr_vec / (dr)**3
								
			dvdt_vec[3*i:3*i+3] = a_vec.copy()			

		return dvdt_vec


	# Euler / Runge-Kutta
	def update_euler(self):
		"""
			Evaluates a Forward-Euler time step with step size self.dt.
		"""

		self.get_state()

		# compute accelerations via the 3N-dimensional vectors
		# for i, obj in enumerate(self.objects):
		drdt_vec = self.v_vec
		dvdt_vec = self.dvdt(self.r_vec)

		# orward Euler step
		self.r_vec += drdt_vec * self.dt
		self.v_vec += dvdt_vec * self.dt
		
		# transfer new state to the individual objects, so the 3D scence can update
		self.set_state()
		self.update_energies()


	# symplectic methods
	def update_verlet(self):
		"""
			Evaluates a velocity Verlet time step with step size self.dt.
		"""
		self.get_state()

		# save for later
		a_vec_old = self.dvdt(self.r_vec)

		self.r_vec += self.v_vec * self.dt + 0.5 * a_vec_old * self.dt **2
		self.v_vec += 0.5 * (a_vec_old + self.dvdt(self.r_vec)) * self.dt

		# transfer new state to the individual objects
		self.set_state()
		self.update_energies()


	def update_forest_ruth(self):
		"""
			Evaluates a Forest-Ruth time step with step size self.dt.
		"""
		K = 1./(2-np.cbrt(2))

		self.get_state()

		self.r_vec += self.v_vec * 0.5 * K * self.dt

		self.v_vec += self.dvdt(self.r_vec) * K * self.dt
		self.r_vec += self.v_vec * 0.5 * (1. - K) * self.dt

		self.v_vec += self.dvdt(self.r_vec) * (1.-2*K) * self.dt
		self.r_vec += self.v_vec * 0.5 * (1.-K) * self.dt

		self.v_vec += self.dvdt(self.r_vec) * K * self.dt
		self.r_vec += self.v_vec * 0.5 * K * self.dt

		# transfer new state to the individual objects
		self.set_state()
		self.update_energies()


	###########################################
	# Physical Quantities
	###########################################

	def total_mass(self):
		"""
			Computes the total mass of the system.
		"""
		M = 0
		for obj in self.objects:
			M += obj.mass

		return M

	def center_of_mass(self):
		"""
			Computes the system's center of mass.
		"""
		R = vector(0,0,0)
		M = self.total_mass()

		for obj in self.objects:
			R += obj.mass/M * obj.pos

		return R

	def total_momentum(self):
		"""
			Computes the total momentum of the system.
		"""
		P = vector(0,0,0)

		for obj in self.objects:
			P += obj.mass * obj.velocity

		return P

	def total_angular_momentum(self):
		"""
			Computes the system's total angular momentum
		"""
		L = vector(0,0,0)

		for obj in self.objects:
			L += obj.mass * cross(obj.pos, obj.velocity)

		return L

	###########################################
	# Coordinate System
	###########################################
	def transform_to_com(self):
		"""
			Transform positions and velocities of the simulation objects to
			the center-of-mass frame.
		"""
		M = self.total_mass()
		R = self.center_of_mass()
		P = self.total_momentum()

		for obj in self.objects:
			obj.pos      -= R
			obj.velocity -= P/M
