#####################################################
# Author: Joyjit Kundu	                            #
# System: Ising Model in two dimension              #
# Version: Python v.3.5                             #
# Programming Approach: Object Oriented Programming #
#####################################################
#!/usr/bin/python3

import numpy as np
from numpy.random import *

seed()

class ising():
	mag = 0.0
	tot_E = 0.0
	def __init__(self, L, T):
		self.L = L
		self.T = T

	# L denotes the linear size of the system and T is the temperature	

	# Initializing the spin configuration

	def initialize(self):
		self.config = np.zeros((self.L, self.L), dtype = int)
		self.lat = np.zeros(self.L*self.L, dtype = int)
		for i in range(self.L):
			for j in range(self.L):
				k = j+i*self.L
				if np.random.uniform(0,1) > 2:
					self.config[i][j] = -1
					self.lat[k] = -1
				else:
					self.config[i][j] = 1
					self.lat[k] = 1
		ising.mag = 0.0
		ising.mag = np.sum(self.lat)

	# Converting the 2d array into an 1d array

	def find_site(self,x,y):
		self.x = (x+self.L) % self.L
		self.y = (y+self.L) % self.L
		return self.x + self.y * self.L 

	# Creating neighbour list for each site

	def cal_neighbour(self):
		self.N = self.L*self.L
		self.ln = np.zeros(self.N, dtype = int)
		self.rn = np.zeros(self.N, dtype = int)
		self.bn = np.zeros(self.N, dtype = int)
		self.tn = np.zeros(self.N, dtype = int)
		for i in range(self.L):
			for j in range(self.L):
				self.site = self.find_site(i,j)
				self.ln[self.site] = self.find_site(i-1,j)
				self.rn[self.site] = self.find_site(i+1,j)
				self.tn[self.site] = self.find_site(i,j-1)
				self.bn[self.site] = self.find_site(i,j+1)

	# Calculates the total energy

	def cal_energy(self):
		for i in range(self.L):
			for j in range(self.L):
				self.site = self.find_site(i,j)
				ising.tot_E = ising.tot_E + self.config[i][j]*(self.lat[self.ln[self.site]]+self.lat[self.rn[self.site]]+self.lat[self.tn[self.site]]+self.lat[self.bn[self.site]])	
		ising.tot_E = ising.tot_E/2.0
		return ising.tot_E

	# Whether to flip a spin or not based on Metropolis criterion

	def spin_flip(self):
		self.x = randint(0, self.L)
		self.y = randint(0, self.L)
		self.s = self.config[self.x, self.y]
		self.site = self.find_site(self.x,self.y)
		self.lener = self.s * (self.lat[self.ln[self.site]]+self.lat[self.rn[self.site]]+self.lat[self.tn[self.site]]+self.lat[self.bn[self.site]])
		self.dE = 2.0*self.lener
		if(self.dE < 0 or np.random.uniform(0,1)< np.exp(-self.dE/self.T)):
			self.dM = -2.0*self.config[self.x, self.y]
			self.config[self.x, self.y] = -self.config[self.x, self.y] 
			self.lat[self.site]=-self.lat[self.site]
			ising.mag += self.dM			
			ising.tot_E = ising.tot_E + self.dE

	# Monte Carlo seep

	def MC_steps(self):
		step = self.L*self.L
		for i in range(step):
			self.spin_flip()

