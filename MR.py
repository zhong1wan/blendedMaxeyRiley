import numpy as np

class MaxeyRiley(object):
	"""docstring for prtcl"""
	def __init__(self, ep, R, flow):
		self.ep = ep
		self.mu = 1/self.ep
		self.R = R
		self.flow = flow
		
	def dynamics(self, x, v, t=0):
		
		u = self.flow.vel(x, t)
		DuDt = self.flow.DuDt(x, t=t)

		a = 3/2*self.R*DuDt - self.mu*(v - u)

		return a

	def order1_manifold_vel(self, x, t=0):

		u = self.flow.vel(x, t)
		DuDt = self.flow.DuDt(x, t=t)

		return u + self.ep*(3/2*self.R - 1)*DuDt

	def RK4step(self, x, v, dt, t=0):

		a1 = self.dynamics(x, v, t)
    
		x2 = x + .5*dt*v
		v2 = v + .5*dt*a1
		a2 = self.dynamics(x2, v2, t+.5*dt)
    
		x3 = x + .5*dt*v2
		v3 = v + .5*dt*a2
		a3 = self.dynamics(x3, v3, t+.5*dt)
    
		x4 = x + dt*v3
		v4 = v + dt*a3
		a4 = self.dynamics(x4, v4, t + dt)
    
		x_new = x + dt*(v + 2*v2 + 2*v3 + v4)/6
		v_new = v + dt*(a1 + 2*a2 + 2*a3 + a4)/6

		return x_new, v_new

	def euler1step(self, x, v, dt, t=0):

		x_new = x + dt*v
		a = self.dynamics(x, v, t)
		v_new = v + dt*a

		return x_new, v_new

