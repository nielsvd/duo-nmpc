from numpy import *

def rk4(f,X,U,DT,M):
	dt = DT/M
	for i in range(0,M):
		k1 = f(X, U)
		k2 = f(X + dt/2 * k1, U)
		k3 = f(X + dt/2 * k2, U)
		k4 = f(X + dt * k3, U)
		X=X+dt/6*(k1 +2*k2 +2*k3 +k4)

	return X

def quadrotor_xplus(x,u,dt):
	f = lambda x: array([x[1],0.,x[3],-9.81,0.])
	G = lambda u: array([[0.,0.],[-sin(x[4]),0.],[0.,0.],[cos(x[4]),0.],[0.,1.]])
	xdot = lambda x, u: f(x) + G(x).dot(u)

	xplus =  rk4(xdot,x,u,dt,4);

	return xplus