#
#    This file is part of duo-nmpc.
#
#    Duo-nmpc
#    Copyright (C) 2018 Niels van Duijkeren, KU Leuven.
#
#    Duo-nmpc is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 3 of the License, or (at your option) any later version.
#
#    Duo-nmpc is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with duo-nmpc; if not, write to the Free Software Foundation,
#    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#

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