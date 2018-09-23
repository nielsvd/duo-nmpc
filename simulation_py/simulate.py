from sqp_interface import *
from numpy import *
from quadrotor import quadrotor_xplus
import matplotlib.pyplot as plt
from timeit import default_timer as timer

# Configuration
nx = 5
nu = 2
dt = 0.02
epsN = 0.1
Nsim = 200

# Parameters:
##  mostly regularization penalties in the economic objective
##  which tend to be necessary for working optimization code
ptr = []
pta = [5e-2, 5e-2, 5e-3] # regxi, regeta, regu

# Initialize controller
## Initial state
x0 = [1,0,0,0.1,0]
## Initial guess of inputs
uig = [9.81, 0]
## Construct controller object, args:
### 1. Initial state (for initializing SQP solution)
### 2. Initial guess for inputs (for initializing SQP solution)
### 3. Number of SQP iterations (same for both NLPs here)
### 4. Minimum eigenvalue clipped to in Hessian regularization scheme
duo_nmpc = DuoNmpc(x0, uig, 3, 1e-4)

# Pre-allocate the simulation, timing and kkt trajectories
X = zeros((nx,Nsim+1))
U = zeros((nu,Nsim))
CT = zeros(Nsim)
KKTtran = zeros(Nsim)
KKTtang = zeros(Nsim)
## Assign initial state
X[:,0] = x0

print()

# Start simulation
for i in range(0,Nsim):	
	## Execute sequential optimization scheme
	ct_start = timer()
	U[:,i] = duo_nmpc.Solve(X[:,i].tolist(),ptr,pta,epsN)
	CT[i] = timer() - ct_start

	## Simulate
	X[:,i+1] = quadrotor_xplus(X[:,i],U[:,i],dt)

	## Print optimality measures
	KKTtran[i] = duo_nmpc.GetTransverseKkt()
	KKTtang[i] = duo_nmpc.GetTangentialKkt()

	## As is common in RTI-type schemes, shift solution
	duo_nmpc.Shift()

# Sample epsilon bounds
## Upper-left entry of inverse of P
Pi11 = 0.999530073198045
## Sample from -pi to pi
theta = linspace(-pi,pi,1000)
## Compute epsilon corresponding to epsN
epsilon = sqrt(Pi11)*epsN
## Compute bounds
XepsLo = sqrt(1 - epsilon)*cos(theta)
XepsHi = sqrt(1 + epsilon)*cos(theta)
YepsLo = sqrt(1 - epsilon)*sin(theta)
YepsHi = sqrt(1 + epsilon)*sin(theta)

# Plot position trajectory
plt.figure()
plt.plot(XepsLo,YepsLo,'r--',XepsHi,YepsHi,'r--',X[0,:],X[2,:])
plt.title('Position')
plt.xlabel('X [m]')
plt.ylabel('Z [m]')

# Plot input trajectory
plt.figure()
plt.subplot(2,1,1)
plt.plot(U[0,:])
plt.title('Thrust')
plt.ylabel('F [N]')
plt.subplot(2,1,2)
plt.plot(U[1,:])
plt.title('Rotation')
plt.ylabel('w [rad/s]')

# Plot computation time
plt.figure()
plt.plot(CT*1000)
plt.title('Computation time')
plt.ylabel('[ms]')

# Plot KKT
plt.figure()
plt.subplot(2,1,1)
plt.semilogy(KKTtran)
plt.title('KKT value - transversal problem')
plt.subplot(2,1,2)
plt.semilogy(KKTtang)
plt.title('KKT value - tangential problem')

plt.show()
