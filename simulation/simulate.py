from sqp_interface import *
from math import sqrt

# Configuration
dt = 0.02
epsN = 0.1
Nsim = 500

# Parameters
ptr = []
pta = [5e-2, 5e-2, 5e-3] # regxi, regeta, regu

# Initialize controller
x0 = [1,0,0,0.1,0]
uig = [9.81, 0]
duo_nmpc = DuoNmpc(x0, uig, 3, 1e-4)

usim = duo_nmpc.Solve(x0,ptr,pta,epsN)

print(usim)