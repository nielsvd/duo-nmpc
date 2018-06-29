function [dims,dyn,tnf] = quadrotor()
import casadi.*


% Time step
dt = SX.sym('dt');


% Dimensions
l = 1;
nu = 2;
nx = 5;
nxi = 2;
neta = 3;
dims = struct('l',l,'nu',nu,'nx',nx,'nxi',nxi,'neta',neta);


% Constants
m = 1;
g = 9.81;
r = 1;
dyn.m = m;
dyn.g = g;
dyn.r = r;

% States
x = SX.sym('x',nx);
u = SX.sym('u',nu);
xi = SX.sym('xi',nxi);
vtr = SX.sym('vtr',l);
eta = SX.sym('eta',neta);
vta = SX.sym('vta',nu-l);


% State-space system
f = [x(2);0;x(4);-g;0];
G = [0,0;-sin(x(5))/m,0;0,0;cos(x(5))/m,0;0,1];
xdot = f +G*u;
dyn.xdot = Function('xdot',{x,u},{xdot});
xplus =  rk4x(dyn.xdot,x,u,dt,1); % Integrate using RK4
dyn.xplus = Function('xplus',{x,u,dt},{xplus});


% Transverse system dynamics
dyn.Atr = [0 1; 0 0];
dyn.Btr = [0; 1];
dyn.Atr_d = [1 dt; 0 1];
dyn.Btr_d = [0.5*dt^2; dt];
dyn.xidot = Function('xidot',{xi,vtr},{dyn.Atr*xi + dyn.Btr*vtr});
dyn.xiplus = Function('xiplus',{xi,vtr,dt},{dyn.Atr_d*xi + dyn.Btr_d*vtr});


% Tangential system dynamics
f0 = [eta(2);
      g*sin(eta(3))/(sqrt(r^2+xi(1))*cos(eta(1)-eta(3))) + (-4*eta(2)*(r^2+xi(1))*xi(2) + (4*eta(2)^2*(r^2+xi(1))^2 + xi(2)^2)*tan(eta(1)-eta(3)))/(4*(r^2+xi(1))^2);
      0];
Gtr = [0;
       -tan(eta(1)-eta(3))/(2*(r^2+xi(1)));
       0];
Gta = [0;
       0;
       1];
dyn.f0 = Function('f0',{xi,eta},{f0});
dyn.Gtr = Function('Gtr',{xi,eta},{Gtr});
dyn.Gta = Function('Gta',{xi,eta},{Gta});
dyn.etadot = Function('etadot',{xi,eta,vtr,vta},{f0 + Gtr*vtr + Gta*vta});
etaplus =  rk4eta(dyn.etadot,xi,eta,vtr,vta,dt,1); % Integrate using RK4
dyn.etaplus = Function('etaplus',{xi,eta,vtr,vta,dt},{etaplus});


% Actuator limits
dyn.lbu = [  0; -pi];
dyn.ubu = [ 15;  pi];


% TNF
tnf.h = @(x,z) x^2 + z^2 - r^2;
Phi_xi = [x(1)^2 + x(3)^2 - r^2; 
          2*x(1)*x(2) + 2*x(3)*x(4)];
Phi_eta = [atan2(-x(1),x(3));
           (-x(3)*x(2)+x(1)*x(4))/(x(1)^2+x(3)^2);
           x(5)];
Phi_inv = [-sin(eta(1))*sqrt(xi(1)+r^2);
           -(2*eta(2)*cos(eta(1))*(xi(1)+r^2) + xi(2)*sin(eta(1)))/(2*sqrt(xi(1)+r^2));
           cos(eta(1))*sqrt(xi(1)+r^2);
           (xi(2)*cos(eta(1)) - 2*eta(2)*sin(eta(1))*(xi(1)+r^2))/(2*sqrt(xi(1)+r^2));
           eta(3)];
tnf.Phi = Function('Phi',{x},{Phi_xi,Phi_eta});
tnf.Phi_inv = Function('Phi_inv',{xi,eta},{Phi_inv});


% Transversely linearizing feedback (Nielsen et al., 2010)
beta_x = [m/(2*x(3)*cos(x(5)) - 2*x(1)*sin(x(5))), 0; 0, 1];
alpha_x = [-(m*(x(2)^2-g*x(3)+x(4)^2))/(x(3)*cos(x(5))-x(1)*sin(x(5))); 0];
tnf.alpha_x = Function('alpha_x',{x},{alpha_x});
tnf.beta_x = Function('beta_x',{x},{beta_x});
alpha = [-(m*(4*(r^2+xi(1))*(-g*cos(eta(1))*sqrt(r^2+xi(1)) + eta(2)^2*(r^2+xi(1))) + xi(2)^2))/(4*(r^2+xi(1))^(3/2)*cos(eta(1)-eta(3))); 0];
beta = [m/(2*sqrt(r^2+xi(1))*cos(eta(1)-eta(3))), 0; 0, 1];
tnf.alpha = Function('alpha',{xi,eta},{alpha});
tnf.beta = Function('beta',{xi,eta},{beta});
tnf.kappa = Function('kappa',{xi,eta,vtr,vta},{alpha + beta*[vtr;vta]});


% Add LQ integrator to dynamics struct
Q = SX.sym('Q',nxi,nxi);
R = SX.sym('R',l,l);
alphadd = 2*x(2)^2 + 2*x(4)^2 + 2*x(1)*xdot(2) + 2*x(3)*xdot(4);
lagcost = Phi_xi.'*Q*Phi_xi + alphadd.'*R*alphadd;
lagcostfun = Function('lagcost',{x,u},{xdot, lagcost});
dyn.lag_cost = Function('lagcost',{x,u},{lagcost});
stage_cost = rk4q(lagcostfun,x,u,dt,1);
dyn.stage_cost = Function('stage_cost',{x,u,Q,R,dt},{stage_cost}); 

end

