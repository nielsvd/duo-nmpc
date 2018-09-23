import sqp_interface.*

% Configuration
nx = 5;
nu = 2;
dt = 0.02;
epsN = 0.1;
Nsim = 200;

% Parameters:
%  mostly regularization penalties in the economic objective
%  which tend to be necessary for working optimization code
ptr = {-1}; % dummy
pta = {5e-2, 5e-2, 5e-3}; % regxi, regeta, regu

% Initialize controller
%  Initial state
x0 = [1;0;0;0.1;0];
%  Initial guess of inputs
uig = [9.81; 0];
%  Construct controller object, args:
%    1. Initial state (for initializing SQP solution)
%    2. Initial guess for inputs (for initializing SQP solution)
%    3. Number of SQP iterations (same for both NLPs here)
%    4. Minimum eigenvalue clipped to in Hessian regularization scheme
duo_nmpc = DuoNmpc(num2cell(x0.'), num2cell(uig.'), 5, 1e-4);

% Pre-allocate the simulation, timing and kkt trajectories
X = zeros(nx,Nsim+1);
U = zeros(nu,Nsim);
CT = zeros(1,Nsim);
KKTtran = zeros(1,Nsim);
KKTtang = zeros(1,Nsim);
%  Assign initial state
X(:,1) = x0;

% Start simulation
for i=1:Nsim	
	% Execute sequential optimization scheme
	ct_start = tic;
	U(:,i) = cell2mat(duo_nmpc.Solve(num2cell(X(:,i).'),ptr,pta,epsN)).';
	CT(i) = toc(ct_start);

	% Simulate
	X(:,i+1) = quadrotor_xplus(X(:,i),U(:,i),dt);

	% Store optimality measures
	KKTtran(i) = duo_nmpc.GetTransverseKkt();
	KKTtang(i) = duo_nmpc.GetTangentialKkt();

	% As is common in RTI-type schemes, shift solution
	duo_nmpc.Shift();
end

% Sample epsilon bounds
%  Upper-left entry of inverse of P
Pi11 = 0.999530073198045;
%  Sample from -pi to pi
theta = linspace(-pi,pi,1000);
%  Compute epsilon corresponding to epsN
epsilon = sqrt(Pi11)*epsN;
%  Compute bounds
XepsLo = sqrt(1 - epsilon)*cos(theta);
XepsHi = sqrt(1 + epsilon)*cos(theta);
YepsLo = sqrt(1 - epsilon)*sin(theta);
YepsHi = sqrt(1 + epsilon)*sin(theta);

% Plot position trajectory
figure(10);
plot(XepsLo,YepsLo,'r--',XepsHi,YepsHi,'r--',X(1,:),X(3,:));
title('Position');
xlabel('$X$ [m]');
ylabel('$Z$ [m]');

% Plot input trajectory
figure(20);
subplot(2,1,1);
plot(U(1,:));
title('Thrust');
ylabel('F [N]');
subplot(2,1,2);
plot(U(2,:));
title('Rotation');
ylabel('w [rad/s]');

% Plot computation time
figure(30);
plot(CT*1000);
title('Computation time');
ylabel('[ms]');

% Plot KKT
figure(40);
subplot(2,1,1);
semilogy(KKTtran);
title('KKT value - transversal problem');
subplot(2,1,2);
semilogy(KKTtang);
title('KKT value - tangential problem');
