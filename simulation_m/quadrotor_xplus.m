function xplus = quadrotor_xplus(x,u,dt)
	f = @(x) [x(2);0;x(4);-9.81;0];
	G = @(u) [0, 0; -sin(x(5)), 0; 0, 0; cos(x(5)), 0; 0, 1];
	xdot = @(x,u) f(x) + G(x)*u;

	xplus =  rk4(xdot,x,u,dt,4);
end