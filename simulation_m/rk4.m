function X = rk4(f,X,U,DT,M)
	dt = DT/M;
	for i=1:M
		k1 = f(X, U);
		k2 = f(X + dt/2 * k1, U);
		k3 = f(X + dt/2 * k2, U);
		k4 = f(X + dt * k3, U);
		X=X+dt/6*(k1 +2*k2 +2*k3 +k4);
	end
end