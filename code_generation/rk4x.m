function X = rk4x(f, X0, U, DT, M)

% Fixed step Runge-Kutta 4 integrator
X = X0;
for j=1:M
    k1 = f(X, U);
    k2 = f(X + DT/2 * k1, U);
    k3 = f(X + DT/2 * k2, U);
    k4 = f(X + DT * k3, U);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
end


end

