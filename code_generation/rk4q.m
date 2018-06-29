function Q = rk4q(f, X0, U, DT, M)

dt = DT/M;

% Fixed step Runge-Kutta 4 integrator
X = X0;
Q = 0;
for j=1:M
    [k1, k1_q] = f(X, U);
    [k2, k2_q] = f(X + dt/2 * k1, U);
    [k3, k3_q] = f(X + dt/2 * k2, U);
    [k4, k4_q] = f(X + dt * k3, U);
    X=X+dt/6*(k1 +2*k2 +2*k3 +k4);
    Q=Q+dt/6*(k1_q +2*k2_q +2*k3_q +k4_q);
end


end

