function ETA = rk4eta(f, XI, ETA0, VTR, VTA,  DT, M)

% Fixed step Runge-Kutta 4 integrator
ETA = ETA0;
for j=1:M
    k1 = f(XI, ETA, VTR, VTA);
    k2 = f(XI, ETA + DT/2 * k1, VTR, VTA);
    k3 = f(XI, ETA + DT/2 * k2, VTR, VTA);
    k4 = f(XI, ETA + DT * k3, VTR, VTA);
    ETA=ETA+DT/6*(k1 +2*k2 +2*k3 +k4);
end


end

