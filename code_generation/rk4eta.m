%
%    This file is part of duo-nmpc.
%
%    Duo-nmpc
%    Copyright (C) 2018 Niels van Duijkeren, KU Leuven.
%
%    Duo-nmpc is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    Duo-nmpc is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with duo-nmpc; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%

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

