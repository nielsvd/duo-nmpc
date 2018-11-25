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

clear;clc;
import casadi.*

% Obtain model, mappings, and the transversely linearizing feedback
[dims,dyn,tnf] = quadrotor();

% Generate code for NLPs
ocpopts.Ntr = 50;
ocpopts.Nta = 40; % Nta > Ntr is currently not supported!
ocpopts.dt = 0.02;
ocpopts.Q = diag([1 1]);
ocpopts.R = 0.001;
ocpopts.P = care(dyn.Atr,dyn.Btr,ocpopts.Q,ocpopts.R);
ocpopts.K = lqr(dyn.Atr,dyn.Btr,ocpopts.Q,ocpopts.R);
ocpopts.e = 0.0033;
ocpopts.M = 5;
% Target directory
target_dir = '../controller_library/codegen';
% Empty funs-folder
delete([target_dir,'/*.o']);
delete([target_dir,'/*.c']);
delete([target_dir,'/*.h']);
% Set-up transverse OCP
generate_transverse_ocp(target_dir,dims,dyn,tnf,ocpopts);
% Set-up tangential OCP
generate_tangential_ocp(target_dir,dims,dyn,tnf,ocpopts);