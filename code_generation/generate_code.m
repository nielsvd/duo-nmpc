clear;clc;
import casadi.*

% Obtain model, mappings, and the transversely linearizing feedback
[dims,dyn,tnf] = quadrotor();

% Generate code for NLPs
ocpopts.Ntr = 50;
ocpopts.Nta = 40;
ocpopts.dt = 0.02;
ocpopts.Q = diag([1 1]);
ocpopts.R = 0.001;
ocpopts.P = care(dyn.Atr,dyn.Btr,ocpopts.Q,ocpopts.R);
ocpopts.K = lqr(dyn.Atr,dyn.Btr,ocpopts.Q,ocpopts.R);
ocpopts.e = 0.0033;
ocpopts.M = 5;
% Target directory
target_dir = '../controller_library/codegen'
% Empty funs-folder
delete([target_dir,'/*.o');
delete([target_dir,'/*.c');
delete([target_dir,'/*.h');
% Set-up transverse OCP
generate_transverse_ocp(target_dir,dims,dyn,tnf,ocpopts);
% Set-up tangential OCP
generate_tangential_ocp(target_dir,dims,dyn,tnf,ocpopts);