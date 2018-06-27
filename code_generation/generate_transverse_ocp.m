function generate_transverse_ocp(target_dir,dims,dyn,tnf,ocpopts)
import casadi.*

% parameters
P = SX.sym('P', 0);

% dims
Ntr = ocpopts.Ntr;
nx = dims.nx;
nu = dims.nu;
ngk = 0;
ngN = 3;
nbxk = 0;
nbxN = 0;
nbuk = nu;
nbuN = 0;

Xk = SX.sym('Xk',nx);
Uk = SX.sym('Uk',nu);
XN = Xk;
UN = SX.sym('UN',0);
[Xik, Etak] = tnf.Phi(Xk);
Vk = dyn.V(Xk,Uk);
HESSk = SX.sym('HESSk',nx+nu,nx+nu);
HESSN = SX.sym('HESSN',nx,nx);
LMUk = SX.sym('LMUk',ngk);
UMUk = SX.sym('UMUk',ngk);
LMUN = SX.sym('LMUk',ngN);
UMUN = SX.sym('UMUk',ngN);
PIk = SX.sym('PIk',nx);

% gk
gk = SX(ngk,1);
dgk = jacobian(gk,[Xk;Uk]);
d2gk = HESSk + jacobian(dgk.'*(UMUk-LMUk), [Xk; Uk]);
gkfun = Function('trans_gk', {Xk, Uk, LMUk, UMUk, HESSk, P}, {d2gk, dgk(:,1:nx), dgk(:,nx+(1:nu)), -gk, -gk});
gkfun.generate('trans_gk.c');
movefile('trans_gk.c',taget_dir);

% gN
gN = [Xik.'*ocpopts.P*Xik; Etak(2); 4*dyn.g*(dyn.r^2+Xik(1)).^(3/2).*sin(Etak(3)) + (Xik(2).^2  - 2*(-ocpopts.K*Xik)'.*(dyn.r^2+Xik(1))).*sin(Etak(1)-Etak(3))];
dgN = jacobian(gN,XN);
d2gN = HESSN + jacobian(dgN.'*(UMUN-LMUN), XN);
gNfun = Function('trans_gN', {XN, UN, LMUN, UMUN, HESSN, P}, {d2gN, dgN(:,1:nx), SX(0,1), -gN, -gN});
gNfun.generate('trans_gN.c');
movefile('trans_gN.c',target_dir);

% jk
jk = dyn.stage_cost(Xk,Uk,ocpopts.Q,ocpopts.R,ocpopts.dt);
[d2jk,djk] = hessian(jk,[Xk;Uk]);
d2jk = HESSk + d2jk;
jkfun = Function('trans_jk', {Xk, Uk, HESSk, P}, {d2jk, djk(1:nx), djk(nx+(1:nu)), jk});
jkfun.generate('trans_jk.c');
movefile('trans_jk.c',target_dir);

% jN
jN = Xik.'*ocpopts.P*Xik;
[d2jN,djN] = hessian(jN,XN);
d2jN = HESSN + d2jN;
jNfun = Function('trans_jN', {XN, UN, HESSN, P}, {d2jN, djN(1:nx), SX(0,1), jN});
jNfun.generate('trans_jN.c');
movefile('trans_jN.c',target_dir);

% fk
fk = dyn.xplus(Xk, Uk, ocpopts.dt);
dfk = jacobian(fk, [Xk;Uk]);
d2fk = HESSk + jacobian(dfk.'*PIk, [Xk; Uk]);
fkfun = Function('trans_fk', {Xk, Uk, PIk, HESSk, P}, {d2fk, dfk(:,1:nx), dfk(:,nx+(1:nu)), fk});
fkfun.generate('trans_fk.c');
movefile('trans_fk.c',target_dir);

% Collect bounds
nboffset = zeros(Ntr+2,1);
ngoffset = zeros(Ntr+2,1);
idxbdata = [];
ubdata = [];
lbdata = [];
ugdata = [];
lgdata = [];
for i=1:Ntr
    % Bounds on inputs
    idxbdata = [idxbdata;(0:nu-1)'];
    ubdata = [ubdata;dyn.ubu];
    lbdata = [lbdata;dyn.lbu];
    nboffset(i+1) = length(ubdata);
end
ugdata = [ugdata; ocpopts.e; 0.001; 0.001];
lgdata = [lgdata; -1; -0.001; -0.001];
numb = length(idxbdata);
numg = length(ugdata);

% Generate header file with bounds and dimensions
dims_str = '#ifndef TRANS_DIMS_BOUNDS_H_\n#define TRANS_DIMS_BOUNDS_H_\n\n';
%dims
dims_str = [dims_str, sprintf('static inline int getN() {return %d;}\n', Ntr)];
dims_str = [dims_str, sprintf('static inline int getNP() {return %d;}\n', size(P,1))];
dims_str = [dims_str, sprintf('static int nx [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nx)];
end
dims_str = [dims_str, sprintf('%d};\n', nx)];
dims_str = [dims_str, sprintf('static int qpnx [%d] = {', Ntr+1)];
dims_str = [dims_str, '0, '];
for i=2:Ntr
    dims_str = [dims_str, sprintf('%d, ', nx)];
end
dims_str = [dims_str, sprintf('%d};\n', nx)];
dims_str = [dims_str, sprintf('static int nu [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nu)];
end
dims_str = [dims_str, sprintf('%d};\n', 0)];
dims_str = [dims_str, sprintf('static int nbx [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nbxk)];
end
dims_str = [dims_str, sprintf('%d};\n', nbxN)];
dims_str = [dims_str, sprintf('static int nbu [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nbuk)];
end
dims_str = [dims_str, sprintf('%d};\n', nbuN)];
dims_str = [dims_str, sprintf('static int nb [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nbxk+nbuk)];
end
dims_str = [dims_str, sprintf('%d};\n', nbxN+nbuN)];
dims_str = [dims_str, sprintf('static int ng [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', ngk)];
end
dims_str = [dims_str, sprintf('%d};\n', ngN)];
dims_str = [dims_str, sprintf('static int ns [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', 0)];
end
dims_str = [dims_str, sprintf('%d};\n', 0)];
%bounds
dims_str = [dims_str, '\n\n'];
dims_str = [dims_str, sprintf('static int idxbdata [%d] = {',numb)];
for i=1:numb-1
    dims_str = [dims_str, sprintf('%d, ', idxbdata(i))];
end
dims_str = [dims_str, sprintf('%d};\n', idxbdata(numb))];
dims_str = [dims_str, sprintf('static double lbdata [%d] = {',numb)];
for i=1:numb-1
    dims_str = [dims_str, sprintf('%d, ', lbdata(i))];
end
dims_str = [dims_str, sprintf('%d};\n', lbdata(numb))];
dims_str = [dims_str, sprintf('static double ubdata [%d] = {',numb)];
for i=1:numb-1
    dims_str = [dims_str, sprintf('%d, ', ubdata(i))];
end
dims_str = [dims_str, sprintf('%d};\n', ubdata(numb))];
dims_str = [dims_str, sprintf('static double ugdata [%d] = {',numg)];
for i=1:numg-1
    dims_str = [dims_str, sprintf('%d, ', ugdata(i))];
end
dims_str = [dims_str, sprintf('%d};\n', ugdata(numg))];
dims_str = [dims_str, sprintf('static double lgdata [%d] = {',numg)];
for i=1:numg-1
    dims_str = [dims_str, sprintf('%d, ', lgdata(i))];
end
dims_str = [dims_str, sprintf('%d};\n', lgdata(numg))];
dims_str = [dims_str, sprintf('static int *idxb [%d] = {',Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('idxbdata+%d, ', nboffset(i))];
end
dims_str = [dims_str, sprintf('idxbdata+%d};\n', nboffset(Ntr+1))];
dims_str = [dims_str, sprintf('static double *ub [%d] = {',Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('ubdata+%d, ', nboffset(i))];
end
dims_str = [dims_str, sprintf('ubdata+%d};\n', nboffset(Ntr+1))];
dims_str = [dims_str, sprintf('static double *lb [%d] = {',Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('lbdata+%d, ', nboffset(i))];
end
dims_str = [dims_str, sprintf('lbdata+%d};\n', nboffset(Ntr+1))];
dims_str = [dims_str, sprintf('static double *ug [%d] = {',Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('ugdata+%d, ', ngoffset(i))];
end
dims_str = [dims_str, sprintf('ugdata+%d};\n', ngoffset(Ntr+1))];
dims_str = [dims_str, sprintf('static double *lg [%d] = {',Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('lgdata+%d, ', ngoffset(i))];
end
dims_str = [dims_str, sprintf('lgdata+%d};\n', ngoffset(Ntr+1))];
dims_str = [dims_str, '\n#endif // TRANS_DIMS_BOUNDS_H_\n'];
fid = fopen([target_dir,'/trans_dims_bounds.h'],'wt');
fprintf(fid, dims_str);
fclose(fid);

% Generate prototype
funtypes = {'trans_jk','trans_gk','trans_fk','trans_jN','trans_gN','trans_fN'};
proto_str = '#ifndef TRANS_FUNCTIONS_H_\n#define TRANS_FUNCTIONS_H_\n\n';
for ft=funtypes
    ft = ft{1};
    proto_str = [proto_str, sprintf('int %s(const double** arg, double** res, int* iw, double* w, void* mem);\n',ft)];
    proto_str = [proto_str, sprintf('const int* %s_sparsity_in(int i);\n',ft)];
    proto_str = [proto_str, sprintf('const int* %s_sparsity_out(int i);\n',ft)];
    proto_str = [proto_str, sprintf('int %s_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);\n',ft)];
    proto_str = [proto_str, '\n'];
end
proto_str = [proto_str, '\n#endif // TRANS_FUNCTIONS_H_\n'];
fid = fopen([target_dir,'/trans_functions.h'],'wt');
fprintf(fid, proto_str);
fclose(fid);



end

