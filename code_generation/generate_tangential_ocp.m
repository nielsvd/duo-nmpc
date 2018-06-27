function generate_tangential_ocp(target_dir,dims,dyn,tnf,ocpopts)
import casadi.*

% stage cost
uref = [9.81; 0];
regxi = SX.sym('regxi',1);
regeta = SX.sym('regeta',1);
regu = SX.sym('regu');
P = [regxi;regeta;regu];
stage_cost = @(x) (-x(3)*x(2) + x(1)*x(4))/(x(1)^2 + x(3)^2);

% dims
Ntr = ocpopts.Ntr;
Nta = ocpopts.Nta;
N = max([Nta,Ntr]);
nx = zeros(N+1,1);
nu = zeros(N+1,1);
ng = zeros(N+1,1);
nx(1) = dims.nx;
nu(1) = dims.nu;
ng(1) = 0;
for i = 2:N+1
    % Primary dynamics
    if i <= Nta
        nx(i) = nx(i) + dims.nx;
        nu(i) = nu(i) + dims.nu;
        % Tunnel constraint
        ng(i) = ng(i) + 1;
    elseif i == Nta+1
        nx(i) = nx(i) + dims.nx;
        % Tunnel constraint
        ng(i) = ng(i) + 1;
    end
    % Secodary dynamics
    if i <= Ntr
        nx(i) = nx(i) + dims.nx + 1;
        nu(i) = nu(i) + dims.nu;
    elseif i == Ntr+1
        nx(i) = nx(i) + dims.nx + 1;
        % Lyapunov constraint (1) + terminal xi constraint (1) + terminal eta constraint (2)
        ng(i) = ng(i) + 4;
    end
end

% inputs
X = cell(N+1,1);U = cell(N+1,1);LMU = cell(N+1,1);UMU = cell(N+1,1);PI = cell(N,1);HESS = cell(N+1,1);
for i = 1:N+1
    X{i} = SX.sym('Xk', nx(i));
    U{i} = SX.sym('Uk', nu(i));
    LMU{i} = SX.sym('LMUk', ng(i));
    UMU{i} = SX.sym('UMUk', ng(i));
    HESS{i} = SX.sym('HESSk', nx(i)+nu(i), nx(i)+nu(i));
end
for i = 1:N
    PI{i} = SX.sym('PIk', nx(i+1));
end

% gk & fk
Xk = X{1};
Uk = U{1};
LMUk = LMU{1};
UMUk = UMU{1};
PIk = PI{1};
HESSk = HESS{1};
gk = [];
f0 = dyn.xplus(Xk(1:dims.nx), Uk(1:dims.nu), ocpopts.dt);
fk = [f0;f0;0];
% Generate constraint function
dgk = jacobian(gk, [Xk; Uk]);
d2gk = HESSk + jacobian(dgk.'*(UMUk-LMUk), [Xk; Uk]);
gkfun = Function('tang_g0', {Xk, Uk, LMUk, UMUk, HESSk, P}, {d2gk, dgk(:,1:nx(1)), dgk(:,nx(1)+(1:nu(1))), -gk, -gk});
gkfun.generate('tang_g0.c');
movefile('tang_g0.c',target_dir);
% Generate rhs of dynamics
dfk = jacobian(fk, [Xk; Uk]);
d2fk = HESSk + jacobian(dfk.'*PIk, [Xk; Uk]);
fkfun = Function('tang_f0', {Xk, Uk, PIk, HESSk, P}, {d2fk, dfk(:,1:nx(1)), dfk(:,nx(1)+(1:nu(1))), fk});
fkfun.generate('tang_f0.c');
movefile('tang_f0.c',target_dir);
% Generate cost
[Xik,Etak] = tnf.Phi(Xk(1:dims.nx));
jk = stage_cost(Xk(1:dims.nx));
jk = jk + regxi*(Xik.')*Xik + regeta*(Etak(2:end).')*Etak(2:end) + regu*(Uk(1:dims.nu)-uref).'*(Uk(1:dims.nu)-uref);
[d2jk, djk] = hessian(jk, [Xk; Uk]);
d2jk = HESSk + d2jk;
jkfun = Function('tang_j0', {Xk, Uk, HESSk, P}, {d2jk, djk(1:nx(1)), djk(nx(1)+(1:nu(1))), jk});
jkfun.generate('tang_j0.c');
movefile('tang_j0.c',target_dir);
for i = 2:N+1
    jk = 0;
    gk = [];
    fk = [];
    Xk = X{i};
    Uk = U{i};
    LMUk = LMU{i};
    UMUk = UMU{i};
    HESSk = HESS{i};
    if i < N+1
        PIk = PI{i};
    end

    % Primary dynamics
    offsetx = 0;
    offsetu = 0;
    if i <= Nta
        [Xik,Etak] = tnf.Phi(Xk(1:dims.nx));
        gk = [gk; Xik.'*ocpopts.P*Xik];
        fk = [fk;dyn.xplus(Xk(1:dims.nx), Uk(1:dims.nu), ocpopts.dt)];
        offsetx = dims.nx;
        offsetu = dims.nu;
        % objective
        jk = jk + stage_cost(Xk(1:dims.nx));
        % Regularize a little
        jk = jk + regxi*(Xik.')*Xik + regeta*(Etak(2:end).')*Etak(2:end) + regu*(Uk(1:dims.nu)-uref).'*(Uk(1:dims.nu)-uref);
    elseif i == Nta+1
        [Xik,Etak] = tnf.Phi(Xk(1:dims.nx));
        gk = [gk; Xik.'*ocpopts.P*Xik];
        offsetx = dims.nx;
        offsetu = 0;
        % objective
        jk = jk + stage_cost(Xk(1:dims.nx));
        % Regularize a little
        jk = jk + regxi*(Xik.')*Xik + regeta*(Etak(2:end).')*Etak(2:end);
    end
    % Secondary dynamics
    idx_gklyap = -1;
    idx_fklyap = -1;
    if i <= Ntr
        fklyap = Xk(offsetx+dims.nx+1) + dyn.stage_cost(Xk(offsetx+(1:dims.nx)), Uk(offsetu+(1:dims.nu)), ocpopts.Q, ocpopts.R, ocpopts.dt);
        % System dynamics
        fk = [fk; dyn.xplus(Xk(offsetx+(1:dims.nx)), Uk(offsetu+(1:dims.nu)), ocpopts.dt)];
        fk = [fk; fklyap];
    elseif i == Ntr+1
        [Xihatk,Etahatk] = tnf.Phi(Xk(offsetx+(1:dims.nx)));
        % Lyapunov constraint
        gk = [gk; Xk(offsetx+dims.nx+1) + Xihatk.'*ocpopts.P*Xihatk];
        % Terminal xi-constraint
        gk = [gk; Xihatk.'*ocpopts.P*Xihatk];
        % Terminal eta-constraint
        gk = [gk; Etahatk(2)];
        gk = [gk; 4*dyn.g*(dyn.r^2+Xihatk(1)).^(3/2).*sin(Etahatk(3)) + (Xihatk(2).^2  - 2*(-ocpopts.K*Xihatk)'.*(dyn.r^2+Xihatk(1))).*sin(Etahatk(1)-Etahatk(3))];
    end
    % Generate constraint function
    dgk = jacobian(gk, [Xk; Uk]);
    d2gk = HESSk + jacobian(dgk.'*(UMUk-LMUk), [Xk; Uk]);       
    gkfun = Function(sprintf('tang_g%d',i-1), {Xk, Uk, LMUk, UMUk, HESSk, P}, {d2gk, dgk(:,1:nx(i)), dgk(:,nx(i)+(1:nu(i))), -gk, -gk});
    gkfun.generate(sprintf('tang_g%d.c',i-1));
    movefile(sprintf('tang_g%d.c',i-1),target_dir);
    % Generate rhs of dynamics
    if i < N+1
        dfk = jacobian(fk, [Xk; Uk]);
        d2fk = HESSk + jacobian(dfk.'*PIk, [Xk; Uk]);       
        fkfun = Function(sprintf('tang_f%d',i-1), {Xk, Uk, PIk, HESSk, P}, {d2fk, dfk(:,1:nx(i)), dfk(:,nx(i)+(1:nu(i))), fk});
        fkfun.generate(sprintf('tang_f%d.c',i-1));
        movefile(sprintf('tang_f%d.c',i-1),target_dir);
    end
    % Generate cost
    [d2jk, djk] = hessian(jk, [Xk; Uk]);
    d2jk = HESSk + d2jk;
    jkfun = Function(sprintf('tang_j%d',i-1), {Xk, Uk, HESSk, P}, {d2jk, djk(1:nx(i)), djk(nx(i)+(1:nu(i))), jk});
    jkfun.generate(sprintf('tang_j%d.c',i-1));
    movefile(sprintf('tang_j%d.c',i-1),target_dir);
end

% Collect bounds
nboffset = zeros(N+2,1);
ngoffset = zeros(N+2,1);
nbx = zeros(N+1);
nbu = zeros(N+1);
nb = zeros(N+1);
idxbdata = [];
ubdata = [];
lbdata = [];
ugdata = [];
lgdata = [];
% Bounds on inputs
idxbdata = [idxbdata;(0:dims.nu-1)'];
ubdata = [ubdata;dyn.ubu];
lbdata = [lbdata;dyn.lbu];
nbx(1) = 0;
nbu(1) = dims.nu;
nb(1) = nbx(1)+nbu(1);
nboffset(2) = length(ubdata);
tunnelbound_idx = zeros(Nta,1);
for i=2:N+1
    offsetx=0;
    offsetu=0;
    % Primary dynamics
    if i <= Nta
        idxbdata = [idxbdata;(0:dims.nu-1)'];
        ubdata = [ubdata;dyn.ubu];
        lbdata = [lbdata;dyn.lbu];
        nbx(i) = nbx(i) + 0;
        nbu(i) = nbu(i) + dims.nu;
        offsetx = dims.nu;
        offsetu = dims.nu;
        % Tunnel constraint
        tunnelbound_idx(i-1) = length(ugdata);
        ugdata = [ugdata; 0];
        lgdata = [lgdata; -1];
    elseif i == Nta+1
        % Tunnel constraint
        tunnelbound_idx(i-1) = length(ugdata);
        ugdata = [ugdata; 0];
        lgdata = [lgdata; -1];
    end
    % Secodary dynamics
    if i <= Ntr
        idxbdata = [idxbdata;offsetu+(0:dims.nu-1)'];
        ubdata = [ubdata; dyn.ubu];
        lbdata = [lbdata; dyn.lbu];
        nbx(i) = nbx(i) + 0;
        nbu(i) = nbu(i) + dims.nu;
    elseif i == Ntr+1
        lyapbound_idx = length(ugdata);
        % Lyapunov constraint (1) + terminal xi constraint (1) + terminal eta constraint (2)
        ugdata = [ugdata; 0; ocpopts.e; 0.001; 0.001];
        lgdata = [lgdata; -1; -1; -0.001; -0.001];
    end
    nboffset(i+1) = length(ubdata);
    ngoffset(i+1) = length(ugdata);
end
numb = length(idxbdata);
numg = length(ugdata);
nb = nbx + nbu;

% Generate header file with bounds and dimensions
dims_str = '#ifndef TANG_DIMS_BOUNDS_H_\n#define TANG_DIMS_BOUNDS_H_\n\n';
%dims
dims_str = [dims_str, sprintf('static inline int getN() {return %d;}\n', N)];
dims_str = [dims_str, sprintf('static inline int getNP() {return %d;}\n', size(P,1))];
dims_str = [dims_str, sprintf('static inline int getNtr() {return %d;}\n', Ntr)];
dims_str = [dims_str, sprintf('static inline int getNta() {return %d;}\n', Nta)];
dims_str = [dims_str, sprintf('static int nx [%d] = {', N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('%d, ', nx(i))];
end
dims_str = [dims_str, sprintf('%d};\n', nx(N+1))];
dims_str = [dims_str, sprintf('static int qpnx [%d] = {', N+1)];
dims_str = [dims_str, '0, '];
for i=2:N
    dims_str = [dims_str, sprintf('%d, ', nx(i))];
end
dims_str = [dims_str, sprintf('%d};\n', nx(N+1))];
dims_str = [dims_str, sprintf('static int nu [%d] = {', N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('%d, ', nu(i))];
end
dims_str = [dims_str, sprintf('%d};\n', nu(N+1))];
dims_str = [dims_str, sprintf('static int nbx [%d] = {', N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('%d, ', nbx(i))];
end
dims_str = [dims_str, sprintf('%d};\n', nbx(N+1))];
dims_str = [dims_str, sprintf('static int nbu [%d] = {', N+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nbu(i))];
end
dims_str = [dims_str, sprintf('%d};\n', nbu(N+1))];
dims_str = [dims_str, sprintf('static int nb [%d] = {', N+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', nbx(i)+nbu(i))];
end
dims_str = [dims_str, sprintf('%d};\n', nbx(N+1)+nbu(N+1))];
dims_str = [dims_str, sprintf('static int ng [%d] = {', N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('%d, ', ng(i))];
end
dims_str = [dims_str, sprintf('%d};\n', ng(N+1))];
dims_str = [dims_str, sprintf('static int ns [%d] = {', Ntr+1)];
for i=1:Ntr
    dims_str = [dims_str, sprintf('%d, ', 0)];
end
dims_str = [dims_str, sprintf('%d};\n', 0)];
%lyapbound_idx
dims_str = [dims_str, '\n\n'];
dims_str = [dims_str, sprintf('static int lyapbound_idx = %d;\n',lyapbound_idx)];
%tunnelbound_idx
dims_str = [dims_str, sprintf('static int tunnelbound_idx [%d] = {', Nta)];
for i=1:Nta-1
    dims_str = [dims_str, sprintf('%d, ', tunnelbound_idx(i))];
end
dims_str = [dims_str, sprintf('%d};\n', tunnelbound_idx(Nta))];
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
dims_str = [dims_str, sprintf('static int *idxb [%d] = {',N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('idxbdata+%d, ', nboffset(i))];
end
dims_str = [dims_str, sprintf('idxbdata+%d};\n', nboffset(N+1))];
dims_str = [dims_str, sprintf('static double *ub [%d] = {',N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('ubdata+%d, ', nboffset(i))];
end
dims_str = [dims_str, sprintf('ubdata+%d};\n', nboffset(N+1))];
dims_str = [dims_str, sprintf('static double *lb [%d] = {',N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('lbdata+%d, ', nboffset(i))];
end
dims_str = [dims_str, sprintf('lbdata+%d};\n', nboffset(N+1))];
dims_str = [dims_str, sprintf('static double *ug [%d] = {',N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('ugdata+%d, ', ngoffset(i))];
end
dims_str = [dims_str, sprintf('ugdata+%d};\n', ngoffset(N+1))];
dims_str = [dims_str, sprintf('static double *lg [%d] = {',N+1)];
for i=1:N
    dims_str = [dims_str, sprintf('lgdata+%d, ', ngoffset(i))];
end
dims_str = [dims_str, sprintf('lgdata+%d};\n', ngoffset(N+1))];
dims_str = [dims_str, '\n#endif // TANG_DIMS_BOUNDS_H_\n'];
fid = fopen([target_dir,'/tang_dims_bounds.h'],'wt');
fprintf(fid, dims_str);
fclose(fid);

% Generate source- and header-file with routine to get function-ptrs
%   source
funtypes = {{'tang_j',N+1},{'tang_g',N+1},{'tang_f',N}};
src_str = '#include "tang_functions.h"\n\n';
header_str = '#ifndef TANG_FUNCTIONS_H_\n#define TANG_FUNCTIONS_H_\n\n';
for ft = funtypes
    ft = ft{1};
    src_str = [src_str, '// Forward declarations\n'];
    for i=0:ft{2}-1
        src_str = [src_str, sprintf('int %s%d(const double** arg, double** res, int* iw, double* w, void* mem);\n',ft{1},i)];
        src_str = [src_str, sprintf('const int* %s%d_sparsity_in(int i);\n',ft{1},i)];
        src_str = [src_str, sprintf('const int* %s%d_sparsity_out(int i);\n',ft{1},i)];
        src_str = [src_str, sprintf('int %s%d_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);\n',ft{1},i)];
    end
    src_str = [src_str, '\n'];
    src_str = [src_str, '// Function picker\n'];
    header_str = [header_str, sprintf('int (*%sk(int stage))(const double **, double **, int *, double *, void *);\n',ft{1})];
    src_str = [src_str, sprintf('int (*%sk(int stage))(const double **, double **, int *, double *, void *) {\n',ft{1})];
    src_str = [src_str, 'int (*out)(const double **, double **, int *, double *, void *);\n'];
    src_str = [src_str, 'switch (stage) {\n'];
    for i=0:ft{2}-1
        src_str = [src_str, sprintf('case %d:\n', i)];
        src_str = [src_str, sprintf('out = &%s%d;\n',ft{1},i)];
        src_str = [src_str, sprintf('break;\n')];
    end
    src_str = [src_str, '}\n'];
    src_str = [src_str, 'return out;\n}\n\n'];
    header_str = [header_str, sprintf('const int *(*%sk_sparsity_in(int stage))(int);\n',ft{1})];
    src_str = [src_str, sprintf('const int *(*%sk_sparsity_in(int stage))(int) {\n',ft{1})];
    src_str = [src_str, 'const int *(*out)(int);\n'];
    src_str = [src_str, 'switch (stage) {\n'];
    for i=0:ft{2}-1
        src_str = [src_str, sprintf('case %d:\n', i)];
        src_str = [src_str, sprintf('out = &%s%d_sparsity_in;\n',ft{1},i)];
        src_str = [src_str, sprintf('break;\n')];
    end
    src_str = [src_str, '}\n'];
    src_str = [src_str, 'return out;\n}\n\n'];
    header_str = [header_str, sprintf('const int *(*%sk_sparsity_out(int stage))(int);\n',ft{1})];
    src_str = [src_str, sprintf('const int *(*%sk_sparsity_out(int stage))(int) {\n',ft{1})];
    src_str = [src_str, 'const int *(*out)(int);\n'];
    src_str = [src_str, 'switch (stage) {\n'];
    for i=0:ft{2}-1
        src_str = [src_str, sprintf('case %d:\n', i)];
        src_str = [src_str, sprintf('out = &%s%d_sparsity_out;\n',ft{1},i)];
        src_str = [src_str, sprintf('break;\n')];
    end
    src_str = [src_str, '}\n'];
    src_str = [src_str, 'return out;\n}\n\n'];
    header_str = [header_str, sprintf('int (*%sk_work(int stage))(int *, int *, int *, int *);\n',ft{1})];
    header_str = [header_str, '\n'];
    src_str = [src_str, sprintf('int (*%sk_work(int stage))(int *, int *, int *, int *) {\n',ft{1})];
    src_str = [src_str, 'int (*out)(int *, int *, int *, int *);\n'];
    src_str = [src_str, 'switch (stage) {\n'];
    for i=0:ft{2}-1
        src_str = [src_str, sprintf('case %d:\n', i)];
        src_str = [src_str, sprintf('out = &%s%d_work;\n',ft{1},i)];
        src_str = [src_str, sprintf('break;\n')];
    end
    src_str = [src_str, '}\n'];
    src_str = [src_str, 'return out;\n}\n\n'];
end
header_str = [header_str, '\n#endif // TANG_FUNCTIONS_H_\n'];

fid = fopen([target_dir,'/tang_functions.c'],'wt');
fprintf(fid, src_str);
fclose(fid);

fid = fopen([target_dir,'/tang_functions.h'],'wt');
fprintf(fid, header_str);
fclose(fid);

end