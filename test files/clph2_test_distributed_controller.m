%% Comparsion with distributed H2 synthesis using different closed-loop parameterization
% IOP, SLS in the standard case


clc; clear all; close all
clear yalmip
% dynamics in discrete time
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]


n = 3;  % number of states
p = 3;  % number of outputs
m = 3;  % number of inputs
r = p;  % number of disturbances

Q = eye(p);
R = eye(m);

%%ASSIGN A,B,C by hand
A = [0.8000         0         0;
   -0.2000    0.6000         0;
   -1.0000    0.4000    1.6000];
B2 = diag([3,-5,-4]);
C2 = diag([1 -4 -4]);

B1 = zeros(n,r); 
D21 = eye(p);     D22 = zeros(p,m);
C1  = [Q^(0.5)*C2;
    zeros(m,n)];
D11 = [Q^(0.5);
    zeros(m,p)];
D12 = [zeros(p,m);
    R^(0.5)];
     

%  generalized state-space model
B = [B1, B2];
C = [C1;C2];
D = [D11 D12;
    D21 D22];


disp('*****Distributed case: comparing the 4 parametrizations\n')
T_struct = tril(ones(m,p));   % lower triangular controller
R_struct = tril(ones(p,p));


N = 4:2:10;    % FIR horizon 


H2optimal_distributed = zeros(length(N),4);
opts.solver  = 'sedumi';
opts.spa     = 1;
opts.T       = T_struct;     % add a SI procedure for these two patterns latter
opts.R       = R_struct;
    

for k = 1:length(N)
    
    opts.N       = N(k);

    % SLS
    opts.type    = 1;
    [Ksls,H2sls,infosls]    = clph2(A,B2,C2,Q,R,opts);
    
    % IOP
    opts.type    = 2;
    [Kiop,H2iop,infoiop]    = clph2(A,B2,C2,Q,R,opts);
    
    % Type 3
    opts.type    = 3;
    [Kpar3,H2par3,infopar3] = clph2(A,B2,C2,Q,R,opts); %%Kpar3 not assigned (couldn't figure out the code) 
    
    
    % Type 4
    opts.type    = 4;
    [Kpar4,H2par4,infopar4] = clph2(A,B2,C2,Q,R,opts); %%Kpar4 not assigned (couldn't figure out the code) 
    
    H2optimal_distributed(k,1:4) = [H2sls,H2iop,H2par3,H2par4];
    
end

H2optimal_distributed

