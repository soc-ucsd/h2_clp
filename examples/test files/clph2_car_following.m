% Example 1: Car-following example


clc; clear all; 
clear yalmip

% dynamics in discrete time
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]


n = 3;  % number of states
p = 2;  % number of outputs
m = 1;  % number of inputs
r = p;  % number of disturbances

Q = eye(p);
R = eye(m);

tau    = 0.5;    % tau in the dynamics
deltaT = 0.1;    % discrete time interval

A   = eye(n) + [0 1 0;
                0 0 1;
                0 0 -1/tau]*deltaT; 
B2  = [0;0;1/tau]*deltaT;   
C2  = [1 0 0;
       0 1 0]; 

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

P = ss(A,B,C,D,deltaT);  % discrete time model -- transfer matrices

%% standard H2 control
[K,CL,gamma,info] = h2syn(P,p,m);

%% H2 synthesis via Closed-loop parameterization

N = 10:10:50;
%N = [N,75,100];

H2optimal = zeros(length(N),4);

for k = 1:length(N)

    opts.N       = N(k);
    opts.solver  = 'mosek';
%    SLS
    opts.type    = 1;
    [Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);

    % IOP
    opts.type    = 2;
    [Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts);
    
    % Type 3 -- mixed sls/iop
    opts.type    = 3;
    [Kty3,H2ty3,info3] = clph2(A,B2,C2,Q,R,opts);
    
    % Type 4 -- mixed sls/iop
    opts.type    = 4;
    [Kty4,H2ty4,info4] = clph2(A,B2,C2,Q,R,opts);
    
    H2optimal(k,:) = [H2iop, H2sls, H2ty3, H2ty4];
   %  H2optimal(k,:) = [0, 0, 0, H2ty4];

end

H2optimal