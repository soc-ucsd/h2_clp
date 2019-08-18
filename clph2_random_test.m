%% Comparsion with centralized H2, IOP, SLS in the standard case


clc; clear all; close all
clear yalmip
% dynamics in discrete time
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

% n = 5; m = 2; p = 2;

n = 2;  % number of states
p = 1;  % number of outputs
m = 1;  % number of inputs
r = p;  % number of disturbances
%q = 2*p;  % number of performance signals

Q = eye(p);
R = eye(m);

A   = randi([-15,15],n,n)/5; B1 = zeros(n,r); B2  = randi([-2,2],n,m);
C2  = randi([-2,2],p,n);  D21 = eye(p);     D22 = zeros(p,m);
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

P = ss(A,B,C,D,[]);  % discrete time model -- transfer matrices

%% standard H2 control
[K,CL,gamma,info] = h2syn(P,p,m);

%% IOP
opts.N       = 8;
opts.type    = 2;
opts.solver  = 'sedumi';
[Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts);

%% SLS
opts.type    = 1;
[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);

[gamma,H2iop,H2sls]