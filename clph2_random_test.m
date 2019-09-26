%% Comparsion with centralized H2, IOP, SLS in the standard case


clc; clear all; close all
clear yalmip
% dynamics in discrete time
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

% n = 5; m = 2; p = 2;

n = 3;  % number of states
p = 3;  % number of outputs
m = 3;  % number of inputs
r = p;  % number of disturbances
%q = 2*p;  % number of performance signals

Q = eye(p);
R = eye(m);



A   = tril(randi([-15,15],n,n)/5); B1 = zeros(n,r); B2  = diag(diag(randi([-5,5],n,m))); %A lower triangular, B,C diagonal
C2  = diag(diag(randi([-5,5],p,n)));  D21 = eye(p);     D22 = zeros(p,m);
C1  = [Q^(0.5)*C2;
    zeros(m,n)];
D11 = [Q^(0.5);
    zeros(m,p)];
D12 = [zeros(p,m);
    R^(0.5)];

%%ASSIGN A,B,C by hand
A = [0.8000         0         0;
   -0.2000    0.6000         0;
   -1.0000    0.4000    1.6000];
B2 = diag([3,-5,-4]);
C2 = diag([1 -4 -4]);
%A =[3.0000         0;
%   -0.4000   -1.4000];

%B2 = [-1     0;
% 0     1];


%C2 = [ 2     0;
%     0    -2];


     

%  generalized state-space model
B = [B1, B2];
C = [C1;C2];
D = [D11 D12;
    D21 D22];



P = ss(A,B,C,D,[]);  % discrete time model -- transfer matrices

%% standard H2 control
%[K,CL,gamma,info] = h2syn(P,p,m);

%% IOP
opts.N       = n;
opts.type    = 2;
opts.solver  = 'sedumi';
[Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts);

%% SLS
%opts.type    = 1;
%[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct);

%[H2iop,H2sls]

disp('*****Centralized case: comparing the 4 parametrizations\n')
T_struct=ones(m,p);
R_struct=ones(p,p);

N = 2:2:10;
for k = 1:length(N)
    
    opts.N       = N(k);
    opts.solver  = 'mosek';
    % SLS
    opts.type    = 1;
    [Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct);
    
    % IOP
    opts.type    = 2;
    [Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct);
    
    %Type 3
    opts.type    = 3;
    [Kpar3,H2par3,infopar3] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct); %%Kpar3 not assigned (couldn't figure out the code) 
    
    
    %Type 4
    opts.type    = 4;
    [Kpar4,H2par4,infopar4] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct); %%Kpar4 not assigned (couldn't figure out the code) 
    
    H2optimal_centralized(k,1:4) = [H2sls,H2iop,H2par3,H2par4];
    
end

disp('*****Distributed case: comparing the 4 parametrizations\n')
T_struct=tril(ones(m,p));
R_struct=tril(ones(p,p));


N = 2 :2:10;
for k = 1:length(N)
    
    opts.N       = N(k);
    opts.solver  = 'mosek';
   % SLS
    opts.type    = 1;
    [Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct);
    
    % IOP
    opts.type    = 2;
    [Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct);
    
    %Type 3
    opts.type    = 3;
    [Kpar3,H2par3,infopar3] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct); %%Kpar3 not assigned (couldn't figure out the code) 
    
    
    %Type 4
    opts.type    = 4;
    [Kpar4,H2par4,infopar4] = clph2(A,B2,C2,Q,R,opts,T_struct,R_struct); %%Kpar4 not assigned (couldn't figure out the code) 
    
    H2optimal_distributed(k,1:4) = [H2sls,H2iop,H2par3,H2par4];
    
end

H2optimal_centralized
H2optimal_distributed




%tf controllers
%z=sym('z');
%Kiop_tf=vpa(simplifyFraction(Kiop.C*inv(z*eye(size(Kiop.A,1))-Kiop.A)*Kiop.B+Kiop.D),3);
%Ksls_tf=vpa(simplifyFraction(Ksls.C*inv(z*eye(size(Ksls.A,1))-Ksls.A)*Ksls.B+Ksls.D),3);