
%% Academic example


n = 2;  % number of states
p = 1;  % number of outputs
m = 1;  % number of inputs
r = p;  % number of disturbances

Q = 1*eye(p);
R = 1*eye(m);

% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

dT = 0.1;  % discrete time interval

A = [1 2; 3 2];
B2 = [0;1];
C2 = [1 0];

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

P = ss(A,B,C,D,dT);  % discrete time model -- transfer matrices

%% standard H2 control
[K,CL,gamma,info] = h2syn(P,p,m);  % there is some issue with this formulation! 2020/03/04


%% Closed-loop optimization
opts.N  = 5;
opts.solver  = 'mosek';
%    SLS
opts.type    = 1;
[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);

% IOP
opts.type    = 2;
[Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts);

%% time-domain simulaiton
close all

eps = 10^(-8);
Kiopr = minreal(Kiop,eps);
Kslsr = minreal(Ksls,eps);


G = ss(A,B2,C2,D22,dT);

T = 5;
Tn = floor(T/dT);   % number of iterations
dx = zeros(n,Tn);
dy = zeros(p,Tn);
du = zeros(m,Tn);

x0 = rand(n,1);  % initial disturbance
%dy(:,1) = 10*rand(p,1);  % initial measurement noise

%% SLS
[x,y,u,kesi] = dynsim(G,Ksls,dT,dx,dy,du,T,x0);
Time = (1:Tn)*dT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('SLS')
        subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1'); title('SLS')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');

% SLSr
[x,y,u,kesi] = dynsim(G,Kslsr,dT,dx,dy,du,T,x0);

Time = (1:Tn)*dT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('SLS reduction')
        subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1'); title('SLS reduction')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');

%% IOP
[x,y,u,kesi] = dynsim(G,Kiop,dT,dx,dy,du,T,x0);
Time = (1:Tn)*dT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('IOP')
        subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1'); title('IOP')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');

% IOPr
[x,y,u,kesi] = dynsim(G,Kiopr,dT,dx,dy,du,T,x0);
Time = (1:Tn)*dT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('IOP reduction')
        subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1'); title('IOP reduction')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');

%% H2 optimal controller  -- this behavior is very strange
% [x,y,u,kesi] = dynsim(G,K,dT,dx,dy,du,T,x0);
% Time = (1:Tn)*dT;
% figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('IOP reduction')
%            subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
% figure; 
% subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1'); title('H2 optimal controller')
% subplot(2,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');





