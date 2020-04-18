
%% Inverted pendulum example
% took from: http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling

% state: 'x' 'x_dot' 'phi' 'phi_dot' 

M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];


n = 4;  % number of states
p = 2;  % number of outputs
m = 1;  % number of inputs
r = p;  % number of disturbances

Q = 1*eye(p);
R = 1*eye(m);


%% discrete-time model
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

dT = 0.1;  % discrete time interval

A = eye(n) + A*dT;
B2 = B*dT;
C2 = C;

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
[K,CL,gamma,info] = h2syn(P,p,m);  % there is some issue with this formulation! 2020/03/04


%% IOP
opts.N  = 40;
opts.solver  = 'mosek';
%    SLS
opts.type    = 1;
[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);

% IOP
opts.type    = 2;
[Kiop,H2iop,infoiop] = clph2(A,B2,C2,Q,R,opts);

%% time-domain simulaiton
eps = 10^(-6);
Kiopr = minreal(Kiop,eps);
Kslsr = minreal(Ksls,eps);


G = ss(A,B2,C2,D22,deltaT);

T = 10;
Tn = floor(T/deltaT);   % number of iterations
dx = zeros(n,Tn);
dy = zeros(p,Tn);
du = zeros(m,Tn);

x0 = [0;0;0.2;0];%0.1*rand(n,1);  % initial disturbance
%dy(:,1) = 10*rand(p,1);  % initial measurement noise
close all
%% SLS
[x,y,u,kesi] = dynsim(G,Ksls,deltaT,dx,dy,du,T,x0);
Time = (1:Tn)*deltaT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('SLS')
           subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure;
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('Position');  title('SLS')
subplot(2,1,2); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('Angle');

% SLS reduction
[x,y,u,kesi] = dynsim(G,Kslsr,deltaT,dx,dy,du,T,x0);
Time = (1:Tn)*deltaT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('SLS reduction')
        subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure;
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('Position');  title('SLS reduction')
subplot(2,1,2); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('Angle');

%% IOP
[x,y,u,kesi] = dynsim(G,Kiop,deltaT,dx,dy,du,T,x0);
Time = (1:Tn)*deltaT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('IOP')
           subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure; title('IOP')
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('Position');
subplot(2,1,2); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('Angle');

% SLS reduction
[x,y,u,kesi] = dynsim(G,Kiopr,deltaT,dx,dy,du,T,x0);
Time = (1:Tn)*deltaT;
figure; subplot(2,1,1); plot(Time,u); xlabel('time (s)'); ylabel('Control input'); title('IOP reduction')
        subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure;
subplot(2,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('Position'); title('IOP reduction')
subplot(2,1,2); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('Angle');



