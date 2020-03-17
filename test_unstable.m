%% open-loop stable systems

clc; clear all; close all
clear yalmip
% dynamics in discrete time
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

% n = 5; m = 2; p = 2;

n = 6;  % number of states
p = 4;  % number of outputs
m = 4;  % number of inputs
r = p;  % number of disturbances
%q = 2*p;  % number of performance signals

Q = eye(p);
R = eye(m);

A   = randi([-5,5],n,n); B2  = randi([-10,10],n,m)/10;
C2  = randi([-10,10],p,n)/10;  

eig(A)

%% SLS

%load data_truefailure4

opts.type    = 1;
opts.solver = 'sedumi';
%opts.solver = 'mosek';
opts.costType = 3;    % no penalty on the cost
opts.eps = 1e-3;      % low solution precision, try 1e-3, this one can lead to unstable controller
                      % 1e-4, 1e-6
opts.N   = 10;
[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);

%% Simulation --
if infosls.problem == 0 % only do time-domain simulation when the solver reports no numerical errors
    T = 5;
    deltaT = 0.2;
    Tn = floor(T/deltaT);   % number of iterations
    dx = zeros(n,Tn);
    dy = zeros(p,Tn);
    du = zeros(m,Tn);

    x0 = 5*rand(n,1);  % initial disturbance


    G = ss(A,B2,C2,[],deltaT);

    Time = (1:Tn)*deltaT;


    [x,y,u,kesi] = dynsim(G,Ksls,deltaT,dx,dy,du,T,x0);
    figure; subplot(2,1,1); plot(Time,u(1,:)); xlabel('time (s)'); ylabel('Control input'); title('SLS')
               subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
    figure; 
    subplot(3,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1');  title('SLS')
    subplot(3,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');
    subplot(3,1,3); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('x_3');

end
