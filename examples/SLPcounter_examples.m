
%% Counter examples for the SLS output feedback controller K = L - MR^(-1)N


clc;clear;close all
load('Counter_examples\data_mosek_failure_with_cost_case6.mat')

opts.solver   = 'mosek';  % high precision
opts.costType = 1;
opts.type     = 1; % SLP
[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);


%% time-domain simulaiton
T = 5;
deltaT = 0.2;
Tn = floor(T/deltaT);   % number of iterations
dx = zeros(n,Tn);
dy = zeros(p,Tn);
du = zeros(m,Tn);
Time = (0:Tn)*deltaT;

x0 = 5*rand(n,1);  % initial disturbance

[x,y,u,kesi] = dynsim(G,Ksls,deltaT,dx,dy,du,T,x0);
figure; subplot(2,1,1); plot(Time,u(1,:)); xlabel('time (s)'); ylabel('Control input'); title('SLS')
           subplot(2,1,2); plot(Time,max(abs(kesi))); xlabel('time (s)'); ylabel('internal state');
figure; 
subplot(3,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1');  title('SLS')
subplot(3,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');
subplot(3,1,3); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('x_3');


