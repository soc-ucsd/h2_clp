%% open-loop stable systems

clc; clear all; close all
clear yalmip
% dynamics in discrete time
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

% n = 5; m = 2; p = 2;
Num = 1;
    
while true
    
    
n = 3;  % number of states
p = 1;  % number of outputs
m = 1;  % number of inputs
r = p;  % number of disturbances
%q = 2*p;  % number of performance signals

Q = eye(p);
R = eye(m);

A   = randi([-5,5],n,n); B2  = randi([-1,1],n,m);
C2  = randi([-1,1],p,n);  

%eig(A)

%% SLS


opts.type    = 1;
%opts.solver = 'sedumi';
opts.solver = 'mosek';
opts.costType = 2;    % 2: penalty for all Y U W Z; 3: no penalty on the cost

opts.N   = 10;
[Ksls,H2sls,infosls] = clph2(A,B2,C2,Q,R,opts);
G = ss(A,B2,C2,[],[]);
CL = closedloop(G,Ksls);

fprintf('Number of trials       : %i\n',Num);
fprintf('Maximum eigenvalue     : %4.2f\n',max(abs(eig(CL.A))));

if infosls.problem == 0 && max(abs(eig(CL.A))) > 1  % test sedumi
    break;
end
Num = Num + 1;

end

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
