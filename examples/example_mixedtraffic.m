
%% 
clc; clear; close all;

% OVM parameters
% These parameters are the same as the paper:
% 
alpha  = 0.6;
beta   = 0.9;
v_max  = 30;
s_st   = 5;
s_go   = 35;

s_star = 20;

% system parameter after linearization  
alpha1 = alpha*v_max/2*pi/(s_go-s_st)*sin(pi*(s_star-s_st)/(s_go-s_st));
alpha2 = alpha+beta;
alpha3 = beta;

% dynamics
P1 = [0 -1;alpha1 -alpha2];
P2 = [0 1;0 alpha3];

A = [P1 zeros(2);P2 P1];
B = blkdiag([0;1],[0;1]);
C = [1 0 0 0;0 0 1 0];

% discretization 
dt = 0.1;
Ad = eye(4)+A*dt;
Bd = B*dt;
Cd = C;

n = 4;        % number of states 
p = 2;        % number of outputs 
m = 2;        % number of inputs 
Q = 1*eye(p); 
R = 1*eye(m);


% design controller
%% Standard H2 synthesis
% Generalized state-space model
% x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% z[t]   = C1x[t] + D11w[t] + D12 u[t]
% y[t]   = C2x[t] + D21w[t] + D22 u[t]

B1  = [Bd zeros(n,p)]; B2 = Bd;
D21 = [zeros(p,m) eye(p)]; D22 = zeros(p,m);
C2  = Cd;
C1  = [Q^(0.5)*C2;
       zeros(m,n)];
D11 = [Q^(0.5)*D21;
       zeros(m,p+m)];
D12 = [zeros(p,m);
        R^(0.5)];
%  generalized state-space model
Bg = [B1, B2];
Cg = [C1;C2];
Dg = [D11 D12;
      D21 D22];

P = ss(Ad,Bg,Cg,Dg,dt);            % discrete time model -- transfer matrices
[K,CL,gamma,info] = h2syn(P,p,m);  % this norm is not consistent with the results from CLP. 

G = ss(Ad,Bd,Cd,[],dt);
sqrt(norm((eye(2)-G*K)^(-1),2)^2 + norm((eye(2)-G*K)^(-1)*G,2)^2 + ...
    norm(K*(eye(2)-G*K)^(-1),2)^2 + norm((eye(2)-K*G)^(-1),2)^2)    % this cost is consistent with our computation


% %% H2 synthesis via Closed-loop parameterization
% 
% 
% N = [10,15,20];%%[10,15, 20, 25, 30, 50, 75, 100];
% 
% H2optimal = zeros(length(N),4);
% Ksls = cell(length(N),1);
% Kiop = cell(length(N),1);
% Kty3 = cell(length(N),1);
% Kty4 = cell(length(N),1);
% for k = 1:length(N)
%     clear Yalmip
%     
%     opts.N        = N(k);
%     opts.solver   = 'mosek';
%     opts.costType = 1;
%     
%     opts.type    = 1;
%     [Ksls{k},H2sls,infosls] = clph2(Ad,Bd,Cd,Q,R,opts);
% 
%     % IOP
%     opts.type    = 2;
%     [Kiop{k},H2iop,infoiop] = clph2(Ad,Bd,Cd,Q,R,opts);
%     
%     % Type 3 -- mixed sls/iop
%      opts.type    = 3;
%      [Kty3{k},H2ty3,info3] = clph2(Ad,Bd,Cd,Q,R,opts);
%     
%     % Type 4 -- mixed sls/iop
%     opts.type    = 4;
%     [Kty4{k},H2ty4,info4] = clph2(Ad,Bd,Cd,Q,R,opts);
%     
%     H2optimal(k,:) = [H2iop, H2sls, H2ty3, H2ty4];
% end

load Ex_mixedtraffic

%% Time domain simulation
    T = 10;
    Tn = floor(T/dt);   % number of iterations
    dx = zeros(n,Tn);
    dy = zeros(p,Tn);
    du = zeros(m,Tn);
    %x0 = 5*rand(n,1);  % initial disturbance
    x0 = [3,0,-2,0]';
    G = ss(Ad,Bd,Cd,[],dt);
    Time = (0:Tn)*dt;

%% SLS
    [x,y,u,kesi_sls] = dynsim(G,Ksls{6},dt,dx,dy,du,T,x0);    
    figure; h1 = plot(Time,u(1,:),'b','linewidth',1.2); hold on
            h2 = plot(Time,u(2,:),'m','linewidth',1.2); 
            xlabel('Time (s)','Interpreter','latex','FontSize',10); 
            ylabel('Control input','Interpreter','latex','FontSize',10);
            h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
            set(h,'FontSize',10,'Interpreter','latex','box','off')
            set(gcf,'Position',[250 150 400 200]);
    figure; 
    subplot(2,1,1); h1 = plot(Time,x(1,:),'b','linewidth',1.2); hold on;
                    h2 = plot(Time,x(3,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',10);
                    ylabel('Spacing error (m)','Interpreter','latex','FontSize',10);
                    h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
                    set(h,'FontSize',10,'Interpreter','latex','box','off')
    subplot(2,1,2); h1 = plot(Time,x(2,:),'b','linewidth',1.2); hold on
                    h2 = plot(Time,x(4,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',10);
                    ylabel('Velocity error (m/s)','Interpreter','latex','FontSize',10);
                    h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
                    set(h,'FontSize',10,'Interpreter','latex','box','off')
    set(gcf,'Position',[250 150 400 350]);
%print(gcf,['../../Figures/timeBlockArrow'],'-painters','-depsc2','-r600')

    
%% IOP
  [x,y,u,kesi] = dynsim(G,Kiop{6},dt,dx,dy,du,T,x0);
    figure; h1 = plot(Time,u(1,:),'b','linewidth',1.2); hold on
            h2 = plot(Time,u(2,:),'m','linewidth',1.2); 
            xlabel('Time (s)','Interpreter','latex','FontSize',10); 
            ylabel('Control input','Interpreter','latex','FontSize',10);
            h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
            set(h,'FontSize',10,'Interpreter','latex','box','off')
            set(gcf,'Position',[250 150 400 200]);
    figure; 
    subplot(2,1,1); h1 = plot(Time,x(1,:),'b','linewidth',1.2); hold on;
                    h2 = plot(Time,x(3,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',10);
                    ylabel('Spacing error (m)','Interpreter','latex','FontSize',10);
                    h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
                    set(h,'FontSize',10,'Interpreter','latex','box','off')
    subplot(2,1,2); h1 = plot(Time,x(2,:),'b','linewidth',1.2); hold on
                    h2 = plot(Time,x(4,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',10);
                    ylabel('Velocity error (m/s)','Interpreter','latex','FontSize',10);
                    h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
                    set(h,'FontSize',10,'Interpreter','latex','box','off')
    set(gcf,'Position',[250 150 400 350]);
