
%% Figure 3 in the paper
close all;
load Ex_mixedtraffic

fontsize = 10;

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
    [x,y,u,kesi_sls] = dynsim(G,Ksls{5},dt,dx,dy,du,T,x0);    
    figure;
    subplot(3,1,1);
    h1 = plot(Time,u(1,:),'b','linewidth',1.2); hold on
            h2 = plot(Time,u(2,:),'m','linewidth',1.2); 
            xlabel('Time (s)','Interpreter','latex','FontSize',fontsize); 
            ylabel('Control input','Interpreter','latex','FontSize',fontsize);
    subplot(3,1,2); h1 = plot(Time,x(1,:),'b','linewidth',1.2); hold on;
                    h2 = plot(Time,x(3,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',fontsize);
                    ylabel('Spacing error (m)','Interpreter','latex','FontSize',fontsize);
    subplot(3,1,3); h1 = plot(Time,x(2,:),'b','linewidth',1.2); hold on
                    h2 = plot(Time,x(4,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',fontsize);
                    ylabel('Velocity error (m/s)','Interpreter','latex','FontSize',fontsize);
                 
   h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','south'); 
  set(h,'orientation','horizontal','FontSize',10,...
      'box','off','Position',[0.2 0.02 0.7 0.03],'Interpreter','latex');   
  set(gca,'TickLabelInterpreter','latex');
  set(gcf,'Position',[250 150 400 400]);
  print(gcf,'FIR30','-painters','-depsc2','-r600')
    
%% IOP
  [x,y,u,kesi_iop] = dynsim(G,Kiop{5},dt,dx,dy,du,T,x0);
  %[x,y,u,kesi] = dynsim(G,Kiop{6},dt,dx,dy,du,T,x0);
  
  
  %% FIR 50
  [x,y,u,kesi_sls] = dynsim(G,Ksls{6},dt,dx,dy,du,T,x0);    
    figure;
    subplot(3,1,1);
    h1 = plot(Time,u(1,:),'b','linewidth',1.2); hold on
            h2 = plot(Time,u(2,:),'m','linewidth',1.2); 
            xlabel('Time (s)','Interpreter','latex','FontSize',fontsize); 
            ylabel('Control input','Interpreter','latex','FontSize',fontsize);
    subplot(3,1,2); h1 = plot(Time,x(1,:),'b','linewidth',1.2); hold on;
                    h2 = plot(Time,x(3,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',fontsize);
                    ylabel('Spacing error (m)','Interpreter','latex','FontSize',fontsize);
                   % h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
                   % set(h,'FontSize',10,'Interpreter','latex','box','off')
    subplot(3,1,3); h1 = plot(Time,x(2,:),'b','linewidth',1.2); hold on
                    h2 = plot(Time,x(4,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',fontsize);
                    ylabel('Velocity error (m/s)','Interpreter','latex','FontSize',fontsize);
                 
   h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','south'); 
  set(h,'orientation','horizontal','FontSize',10,...
      'box','off','Position',[0.2 0.02 0.7 0.03],'Interpreter','latex');
  set(gca,'TickLabelInterpreter','latex');
  set(gcf,'Position',[250 150 400 400]);
  print(gcf,'FIR50','-painters','-depsc2','-r600')
  
  %% FIR 75
  [x,y,u,kesi_sls] = dynsim(G,Ksls{7},dt,dx,dy,du,T,x0);    
    figure;
    subplot(3,1,1);
    h1 = plot(Time,u(1,:),'b','linewidth',1.2); hold on
            h2 = plot(Time,u(2,:),'m','linewidth',1.2); 
            xlabel('Time (s)','Interpreter','latex','FontSize',fontsize); 
            ylabel('Control input','Interpreter','latex','FontSize',fontsize);
    subplot(3,1,2); h1 = plot(Time,x(1,:),'b','linewidth',1.2); hold on;
                    h2 = plot(Time,x(3,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',fontsize);
                    ylabel('Spacing error (m)','Interpreter','latex','FontSize',fontsize);
                   % h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','Northeast');
                   % set(h,'FontSize',10,'Interpreter','latex','box','off')
    subplot(3,1,3); h1 = plot(Time,x(2,:),'b','linewidth',1.2); hold on
                    h2 = plot(Time,x(4,:),'m','linewidth',1.2); 
                    xlabel('Time (s)','Interpreter','latex','FontSize',fontsize);
                    ylabel('Velocity error (m/s)','Interpreter','latex','FontSize',fontsize);
                 
   h = legend([h1,h2],'Vehicle 1','Vehicle 2','Location','south'); 
  set(h,'orientation','horizontal','FontSize',10,...
      'box','off','Position',[0.2 0.02 0.7 0.03],'Interpreter','latex');          
  set(gcf,'Position',[250 150 400 400]);
  set(gca,'TickLabelInterpreter','latex');
  print(gcf,'FIR75','-painters','-depsc2','-r600')
