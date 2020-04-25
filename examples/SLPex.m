function SLPex
% Counter Example in the paper:
% Y. Zheng, L. Furieri, M. Kamgarpour, N, Li, Convex characterizations of
% stabilizing controllers: System-level, Input-output parameterizations and
% beyond
%
 
close all

load('slsex1.mat')   % load data
syms z;
T = 4;
n = size(A,1);

% --------------------------------------
%           Symbolic operations
% --------------------------------------

% step 1: check SLP consitriants  
Delta1 =  simplify((z*eye(n) - A)*R - B*M - eye(n)) % constraint 1
Delta2 = simplify((z*eye(n) - A)*N - B*L)          % constraint 2
Delta3 = simplify(R*(z*eye(n) - A) - N*C - eye(n)) % constraint 3
Delta4 = simplify(M*(z*eye(n) - A) - L*C)          % constraint 4

fprintf('\n Results from symbolic operations.\n')

% find hatDelta 

hDelta = simplify(Delta3 + R*(eye(n) + Delta1)^(-1)*(B*Delta4 - (z*eye(n) - A)*Delta3));
[~, d] = numden(det((eye(n) + hDelta)^(-1)));
hDeltaeig = roots(sym2poly(d))    % this is a subset of the closed-loop eigenvalues

% step 2: check controller & closed-loop    
K = simplify(L - M*R^(-1)*N);               % controller
CL = simplify((z*eye(n) - A - B*K*C)^(-1)); % closed-loop systems

% the following indentify should be zero
simplify( (eye(n) + hDelta)^(-1)*R - CL )

[~, d] = numden(det(CL));
Cleig = roots(sym2poly(d));    % this is a subset of the closed-loop eigenvalues

maxeig = max(abs(Cleig));
fprintf('The closed-loop eigenvalues using K = L - MR^(-1)N has a maximum norm %f .\n', maxeig)

% another controller realization K = L*(CN + I)^(-1)
K   = simplify(L*(C*N + eye(n))^(-1));                   % controller
CL2 = simplify((z*eye(n) - A - B*K*C)^(-1));             % closed-loop systems
[~, d] = numden(det(CL2));
maxeig = max(abs(roots(sym2poly(d))));
fprintf('The closed-loop eigenvalues using K = L(CN + I)^(-1) has a maximum norm %f .\n', maxeig)

% another controller K = M*R^(-1)
K   = simplify(M*R^(-1));                    % controller
CL1 = simplify((z*eye(n) - A - B*K*C)^(-1)); % closed-loop systems
[~, d] = numden(det(CL1));
maxeig = max(abs(roots(sym2poly(d))));
fprintf('The closed-loop eigenvalues using K = MR^(-1) has a maximum norm %f .\n %f .\n', maxeig)


% --------------------------------------
%      Convert into Transfer matrices
% --------------------------------------
fprintf('\n Results from transfer matrix operations.\n')

Rnum = simplify(R*z^(T+1),10);    % numerator of each component
Mnum = simplify(M*z^(T+1),10);
Nnum = simplify(N*z^(T+1),10);
Lnum = simplify(L*z^(T+1),10);

Rfir = subs(Rnum, z, 0);                       % FIR components
Mfir = subs(Mnum, z, 0);
Nfir = subs(Nnum, z, 0);
Lfir = subs(Lnum, z, 0);
for k = 1:T+1
    Rfir = [subs(diff(Rnum,k)/(factorial(k)), z, 0),Rfir];
    Mfir = [subs(diff(Mnum,k)/(factorial(k)), z, 0),Mfir];
    Nfir = [subs(diff(Nnum,k)/(factorial(k)), z, 0),Nfir];
    Lfir = [subs(diff(Lnum,k)/(factorial(k)), z, 0),Lfir];
end

Rfir = double(Rfir);
Mfir = double(Mfir);
Nfir = double(Nfir);
Lfir = double(Lfir);

Kstate  = Kslsreal(Lfir,Mfir,Rfir,Nfir,T+1);
G       = ss(A,B,C,[],[]);       % plant
CLstate = closedloop(G,Kstate);  % plant
Cleigs  = eig(CLstate.A);
maxeig = max(abs(Cleigs));
fprintf('The closed-loop eigenvalues using K = L - MR^(-1)N has a maximum norm %f .\n', maxeig);

% sort(Cleig)
% sort(Cleigs)

% another controller realization
Nfir = C*Nfir;
Nfir(:,1:n) = Nfir(1:n) + eye(n);
K1 = Kstatereal(Lfir,Nfir,T+1);  % state space realizaiton K = L(CN+I)^(-1)
CLstate1 = closedloop(G,K1);  % plant
Cl1eigs  = eig(CLstate1.A);
fprintf('The closed-loop eigenvalues using K = L(CN+I)^(-1) has a maximum norm %f .\n',  max(abs(Cl1eigs)));

% another controller realization K = MR^(-1)
Mfir = Mfir(:,n+1:end);
Rfir = Rfir(:,n+1:end);
K2 = Kstatereal(Mfir,Rfir,T);  % state space realizaiton K = (zM)(zR)^(-1)
CLstate2 = closedloop(G,K2);   % plant
Cl2eigs  = eig(CLstate2.A);
fprintf('The closed-loop eigenvalues using K = MR^(-1) has a maximum norm %f .\n',  max(abs(Cl2eigs)));


% -----------------------------------------------------
%      Time-domain simulation
% -----------------------------------------------------
dT = 1;
T  = 60;
Tn = floor(T/dT);   % number of iterations
dx = zeros(n,Tn);
dy = zeros(n,Tn);
du = zeros(n,Tn);

x0 = rand(n,1);    % initial disturbance

[x,~,~,~] = dynsim(G,Kstate,dT,dx,dy,du,T,x0);
Time = 1:Tn+1;
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('Iteration'); ylabel('x_1'); %title('IOP')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('Iteration'); ylabel('x_2');

[x,~,~,~] = dynsim(G,K1,dT,dx,dy,du,T,x0);
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('Iteration'); ylabel('x_1'); %title('IOP')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('Iteration'); ylabel('x_2');

[x,~,~,~] = dynsim(G,K2,dT,dx,dy,du,T,x0);
figure; 
subplot(2,1,1); plot(Time,x(1,:)); xlabel('Iteration'); ylabel('x_1'); %title('IOP')
subplot(2,1,2); plot(Time,x(2,:)); xlabel('Iteration'); ylabel('x_2');

end