function [K,H2,info] = clph2(A,B,C,Q,R,userOpts,T_struct,R_struct)
%
% close-loop parameterization of H2 optimal control for discrete-time systems
%
% Input variables
%      (A,B,C):    system dynamics in discrete time
%      Q:    performance weights on output y
%      R:    performance weights on input u
%
% userOpts is a structure and contains the following options
%      N:     Oder of FIR approximation (default:8)
%      solver: sedumi, sdpt3, csdp or mosek (default)

% The plant dynamics are:
%                           x = Ax_t + Bu_t
%                           y = Cx_t
% The orginal problem is as follows
%               min_{K} max_{G} ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
%                  subject to     K internally stabilizes G
%                                 K \in S
% where Y, U denote the closed-loop transfer matrices from w to y and u,
% respectively.
%
% We solve the problem using closed-loop parameterizations, one of them is
% as follows
%
%   min_{beta \in [0,1/gamma)}  min_{Y,U,W,Z}  ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
%               subjec to      [I -G][Y W]
%                                    [U Z]  = [I 0]   (1)
%                                 [Y W][-G] = [I]
%                                 [U Z][I]    [0]     (2)
%                              Y,U,W,Z \in FIR(N)     (3)


% Rely on YALMIP to reformulate the above problem into an SDP, then call
% Mosek/SeDuMi to get a soluton
%
% Authors: Y. Zheng & L. Furieri

% ============================================
%            set uers' options
% ============================================
opts = clpOpts;
if nargin < 5
    error('Robustiop requires three inputs: A, B, C, Q, R; check your inputs.')
elseif nargin > 3
    opts = setUserOpts(opts,userOpts);
end
N     = opts.N;       % FIR horizon
Type  = opts.type;    % which parameterization is used


% ========================================================================
%           Dynamics in symbolic form
% ========================================================================
% To do: we need to avoid this symbolic operation
z  = sym('z');
n  = size(A,1);
G  = C*(z*eye(n)- A)^(-1)*B;  % symbolic form
[p,m] = size(G);      % system dimension

% ========================================================================
%           Controller design using IOP
% ========================================================================
fprintf('==========================================================\n')
fprintf('       H2 synthesis via closed-loop paramterization       \n')
fprintf('==========================================================\n')
fprintf('Number of inputs       : %i\n',m);
fprintf('Number of outputs      : %i\n',p);
fprintf('FIR horizon            : %i\n',N);
fprintf('Parameterization No.   : %i\n',Type);
fprintf('Numerical solver       : %s\n',opts.solver);
fprintf('==========================================================\n');

% ========================================================================
% Start CLP: closed-loop parameterization
% ========================================================================

switch Type
    case 1    % sls
        % Define FIR variables (3)
        CRv = sdpvar(n,n*(N+1));          % decision variables for R
        CMv = sdpvar(m,n*(N+1));          % decision variables for M
        CNv = sdpvar(n,p*(N+1));          % decision variables for N
        CLv = sdpvar(m,p*(N+1));          % decision variables for L
        
        % achivability constraint (1)-(3)
        fprintf('Step 1: Encoding the achievability constraints ...\n')
        const = slscons(A,B,C,CRv,CMv,CNv,CLv,N);
        sls_sparsity_cons
        %const=[const, const_sparsity];
        
        
        % H2 cost
        fprintf('Step 2: Encoding the H2 cost ...\n')
        cost   = slscost(CNv,CLv,Q,R,C,N);
        
    case 2    % iop
        % decision variables
        CYv = sdpvar(p,p*(N+1));          % decision variables for Y
        CUv = sdpvar(m,p*(N+1));          % decision variables for U
        CWv = sdpvar(p,m*(N+1));          % decision variables for W
        CZv = sdpvar(m,m*(N+1));          % decision variables for Z
        
        % achivability constraint (1)-(3)
        fprintf('Step 1: Encoding the achievability constraints ...\n')
        const = iopcons(G,CYv,CUv,CWv,CZv,N,z);
        iop_sparsity_cons
        
        % H2 cost
        fprintf('Step 2: Encoding the H2 cost ...\n')
        cost   = iopcost(CYv,CUv,Q,R,N);
        
    case 3
        YXv = sdpvar(p,n*(N+1));          % decision variables for Y
        YYv = sdpvar(p,p*(N+1));          % decision variables for U
        UXv = sdpvar(m,n*(N+1));          % decision variables for W
        UYv = sdpvar(m,p*(N+1));          % decision variables for Z
        
        % achivability constraint (1)-(3)
        fprintf('Step 1: Encoding the achievability constraints ...\n')
        const = par3cons(G,YXv,YYv,UXv,UYv,N,z,A,C);
        par3_sparsity_cons
        
        % H2 cost
        fprintf('Step 2: Encoding the H2 cost ...\n')
        cost   = par3cost(YYv,UYv,Q,R,N);
        
    case 4
        XYv = sdpvar(n,p*(N+1));          % decision variables for Y
        XUv = sdpvar(n,m*(N+1));          % decision variables for U
        UYv = sdpvar(m,p*(N+1));          % decision variables for W
        UUv = sdpvar(m,m*(N+1));          % decision variables for Z
        
        % achivability constraint (1)-(3)
        fprintf('Step 1: Encoding the achievability constraints ...\n')
        const = par4cons(G,XYv,XUv,UYv,UUv,N,z,A,B);
        par4_sparsity_cons
        
        % H2 cost
        fprintf('Step 2: Encoding the H2 cost ...\n')
        cost   = par4cost(XYv,UYv,Q,R,N,C,n);
        
    otherwise
        error('Unknown parameterizaiton')
end


% Get a solution
fprintf('Step 3: call a solver to get a solution...\n\n');
yalmipOpts = sdpsettings('solver',opts.solver,'verbose',opts.verbose);
sol        = optimize(const,cost,yalmipOpts);


if sol.problem == 0
    fprintf('\n Problem is solved, and an h2 stabilizing controller is found\n')
    H2 = sqrt(value(cost));  % value of the H2 norm!
    fprintf(' H2 norm of the normal closed-loop system is %6.4f \n\n', H2);
else
    fprintf('\n Numerical issue ... the problem may be not feasible ...\n\n')
    H2 = sqrt(value(cost));
end


% ========================================================================
%                            Post process
% ========================================================================

info.H2 = H2;  %  H2 norm of the normial system

switch Type
    case 1    % sls
        info.var.R = value(CRv);
        info.var.M = value(CMv);
        info.var.N = value(CNv);
        info.var.L = value(CLv);
        temp = C*info.var.N;
        temp(:,1:p) = temp(1:p) + eye(p);
        K = Kstatereal(info.var.L,temp,N);  % state space realizaiton K = L(CN+I)^(-1)
        
        %  closed-loop responses from optimization
        z = tf('z');
        Rt = zeros(n,n); Mt = zeros(m,n); Nt = zeros(n,p); Lt = zeros(m,p);
        for k = 1:N+1
            Rt = Rt + info.var.R(:,n*(k-1)+1:n*k)*z^(1-k);  % FIR
            Mt = Mt + info.var.M(:,n*(k-1)+1:n*k)*z^(1-k);  % FIR
            Nt = Nt + info.var.N(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
            Lt = Lt + info.var.L(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
        end
        info.cl.R = Rt;
        info.cl.M = Mt;
        info.cl.N = Nt;
        info.cl.L = Lt;
        
    case 2    % iop
        
        info.var.Y  = round(value(CYv),opts.eps);
        info.var.U  = round(value(CUv),opts.eps);
        info.var.W  = round(value(CWv),opts.eps);
        info.var.Z  = round(value(CZv),opts.eps);
        K           = Kstatereal(info.var.U,info.var.Y,N);  % state-space realization
        info.Ks     = K;
        %  closed-loop responses from optimization
        z = tf('z');
        Y = zeros(p,p); U = zeros(m,p);
        for k = 1:N+1
            Y = Y + info.var.Y(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
            U = U + info.var.U(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
        end
        info.cl.Y = Y;
        info.cl.U = U;
        
    case 3
        disp('to do')
        K=0;
        
    case 4
        disp('to do')
        K=0;
        
    otherwise
        error('Unknown parameterizaiton')
end


end
