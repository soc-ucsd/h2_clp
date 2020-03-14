function [K,h2,info] = iop(A,B,C,Q,R,K0,userOpts)
% A, B, C: system dynamics
% Q, R perform index
% K0, an initial stabilizing controller  -- assume to be static for
% simplicity

% ============================================
%            set uers' options
% ============================================
    opts = clpOpts;
    if nargin < 6
        error('Robustiop requires three inputs: A, B, C, Q, R, K0; check your inputs.')
    elseif nargin > 6
        opts = setUserOpts(opts,userOpts);
    end

    
    [n,m] = size(B);
    [p,~] = size(C);

    G     = ss(A,B,C,[],[]);
    %% state-space realization of T1, T2, T3
    T1 = ss(A + B*K0*C,B*K0,C,eye(p),[]);
    T2 = ss(A + B*K0*C,B,C,zeros(p,m),[]);
    T3 = ss(A + B*K0*C,B,K0*C,eye(m),[]);

    % minreal(T1 - (eye(p) - G*K0)^(-1))
    % minreal(T2 - (eye(p) - G*K0)^(-1)*G)
    % minreal(T3 - (eye(m) - K0*G)^(-1))


    %% IOP formulation 
    N     = opts.N;       % FIR horizon
  % decision variables
    CYv = sdpvar(p,p*(N+1));          % decision variables for Y
    CUv = sdpvar(m,p*(N+1));          % decision variables for U
    CWv = sdpvar(p,m*(N+1));          % decision variables for W
    CZv = sdpvar(m,m*(N+1));          % decision variables for Z

    % To do: we need to avoid this symbolic operation
    z  = sym('z');
    hG  = T2.C*(z*eye(size(T2.A,1))- T2.A)^(-1)*T2.B;  % symbolic form

    % achivability constraint (1)-(3)
    Step = 1;
    fprintf('Step %d: Encoding the achievability constraints ...\n',Step)
    const = iopcons(hG,CYv,CUv,CWv,CZv,N,z,opts);
    Step  = Step +1;
   
    % H2 cost
    fprintf('Step %d: Encoding the H2 cost ...\n',Step)
    [cost,const] = iopcost2(const,CYv,CUv,CWv,CZv,Q,R,T1,T2,T3,K0,N,opts);
    Step = Step +1;

    % Get a solution
    fprintf('Step %d: call a solver to get a solution...\n\n',Step);
    yalmipOpts = sdpsettings('solver',opts.solver,'verbose',opts.verbose);
    sol        = optimize(const,cost,yalmipOpts);



    % ========================================================================
    %                            Post processing
    % ========================================================================
    h2 = sqrt(value(cost));
    info.H2      = h2;  %  H2 norm of the normial system
    info.problem = sol.problem;
    if sol.problem == 0
        [info,K] = ioppost(CYv,CUv,CWv,CZv,N,m,p,info);  
    else
         K = [];    
    end
    

end


% ========================================================================
%           Nested functions
% ========================================================================

function [cost,const] = iopcost2(const,CYv,CUv,CWv,CZv,Q,R,T1,T2,T3,K0,N,opts)

    p = size(K0,2);
    %------------------------------
    % Data in state space
    %------------------------------
    hatIp = zeros(N*p,p);  hatIp(1:p,1:p) = eye(p);
    Zp   = diag(ones(N-1,1), -1);     % downshift operator
    Zp   = kron(Zp,eye(p));
    hatY = CYv(:,p+1:end);
    Y0   = CYv(:,1:p);
    hatU = CUv(:,p+1:end);
    U0   = CUv(:,1:p);                   % this one must be zero

    
    %% cost in |Q Y T1|
    
    A1 = [Zp hatIp*T1.C;
          zeros(size(T1.A,1),N*p) T1.A];
    B1 = [hatIp*T1.D; T1.B];
    C1 = Q^(1/2)*[hatY Y0*T1.C];
    D1 = Y0*T1.D;
    
    [n1,m1] = size(B1); % n state dimension, m inpute dimension; p output dimension
    [p1,~] = size(C1);
    
    X1 = sdpvar(n1); gamma1 = sdpvar(1);
    Z1 = sdpvar(p1);
    
    cons1 = [X1 zeros(n1,m1) A1'*X1;
             zeros(m1,n1) gamma1*eye(m1) B1'*X1;
             X1*A1  X1*B1  X1];
    cons2 = [X1 zeros(n1,m1) C1';
             zeros(m1,n1) gamma1*eye(m1) D1';
             C1  D1  Z1];
    
    const = [const, cons1 >= 0, cons2>=0, trace(Z1) <= gamma1];  % constraints
    
    cost = gamma1^2;
    
%     %% cost in |R^{1/2} (K0Y + U)T1|
%     A2 = A1;
%     B2 = B1;
%     C2 = R^(1/2)*[K0*hatY+hatU (K0*Y0+U0)*T1.C];
%     D2 = [(K0*Y0+U0)*T1.D];
%     
%     [n1,m1] = size(B2); % n state dimension, m inpute dimension; p output dimension
%     [p1,~]  = size(C2);
% 
%     X2 = sdpvar(n1); gamma2 = sdpvar(1);
%     Z2 = sdpvar(p1);
%     
%     cons3 = [X2 zeros(n1,m1) A2'*X2;
%              zeros(m1,n1) gamma2*eye(m1) B2'*X2;
%              X2*A2  X2*B2  X2];
%     cons4 = [X2 zeros(n1,m1) C2';
%              zeros(m1,n1) gamma2*eye(m1) D2';
%              C2  D2  Z2];
%     
%     const = [const, cons3 >= 0, cons4>=0, trace(Z2) <= gamma2];  % constraints
%     cost = cost + gamma2^2;                                      % cost
%     



end
