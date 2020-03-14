function h2 = h2compute(A,B,C,D)
% Compute h2 norm of a discrete-time system 
%  x_{k+1} = A x_k + B u_k
%  y_k     = C x_k + D u_k
%  see the formula in 
%  Carsten W Scherer. An efficient solution to 
%     multi-objective control problems with lmi objectives.
%     Systems & control letters, 40(1):43¨C57, 2000.


    [n,m] = size(B); % n state dimension, m inpute dimension; p output dimension
    [p,~] = size(C);
    
    X     = sdpvar(n);
    gamma = sdpvar(1);
    Z     = sdpvar(p);
    
    cons1 = [X zeros(n,m) A'*X;
             zeros(m,n) gamma*eye(m) B'*X;
             X*A  X*B  X];
    cons2 = [X zeros(n,m) C';
             zeros(m,n) gamma*eye(m) D';
             C  D  Z];
    
    const = [cons1 >= 0, cons2>=0, trace(Z) <= gamma];  % constraints
    cost = gamma^2;                                      % cost
    
    yalmipOpts = sdpsettings('solver','mosek');
    sol        = optimize(const,cost,yalmipOpts);
    
    h2 = value(cost);

end

