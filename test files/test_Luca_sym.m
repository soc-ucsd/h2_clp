clear all;


%% We need a nihilpotent A to satisfy the constraint M,N in z^-1RH_inf, L in RH_inf
syms z
n = 2;
T = 4;
I = eye(n);
B = eye(n);
C = eye(n);


%1) Create random FIR R, Delta1, Delta3, Delta4. Ensure that FIR/strictly
%proper requirements are satisfied by using some tricks. My tricks requires A
%nihilpotent, B=C=I
%2) Check if CL = R(I+\hat{Delta})^-1 is stable.

found = 0;

while found == 0
    
    A = zeros(n);
    R_tilde = sym(ones(n,n)); 
    Delta3  = sym(ones(n,n));
     
    for i = 1:n
        for j = 1:n
            for k = 1:T-1
                R_tilde(i,j) = R_tilde(i,j)*(z-randi([-6,6]));
                Delta3(i,j) = Delta3(i,j)*(z-randi([-2,2]));
            end
            
            %make them strictly proper.
            R_tilde(i,j) = R_tilde(i,j)/z^T;
            Delta3tf(i,j) = 0.001*Delta3(i,j)/z^T;  %very small
        end
    end
    
    Delta3 = Delta3tf.*I; %for simplicity, Delta3 is diagonal.
    
    R = simplify(inv(z*I-A)*(I+R_tilde));  %Trick: By defining R as such, we ensure that Mtf = R_tilde - Delta1, and so Mtf is guaranteed to lie in z^-1RH_inf and FIR. We need A nihilpotent to ensure Rtf is itself FIR.
    M = R_tilde;
    N = simplify(R * (z*I - A) - I - Delta3);
    L = simplify(M * (z*I - A));
   
    %% constraint 
    simplify((z*eye(n) - A)*R - B*M - eye(n)) % constraint 1
    simplify((z*eye(n) - A)*N - B*L)          % constraint 2
    simplify(R*(z*eye(n) - A) - N*C - eye(n)) % constraint 3
    simplify(M*(z*eye(n) - A) - L*C)          % constraint 4

    
    K = simplify(L - M*R^(-1)*N);               % controller
    CL = simplify((z*eye(n) - A - B*K*C)^(-1)); % closed-loop systems
  
    [~, d] = numden(det(CL));
    Cleig = roots(sym2poly(d));    % this is a subset of the closed-loop eigenvalues
    
    maxeig = max(abs(Cleig))
    
    if maxeig > 1
        fprintf('\n\nThe maximum eigenvalue of CL has norm %f .\n', maxeig)
        break;
    end
    
end

%% Check the formula in Theorem 6 of our paper 

Delta1 = zeros(n);
Delta4 = zeros(n);

hDelta = Delta3 + R*(eye(n) + Delta1)^(-1)*(B*Delta4 - (z*eye(n) - A)*Delta3);

CLn = (eye(n) + hDelta)^(-1)*R*(eye(n) + Delta1)^(-1);  % see eq.(41) in our paper

simplify(CLn - CL)    % this one should be zero

% --------------------------------------
%           Symbolic operations
% --------------------------------------

% step 1: check SLP consitriants  
simplify((z*eye(n) - A)*R - B*M - eye(n)) % constraint 1
simplify((z*eye(n) - A)*N - B*L)          % constraint 2
simplify(R*(z*eye(n) - A) - N*C - eye(n)) % constraint 3
simplify(M*(z*eye(n) - A) - L*C)          % constraint 4

% step 2: check controller & closed-loop    
K = simplify(L - M*R^(-1)*N);               % controller
CL = simplify((z*eye(n) - A - B*K*C)^(-1)); % closed-loop systems

[~, d] = numden(det(CL));
Cleig = roots(sym2poly(d));    % this is a subset of the closed-loop eigenvalues

maxeig = max(abs(Cleig));
fprintf('\n\nThe maximum eigenvalue of CL has norm %f .\n', maxeig)
