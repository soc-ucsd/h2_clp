clear all;

syms z

n = 3;
T = 3;
I = eye(3);

%% We need a nihilpotent A to satisfy the constraint M,N in z^-1RH_inf, L in RH_inf
B = eye(3);
C = eye(3);


%1) Create random FIR R, Delta1, Delta3, Delta4. Ensure that FIR/strictly
%proper requirements are satisfied by using some tricks. My tricks requires A
%nihilpotent, B=C=I
%2) Check if CL = R(I+\hat{Delta})^-1 is stable.

found = 0;

while(found == 0)
    
    A = [0 0 0;
    randi([-20,20]) 0 0;
    randi([-20,20]) randi([-20,20]) 0]

    R_tilde = sym(ones(n,n)); 
    Delta3  = sym(ones(n,n));
     
    for i = 1:n
        for j = 1:n
            for k = 1:T-1
                R_tilde(i,j) = R_tilde(i,j)*(z-randi([-5,5]));
                Delta3(i,j) = Delta3(i,j)*(z-randi([-1,1]));
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
    simplify((z*eye(3) - A)*R - B*M - eye(3)) % constraint 1
    simplify((z*eye(3) - A)*N - B*L)          % constraint 2
    simplify(R*(z*eye(3) - A) - N*C - eye(3)) % constraint 3
    simplify(M*(z*eye(3) - A) - L*C)          % constraint 4

    
    K = simplify(L - M*R^(-1)*N);               % controller
    CL = simplify((z*eye(3) - A - B*K*C)^(-1)); % closed-loop systems
  
    [~, d] = numden(det(CL));
    Cleig = roots(sym2poly(d));    % this is a subset of the closed-loop eigenvalues
    
    maxeig = max(abs(Cleig))
    
    if maxeig > 1
        fprintf('\n\nThe maximum eigenvalue of CL has norm %f .\n', maxeig )
        break;
    end
    
end



