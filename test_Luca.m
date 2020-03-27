clear all;
clc;

z = tf('z');

n = 3;
T = 3;
I = eye(3);

%% We need a nihilpotent A to satisfy the constraint M,N in z^-1RH_inf, L in RH_inf
A = [0 0 0;
    randi([-20,20]) 0 0;
    randi([-20,20]) randi([-20,20]) 0]
B= eye(3)
C=eye(3)


%1) Create random FIR R, Delta1, Delta3, Delta4. Ensure that FIR/strictly
%proper requirements are satisfied by using some tricks. My tricks requires A
%nihilpotent, B=C=I
%2) Check if CL = R(I+\hat{Delta})^-1 is stable.
found = 0;


while(found == 0)
    R_tilde = tf(ones(n,n));
    Delta1 = tf(ones(n,n));
    Delta3 = tf(ones(n,n));
    Delta4 = tf(ones(n,n));
    
    for(i = 1:n)
        for(j = 1:n)
            for(k = 1:T-1)
                R_tilde(i,j) = R_tilde(i,j)*(z-randi([-5,5]));
                Delta1(i,j) = Delta1(i,j)*(z-randi([-5,5]));
                Delta3(i,j) = Delta3(i,j)*(z-randi([-5,5]));
                Delta4(i,j) = Delta4(i,j)*(z-randi([-5,5]));
            end
            
            %make them strictly proper.
            R_tilde(i,j) = R_tilde(i,j)/z^T;
            Delta1(i,j) = 0*Delta1(i,j)/z^(T);
            Delta3tf(i,j) = 0.0001*Delta3(i,j)/z^T;  %very small
            Delta4(i,j) = 0*Delta4(i,j)/z^T;
        end
    end
    
    Delta3tf = Delta3tf.*I; %for simplicity, Delta3 is diagonal.
    Delta3 = ss(Delta3tf);
    
    Rtf = inv(z*I-A)*(I+R_tilde);  %Trick: By defining R as such, we ensure that Mtf = R_tilde - Delta1, and so Mtf is guaranteed to lie in z^-1RH_inf and FIR. We need A nihilpotent to ensure Rtf is itself FIR.
    Mtf = R_tilde - Delta1;
    Ntf = Rtf * (z*I - A) - I - Delta3;
    Ltf = Mtf * (z*I - A) - Delta4;
    
    R = ss(Rtf);
    M = ss(Mtf);
    N = ss(Ntf);
    L = ss(Ltf);
    
    hat_Delta =ss(Delta3 - R * (z*I-A)*Delta3);
    inv_hat_Delta = inv(I+hat_Delta);
    %CL =
    %minreal(ss(inv(I+Delta3+R*inv(I+Delta1)*(Delta4-(z*I-A)*Delta3))*R*inv(I+Delta1)));
    %COMPLETE EXPRESSION
    
    CL = minreal(inv_hat_Delta*R); %only Delta3R
    
    eigs_CL = eig(CL.A);
    for(i=1:size(eigs_CL,1))
        eigs_CL(i) = norm(eigs_CL(i));
    end
    max_eig_CL = max(eigs_CL);
    
    if(max_eig_CL > 1)
        found = 1;
    end
    fprintf('\n\nThe maximum eigenvalue of CL has norm %f .\n', max_eig_CL )
end

fprintf('\nFOUND UNSTABLE\n')

[SV,W]=sigma(Rtf);
fprintf('\nNote that the norm of R is %f\n', max(max(SV)))



[SV,W]=sigma(Delta3tf-Rtf * (z*I-A)*Delta3tf);



fprintf('The norm of hat_delta is %f .\n', max(max(SV)))


eigs_inv_hat_Delta = eig(inv_hat_Delta.A);
for(i=1:size(eigs_inv_hat_Delta,1))
    eigs_inv_hat_Delta(i) = norm(eigs_inv_hat_Delta(i));
end
max_eig_eigs_inv_hat_Delta = max(eigs_inv_hat_Delta);

fprintf('The max eigenvalue of (I+hat_delta)^-1 is %f .\n', max_eig_eigs_inv_hat_Delta)

if(max(max(SV)) > 10000000)
    fprintf('NORM OF hat_Delta is infinite. Numerical issues. *RUN AGAIN*\n')
end



