
%% test files

n = 1;
m = 1;
p = 1;

loop = 1;
while true
    A  = randi([-15,15],n,n)/5; B  = randi([-5,5],n,m);
    C  = randi([-5,5],p,n);    % D  = randi([-2,2],p,m);

    K0 = rand(m,p);
    
    if max(abs(eig(A + B*K0*C))) < 1
        break;
    end
    loop = loop + 1
end

abs(eig(A))
abs(eig(A+B*K0*C))

opts.N = 5;
[K,h2,info] = iop(A,B,C,Q,R,K0,opts);

% G = ss(A,B,C,[],[]);  
% T1 = ss(A + B*K0*C,B*K0,C,eye(p),[]);
% T2 = ss(A + B*K0*C,B,C,zeros(p,m),[]);
% T3 = ss(A + B*K0*C,B,K0*C,eye(m),[]);
% 
% minreal(T1 - (eye(p) - G*K0)^(-1))
% minreal(T2 - (eye(p) - G*K0)^(-1)*G)
% minreal(T3 - (eye(m) - K0*G)^(-1))



