
function cost = slscost(CNv,CUv,Q,R,C,T)
    %%H2cost in SLS
    %cost = 0;
    p = size(C,1);
    n = size(Q,1);
    m = size(R,1);
    
    cost = trace((CNv(:,1:p)'*C'+eye(n))*Q*(C*CNv(:,1:p)+eye(n))) + trace(CUv(:,1:m)'*R*CUv(:,1:m));
    Qc = C'*Q*C;
    for t = 2:T + 1
        cost = cost + trace(CNv(:,(t-1)*p+1:t*p)'*Qc*CNv(:,(t-1)*p+1:t*p)) + ...
                trace(CUv(:,(t-1)*m+1:t*m)'*R*CUv(:,(t-1)*m+1:t*m));
    end
end