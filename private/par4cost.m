function cost = par4cost(XYv,UYv,Q,R,N,C,n)
    %%H2cost
    cost = 0;
    [p,~] = size(Q);
    [m,~] = size(R);
    
    

    for t = 2: N+1
        cost = cost + trace(XYv(:,(t-1)*n+1:t*p)'*C'*Q*C*XYv(:,(t-1)*n+1:t*p))+trace(UYv(:,(t-1)*m+1:t*m)'*R*UYv(:,(t-1)*m+1:t*m));
    end
    cost = cost + trace((C*XYv(:,(1-1)*n+1:1*p)+eye(p))'*Q*(C*XYv(:,(1-1)*n+1:1*p)+eye(p)))+trace(UYv(:,(1-1)*m+1:1*m)'*R*UYv(:,(1-1)*m+1:1*m)); %%left out  because of the identity, added again separately
end
