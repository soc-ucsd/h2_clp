function cost = par3cost(YYv,UYv,Q,R,N)
    %%H2cost
    cost = 0;
    [p,~] = size(Q);
    [m,~] = size(R);

    for t = 1: N+1
        cost = cost + trace(YYv(:,(t-1)*p+1:t*p)'*Q*YYv(:,(t-1)*p+1:t*p))+trace(UYv(:,(t-1)*m+1:t*m)'*R*UYv(:,(t-1)*m+1:t*m));
    end
end
