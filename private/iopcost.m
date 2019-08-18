function cost = iopcost(CYv,CUv,Q,R,N)
    %%H2cost
    cost = 0;
    [p,~] = size(Q);
    [m,~] = size(R);

    for t = 1: N+1
        cost = cost + trace(CYv(:,(t-1)*p+1:t*p)'*Q*CYv(:,(t-1)*p+1:t*p))+trace(CUv(:,(t-1)*m+1:t*m)'*R*CUv(:,(t-1)*m+1:t*m));
    end
end