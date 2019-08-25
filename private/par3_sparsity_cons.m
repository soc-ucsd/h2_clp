%encode sparsity constraints

%{
%sparsity of UY
for(i=1:m)
    for(j=1:p)
        for(t=1:N+1)
            UYv(:,[(t-1)*m+1:t*p])=UYv(:,[(t-1)*m+1:t*p]).*T_struct;
        end
    end
end

%sparsity of YY

for(i=1:p)
    for(j=1:p)
        for(t=1:N+1)
            YYv(:,[(t-1)*p+1:t*p])=YYv(:,[(t-1)*p+1:t*p]).*R_struct;
        end
    end
end
%}



%sparsity of UY
UYs = sym('UY',[m p*(N+1)]);

UY_tf = zeros(m,p);
for t = 1:N+1
    UY_tf = UY_tf + UYs(:,[(t-1)*m+1:t*p])/z^(t-1);
end

for i=1:m      %ach1
    for j=1:p
        if(T_struct(i,j)==0)
            fprintf(' sparsity of UY:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/p );
            [num,~] = numden(UY_tf(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,vec(UYs));
            A_eqs          = double(A_eq);
            b_eqs          = double(b_eq);
            const    = [const, A_eqs*vec(UYv)== b_eqs];
        end
    end
end



%SPARSITY OF YY
YYs = sym('YY',[p p*(N+1)]);

YY_tf = zeros(p,p);
for t = 1:N+1
    YY_tf = YY_tf + YYs(:,[(t-1)*p+1:t*p])/z^(t-1);
end



for i=1:p      
    for j=1:p
        if(R_struct(i,j)==0)
            fprintf(' sparsity of Y:  Percentage %6.4f \n', 100*(p*(i-1)+j)/p/p );
            [num,~] = numden(YY_tf(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,vec(YYs));
            A_eqs          = double(A_eq);
            b_eqs          = double(b_eq);
            const    = [const, A_eqs*vec(YYv)== b_eqs];
        end
    end
end
