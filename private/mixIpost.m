function [info,K] = mixIpost(YXv,YYv,UXv,UYv,N,m,p,info)
% iop - post processing

        info.var.YX  = (value(YXv));
        info.var.YY  = (value(YYv));
        info.var.UX  = (value(UXv));
        info.var.UY  = (value(UYv));
        K           = Kstatereal(info.var.UY,info.var.YY,N);  % state-space realization
        info.Ks     = K;
        
        %  closed-loop responses from optimization
%         z = tf('z');
%         Y = zeros(p,p); U = zeros(m,p);
%         W = zeros(p,m); Z = zeros(m,m);
%         for k = 1:N+1
%             Y = Y + info.var.Y(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
%             U = U + info.var.U(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
%             W = W + info.var.Y(:,m*(k-1)+1:m*k)*z^(1-k);  % FIR
%             Z = Z + info.var.U(:,m*(k-1)+1:m*k)*z^(1-k);  % FIR
%         end
%         info.cl.Y = Y;
%         info.cl.U = U;
%         info.cl.W = W;
%         info.cl.Z = Z;
end

