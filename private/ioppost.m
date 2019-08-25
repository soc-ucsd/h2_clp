function [info,K] = ioppost(CYv,CUv,CWv,CZv,N,m,p,info)
% iop - post processing

        info.var.Y  = (value(CYv));
        info.var.U  = (value(CUv));
        info.var.W  = (value(CWv));
        info.var.Z  = (value(CZv));
        K           = Kstatereal(info.var.U,info.var.Y,N);  % state-space realization
        info.Ks     = K;
        %  closed-loop responses from optimization
        z = tf('z');
        Y = zeros(p,p); U = zeros(m,p);
        for k = 1:N+1
            Y = Y + info.var.Y(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
            U = U + info.var.U(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
        end
        info.cl.Y = Y;
        info.cl.U = U;
end
