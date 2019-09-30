function [Nr,Mr,Vl,Ul,Ml,Nl,Ur,Vr] = doubly_coprime_factorization_discrete(A,B,C)
% Computes a doubly coprime factorization of G using "?A Connection Between
% State-Space and Doubly Coprime Fractional Representations" by C. N. Nett,
% C. A. Jacobson and M. J. Barls, TAC 1984


n  = size(A,1);
m  = size(B,2);
p  = size(C,1);
polesK=rand([1,n])-rand([1,n]);
polesF=rand([1,n])-rand([1,n]);
K  = place(A,B,polesK);
F  = place(A',C',polesF)';

z=sym('z');
Nr  = C*inv(z*eye(n)-A+B*K)*B;
Mr  = eye(m)-K*inv(z*eye(n)-A+B*K)*B;
Vl  = -K*inv(z*eye(n)-A+F*C)*F;
Ul  = eye(m)+K*inv(z*eye(n)-A+F*C)*B;
Ml  = eye(m)-C*inv(z*eye(n)-A+F*C)*F;
Nl  = C*inv(z*eye(n)-A+F*C)*B;
Ur  = eye(m)+C*inv(z*eye(n)-A+B*K)*F;
Vr  = -K*inv(z*eye(n)-A+B*K)*F;

%Xl=V;
%Yl=-U;
%Ur=M;
%Yr=-U_tilda;
%Vl=N_tilda;
%Ul=M_tilda;
%Vr=N;
%Xr=V_tilda;

end

