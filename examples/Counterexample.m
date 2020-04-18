
%% Failure example
clc;clear
syms z
syms d

%z = tf('z');

R = eye(2)/z + randi([-2 2],2,2)/z^2;% + randi([-2 2],2,2)/z^3;

delta = randi([-2 2],2,2)*d;

A = randi([-2 2],2,2);

simplify(R^(-1)*delta*R)

simplify((eye(2) - (z*eye(2) -A)*delta*R + R^(-1)*delta*R)^(-1))

tmp = simplify(R*(eye(2) - (z*eye(2) -A)*delta*R + R^(-1)*delta*R)^(-1))




% tmp0 = minreal(R^(-1)*delta*R);
% pole(tmp0)
% tzero(tmp0)
% 
% 
% A = eye(2);
% 
% 
% tmp = (eye(2) - (z*eye(2)-A)*delta*R + R^(-1)*delta*R);
% tmp = minreal(tmp^(-1));
% pole(tmp)
% 
% hR = minreal(R*tmp);
% pole(hR)
