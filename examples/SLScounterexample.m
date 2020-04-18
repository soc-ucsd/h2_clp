
%% counterexample

A = [0     0     0
     4     0     0
     3    -5     0];
 
B = eye(3);
C = eye(3);

z = tf('z');

%pvar z
 %z

syms z
 
inA = [1/z 0 0; 4/z^2 1/z 0; (3*z-20)/z^3 -5/z^2 1/z]; %% inverse of (zI - A)

% R = [(z^2 + z - 2)/z^3       (z + 2)/z^3   (z - 4)/z^3;
%     (5*z^2 + 7*z - 8)/z^4     (z^3 + z^2 + 2*z + 8)/z^4 (z^2 + z - 16)/z^4;
%   (z^3 - 23*z^2 - 41*z + 40)/z^5  (-4*z^3 + z^2 - 4*z - 40)/z^5 (z^4 + z^3 - 6*z^2 - 17*z + 80)/z^5];

M = [ (z - 2)/z^2, (z + 2)/z^2, (z-4)/z^2;
      (z + 3)/z^2 ,(z - 2)/z^2,(z - 3)/z^2;
     (z - 1)/z^2 ,(z + 3)/z^2,(z - 4)/z^2];

R = simplify(inv(z*eye(3)-A)*(eye(3)+M));
 

Delta3 = [(0.0001*z + 0.0001)/z^2 0 0;
 0  (0.0001*z + 0.0002)/z^2 0 ;
 0 0   (0.0001*z - 0.0003)/z^2];

N = simplify(R * (z*eye(3) - A)) - eye(3)- Delta3; 
L = simplify(M * (z*eye(3) - A));

simplify((z*eye(3) - A)*R - B*M - eye(3)) % constraint 1
simplify((z*eye(3) - A)*N - B*L) % constraint 2
simplify(R*(z*eye(3) - A) - N*C - eye(3)) % constraint 3
simplify(M*(z*eye(3) - A) - L*C) % constraint 4

K = simplify(L - M*R^(-1)*N)     % controller

% closed-loop system 
CL = simplify((z*eye(3) - A - B*K*C)^(-1))


