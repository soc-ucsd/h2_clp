
%% Example for state-feedback 
clc;
n = 5;
% system
A = rand(n);

B = zeros(5,1);
B(end) = rand(1);

C = eye(n);

[rank(ctrb(A,B)),rank(obsv(A,C))]

% dynamical controller
Ak = A';
Bk = C';
Ck = B';
Dk = rand(1,n)

% closed-loop matrix Acl

Acl = [A+B*Dk*C  B*Ck;
       Bk*C     Ak];

eig(Acl)
hB = [eye(n);zeros(n)]
ctrb(Acl,hB)
rank(ctrb(Acl,hB))
