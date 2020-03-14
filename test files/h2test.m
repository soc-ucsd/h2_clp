
% test h2 norm computation

n = 4;
m = 2;
p = 1;

A  = randi([-15,15],n,n)/5; B  = randi([-2,2],n,m);
C  = randi([-2,2],p,n);     D  = randi([-2,2],p,m);

A  = A./max(abs(eig(A)))/1.1;  % stable plant

sys1 = ss(A,B,C,D,[]);


h2 = h2compute(A,B,C,D);

[norm(sys1,2), h2]
