F = 0.1
l = 1
g = 9.8
M = 1.5
m = 0.1
A = [0 1 0 0; 0 -F/M -m*g/M 0; 0 0 0 1; 0 F/(l*M) g*(M+m)/(l*M) 0];
B = [0; 1/M; 0; -1/(l*M)];
C = [1 0 0 0; 0 0 1 0];
D = [0; 0];
P = syslin("c",A,B,C,D)
Q = C'*C
R = 1;
K = lqr(P,Q,R)
Ac = [(A-B*K)]
Bc = [B]
Cc = [C]
Dc = [D]
S = syslin("c",Ac,Bc,Cc,Dc)
h = 1e-4;
simTime = 50;
t = 0:h:simTime;
u = 10*ones(t);
y=csim(u,t,S,[0; 0; 0; 0]);
clf;
plot(t',y');
xlabel(_("Tiempo [seg]"))
L = legend(["$\delta(t)$","$\phi(t)$"]);
L.font_size = 4

