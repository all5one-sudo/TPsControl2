using LinearAlgebra, Plots, ControlSystems, Polynomials
Laa = 5e-3;
J = 0.004;
Ra = 0.2;
Bm = 0.005;
Ki = 6.5e-5;
Km = 0.055;
A = [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0]
B = [1/Laa 0; 0 -1/J; 0 0]
C = [0 0 1]
D = [0 0]
M = [B[:, 1] A * B[:, 1] A^2 * B[:, 1]]
charPoly = fromroots(A)
W = [charPoly[1] charPoly[2] 1; charPoly[1] 1 0; 1 0 0]
T = M*W
Mca = inv(T)*A*T