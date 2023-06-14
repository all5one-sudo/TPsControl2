%Caso de estudio 3. Sistema no lineal de cuatro variables de estado
close all;clear all; clc;
color = 'r';

% Se declaran los parámetros del sistema
m = 0.1;
Fricc = 0.1;
F = 0.1;
l = 1.6; 
g = 9.8;
M = 1.5;

% Modelos en Tiempo Continuo
% Matrices de lazo abierto en el equilibrio estable
Ac=[0       1               0       0;  %x1=delta - desplazamiento
    0    -Fricc/M         -m*g/M    0;  %x2=delta_p
    0       0               0       1;  %x3=phi - angulo
    0  -Fricc/(l*M)  -g*(m+M)/(l*M) 0];  %x4=phi_p
    
Ac = [0 1 0 0; 0 -F/M -m*g/M 0; 0 0 0 1; 0 -F/(l*M) -g*(m+M)/(l*M) 0]
Bc = [0; 1/M; 0; 1/(l*M)]
Cc = [1 0 0 0; 0 0 1 0];
Dc = [0];


%Condiciones iniciales
phi(1) = pi;               %Ángulo inicial
deltaRef = 10;               %Posición de referencia           
flag = 0;

%Tiempos
Ts = 1e-2;          
T = 15;            
At = 1e-4;          
Kmax = T/Ts;

%Pasaje de tiempo continuo a tiempo discreto
sys = ss(Ac, Bc, Cc, Dc);
dSys=c2d(sys,Ts,'zoh');

%TIEMPO DISCRETO
%Matrices de lazo abierto en el equilibrio estable
A = dSys.a;
B = dSys.b;

%Integrador se mide desplazamiento
Cref = Cc(1,:);
Aamp1=[A,zeros(4,1);-Cref*A,eye(1)];
Bamp1=[B;-Cref*B];

% %Matriz de Ma y Mc
% Ma = [BB AA*BB AA^2*BB AA^3*BB AA^4*BB AA^5*BB]; Alcanzabilidad = rank(Ma);
% Mc = [BB AA*BB AA^2*BB AA^3*BB AA^4*BB AA^5*BB AA^6]; Controlabilidad = rank(Mc); 


%parametros DLQR   de 0 a 10 m
dd = [.1 1e-2 1 .1 0.0000093]; %Desplazamiento, Velocidad, Angulo, Velocidad angular, Integrador
QQ = diag(dd);
%%RR = .0095;    
RR = 1.9e-4;

KK = dlqr(Aamp1,Bamp1,QQ,RR);
K = KK(1:4);
KI = -KK(5);

% Controlador  de 10 a 0 m
%matrices para m*10
m1 = m*10;

Ac2=[0       1               0       0;  %x1=delta - desplazamiento
    0    -Fricc/M         -m1*g/M    0;  %x2=delta_p
    0       0               0       1;  %x3=phi - angulo
    0  -Fricc/(l*M)  -g*(m1+M)/(l*M) 0];  %x4=phi_p

sys2 = ss(Ac2, Bc, Cc, Dc);
dSys2=c2d(sys2,Ts,'zoh'); 

%Matrices de lazo abierto en el equilibrio estable
A2 = dSys2.a;
B_m2 = dSys2.b;
C_m2 = dSys2.c;

Aamp2=[A2,zeros(4,1);-Cref*A2,eye(1)];% para el integrador de m2

% parametros DLQR
dd_m2 = [.1 1e-3 1e-3 .1 0.001]; %Desplazamiento, Velocidad, Angulo, Velocidad angular, Integrador
Q2 = diag(dd);
R2 = .008;                    

K2 = dlqr(Aamp2,Bamp1,Q2,R2);
Kp2 = K2(1:4);
Kint2 = -K2(5);


%Observador
Ao = A';
Bo = Cc';
Co = B';


%parametros DLRQ- Observador
do = [0.001 1000 0.5 0.0001]; %Desplazamiento, Velocidad, Ã?ngulo, Velocidad angular
Qo = diag(do); 
Ro = diag([80 10000]);

Kko = dlqr(Ao,Bo,Qo,Ro);
Ko=Kko';
t = 0;

x = [0;0;phi(1);0];
delta = x(1);
deltaP = x(2);
phi = x(3);
omega = x(4);

phiPP(1) = 0;
h = Ts/20;
i = 1;
deltaRef = 10; 
flag = 0;
v(1) = 0;
xHat = [0;0;pi;0];
xOp=[0 0 pi 0]';
reference(1) = 10;

for index=1:Kmax
        
    yOut=Cc*x;  %Salida de dos componentes
    yOutObs=Cc*(xHat-xOp);
    v(index+1)=v(index)+deltaRef-yOut(1);
    
    %Ley de control
    u1(index)=-K*(x-xOp)+KI*v(index+1); %color = '';%Sin observador
    %u1(index)=-K*xHat+KI*v(index+1);color = ''; %Con observador
    
    %Zona Muerta
    deadZone=1;
    if(abs(u1(index))<deadZone)
        u1(index)=0;               
    else
        u1(index)=sign(u1(index))*(abs(u1(index))-deadZone);
    end
    %-----------------------------------------------------
    
    for j=1:Ts/h
        
        u(i)=u1(index);
        p_pp=(1/(M+m))*(u(i)-m*l*phiPP*cos(phi(i))+m*l*omega(i)^2*sin(phi(i))-F*deltaP(i));
        phiPP=(1/l)*(g*sin(phi(i))-p_pp*cos(phi(i)));
        deltaP(i+1)=deltaP(i)+h*p_pp;
        delta(i+1)=delta(i)+h*deltaP(i);
        omega(i+1)=omega(i)+h*phiPP;
        phi(i+1)=phi(i)+h*omega(i);
        if(delta(i)>=9.99)
            if(flag==0)
                deltaRef=0;
                m=m*10;
                flag=1;
                K=Kp2;
                KI=Kint2;
            end
        end
        i=i+1;
        reference(i) = deltaRef;
    end
    x=[delta(i-1); deltaP(i-1); phi(i-1); omega(i-1)];
    xHat=A*xHat+B*u1(index)+Ko*(yOut-yOutObs)+xOp;
end

u(i)=u1(index);
t=0:h:T;

figure(1);
subplot(3,2,1); grid on; hold on;
plot(t,omega,color,'LineWidth',1.5);grid on; title('Velocidad angular \omega');

subplot(3,2,2); grid on; hold on;
plot(t,phi,color,'LineWidth',1.5); title('Ángulo \phi');xlabel('Tiempo');

subplot(3,2,3); grid on; hold on;
plot(t,delta,color,'LineWidth',1.5);hold on;plot(t,reference,'b--','LineWidth',0.5);
legend('Medicion','Referencia')
title('Posición grúa \theta');xlabel('Tiempo');

subplot(3,2,4); grid on; hold on;
plot(t,deltaP,color,'LineWidth',1.5);title('Velocidad de grúa \theta_p');

subplot(3,1,3); grid on; hold on;
plot(t,u,color,'LineWidth',1.5);title('Acción de control u');xlabel('Tiempo en Seg.');
 
 figure(2);
 subplot(2,1,1);grid on; hold on;
 plot(phi,omega,'');
 title('Ángulo vs Velocidad angular');
 xlabel('Ángulo');ylabel('Velocidad angular');
 
 subplot(2,1,2);grid on; hold on;
 plot(delta,deltaP, '');
 title('Distancia vs velocidad');
 xlabel('Distancia');ylabel('Velocidad');