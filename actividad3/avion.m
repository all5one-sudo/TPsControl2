clc, clear all, close all;

% Caso de estudio 2. Sistema lineal de cuatro variables de estado.

% Parametros del sistema
a = 0.07;
b = 5;
c = 150;
w = 9;

%x1=alpha x2=phi x3=phi_p x4=h
A = [-a a 0 0; 0 0 1 0; w^2 -w^2 0 0; c 0 0 0];
B = [0; 0; b*w^2; 0];
C = [0 0 0 1;0 1 0 0];

% Polos deseados para el controlador
p1 = -15+15i;
p2 = -15-15i;
p3 = -0.5+0.5i;
p4 = -0.5-0.5i;
poles = [p1 p2 p3 p4];

 K=acker(A,B,[p1 p2 p3 p4]);

Q=diag([1 1000000 1 1]);R=1000000/1;
%K=lqr(A,B,Q,R);

%K = place(A,B,poles);

G=-inv(C(1,:)*inv(A-B*K)*B);
%Observador--------------------------------------
C=[0 0 0 1;0 1 0 0];
Ao=A';
Bo=C';
Co=B';

Qo=100*diag([1 100 1 100]);    Ro=100000000;

Ko=lqr(Ao,Bo,Qo,Ro);

% Ko=place(Ao,Bo,[-100 -200 -300 -400]);

%------------------------------------------------

%ganancia de prealimentacion:
G=-inv(C(1,:)*inv(A-B*K)*B);

%Variables
T=100; 
h=1e-4; 
t=0:h:(T-h);
pasos=T/h;

%referencias
%ref=-100;
 ref=100;

%condiciones iniciales
alpha(1)=0;
phi(1)=0;
phi_p(1)=0;
high(1)=500;
% high(1)=-500;
u(1)=0;
uu(1)=0;
estados=[alpha(1);
    phi(1);
    phi_p(1);
    high(1)];

    
Xop=[0 0 0 0]';

x=[alpha(1);phi(1);phi_p(1); high(1)];


alpha_obs(1)=0;
phi_obs(1)=0;
phip_obs(1)=0;
high_obs(1)=400;
xobs=[alpha_obs(1);phi_obs(1);phip_obs(1); high_obs(1)];
estados_obs=[alpha_obs(1);
    phi_obs(1);
    phip_obs(1);
    high_obs(1)];

for i=1:pasos
    
    u(i)= -K*estados_obs+G*ref;
    
    Alin=0.1;
    if abs(u(i))<Alin
        uu(i)=0;
    else
        % si le resto la alinealidad obtengo la zona muerta en el actuador
        % y linealidad. Sino le resto es un controlador todo o nada
        uu(i)=sign(u(i))*(abs(u(i)));
    end
    
    %Variables del sistema lineal
    alpha(i)= x(1);
    phi(i)= x(2);
    phi_p(i)= x(3);
    high(i)=x(4);
    
    %Sistema lineal
    xp=A*x+B*uu(i);
    x=x+h*xp;
    
    %------con Observador----------------------
  
    alpha_obs(i)= xobs(1);
    phi_obs(i)= xobs(2);
    phip_obs(i)= xobs(3);
    high_obs(i)= xobs(4);
    
    ysal_o = C * estados_obs;
    ysal   = C * estados;
    
    e=ysal-ysal_o;
    
    x_antp     = A*xobs+B*uu(i)+Ko'*e;%complicacion
    xobs       = xobs + x_antp*h;
    
    %-------------------------------------------
    estados=[alpha(i);
    phi(i);
    phi_p(i);
    high(i)];

    estados_obs=[alpha_obs(i);
    phi_obs(i);
    phip_obs(i);
    high_obs(i)];

end


%_______________SIN OBS---------------------------------------------------
%condiciones iniciales

u_so(1)=0;
estados_so=[alpha(1);
    phi(1);
    phi_p(1);
    high(1)];
    

x_so=[alpha(1);phi(1);phi_p(1); high(1)];

for i=1:pasos
    
    u_so(i) = -K*estados_so+G*ref;
    %Variables del sistema lineal
    alpha_so(i)= x_so(1);
    phi_so(i)= x_so(2);
    phi_p_so(i)= x_so(3);
    high_so(i)=x_so(4);
    
    %Sistema lineal
    xp_so=A*x_so+B*u_so(i);
    x_so=x_so+h*xp_so;
    estados_so=[alpha_so(i);
    phi_so(i);
    phi_p_so(i);
    high_so(i)];
end

    
%----------------------------------------------------------------------

subplot(3, 2, 1);
hold on
plot(t,alpha);
plot(t,alpha_so);
hold off
title('Angulo con la horizontal \alpha');
legend({'Con observador','Sin observador'})
xlabel('Tiempo (seg.)');
ylabel('rad');
grid on;

subplot(3, 2, 2);
hold on
plot(t,phi);
plot(t,phi_so);
hold off
title('Angulo de cabeceo \phi');
legend({'Con observador','Sin observador'})
xlabel('Tiempo (seg.)');
ylabel('Velocidad (m/s)');
grid on;

subplot(3, 2, 3);
hold on
plot(t,phi_p);
plot(t,phi_p_so);
hold off
title('velocidad de angulo de cabeceo \phi_p');
legend({'Con observador','Sin observador'})
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad/s)');
grid on;

subplot(3, 2, 4);
hold on
plot(t,high);
plot(t,high_so);
hold off
title('Altura h');
legend({'Con observador','Sin observador'})
xlabel('Tiempo (seg.)');
ylabel('metros');
grid on;

subplot(3, 1, 3);
hold on
plot(t,uu);
plot(t,u_so);
hold off
title('Accion de control u_t');
legend({'Con observador','Sin observador'})
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;