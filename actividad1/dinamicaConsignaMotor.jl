using ControlSystems
using Plots

# Constantes del motor dadas en la consigna
Laa = 366e-6;
J = 5e-9;
Ra = 55.6;
Bm = 0;
Ki = 6.49e-3;
Km = 6.53e-3;

# Matriz de Funciones de Transferencia MIMO
G = [tf([J, Bm],[J*Laa, (Bm*Laa+J*Ra), (Bm*Ra+Ki*Km)]) tf([Km],[J*Laa, (Bm*Laa+J*Ra), (Bm*Ra+Ki*Km)]);
     tf([Ki],[J*Laa, (Bm*Laa+J*Ra), (Bm*Ra+Ki*Km)]) tf([-Laa, -Ra],[J*Laa, (Bm*Laa+J*Ra), (Bm*Ra+Ki*Km)]);
     tf([Ki],[J*Laa, (Bm*Laa+J*Ra), (Bm*Ra+Ki*Km), 0]) tf([-Laa, -Ra],[J*Laa, (Bm*Laa+J*Ra), (Bm*Ra+Ki*Km), 0])];

t = 0:0.00001:0.6  # Vector de tiempo de 0 a 1 segundos con un paso de 0.01 segundos
ran = 60001 # Longitud del Vector
# Se crea la entrada de tensión
inputVa = zeros(Float64, ran)
for i in 1:ran
    if t[i] < 0.025
        inputVa[i] = 0.0
    elseif t[i] >= 0.025 && t[i] < 0.15
        inputVa[i] = 12.0
    else
        inputVa[i] = -12.0
    end
end
# Ahora se crea la entrada de torque
inputTL = zeros(Float64, ran)
for i in 1:ran
    if t[i] < 0.15
        inputTL[i] = 0;
    else
        inputTL[i] = -1.04e-3;
    end
end

# Se crea el vector entrada
u = [inputVa inputTL]'
# Se realiza la simulación lineal
y, t, _ = lsim(G, u, t)

# Se crea la figura
p = plot(layout = (4, 1))  

# Ahora se tienen los subplots
plot!(p[1], t, y[2,:], label="Velocidad angular")
plot!(p[2], t, y[1,:], label="Corriente de Armadura")
plot!(p[3], t, inputVa, label="Entrada de tensión")
plot!(p[4], t, inputTL, label="Entrada de torque")