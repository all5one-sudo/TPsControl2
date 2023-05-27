using ControlSystemsBase, RobustAndOptimalControl, ForwardDiff, LinearAlgebra, Plots

function cartpole(x, u)
    mc, mp, l, g = 1.0, 0.2, 0.5, 9.81

    q  = x[1:2]
    qd = x[3:4]

    s = sin(q[2])
    c = cos(q[2])

    H = [mc+mp mp*l*c; mp*l*c mp*l^2]
    C = [0.1 -mp*qd[2]*l*s; 0 0]
    G = [0, mp * g * l * s]
    B = [1, 0]

    qdd = -H \ (C * qd + G - B * u[1])
    return [qd; qdd]
end

nu = 1    # number of control inputs
nx = 4    # number of states
ny = 2    # number of outputs (here we assume that the cart position and the pendulum angle are measurable)

x0 = [0, π, 0, 0]
u0 = [0]

Ac = ForwardDiff.jacobian(x->cartpole(x, u0), x0)
Bc = ForwardDiff.jacobian(u->cartpole(x0, u), u0)
Cc = [1 0 0 0; 0 1 0 0]
Λ = Diagonal([0.4, deg2rad(25)]) # Maximum output ranges
Cc = Λ\Cc # This normalizes expected outputs to be ∈ [-1, 1], a good practice for MIMO systems

sys = ss(Ac, Bc, Cc, 0)

pzmap(sys)

P = sminreal(sys[2,1]) # Position state goes away, not observable

pzmap(P)

C = pid(1,0,0)

w = exp10.(LinRange(-2.5, 3, 500))
function pid_marginplot(C)
    f1 = marginplot(P*C, w)
    vline!([2*4.85], sp=1, lab="Fundamental limitation", l=(:dash, :black))
    ylims!((1e-3, 1e2), sp=1)
    f2 = nyquistplot(P*C)
    plot(f1, f2)
end
pid_marginplot(C)

C = pid(20, 0, 0)
pid_marginplot(C)

C = pid(20, 0, 0.2, Tf=0.01)
pid_marginplot(C)

isstable(feedback(P*C))

C = pid(20, 1.25, 0.2, Tf=0.01)
pid_marginplot(C)

isstable(minreal(feedback(P*C)))

plot(step(feedback(P,C), 8), ylab="ϕ")

f1 = gangoffourplot(P,C,w, Ms_lines=[1.4], Mt_lines=[1.5])
f2 = nyquistplot(P*C, Ms_circles=[1.4], Mt_circles=[1.5], ylims=(-2, 2), xlims=(-4,1))
plot(f1, f2, size=(1000,800))

Pe = ExtendedStateSpace(sys, C2 = sys.C[2:2, :]) # Indicate that we can only measure the pendulum angle
Gecl = feedback(Pe, ss(C)) |> minreal
plot(step(Gecl, 8), ylab=["Cart pos" "ϕ"])

dampreport(sys)

desired_poles = [-4.85, -4.85, -5, -5]
L = place(sys, desired_poles, :c)

Bw = [0 0; 0 0; 1 0; 0 1]
R1 = Bw*I(2)*Bw'
R2 = 0.0001I(ny)
K = kalman(sys, R1, R2)

controller = observer_controller(sys, L, K)
@assert isstable(controller)
@assert isstable(feedback(sys * controller))

nyquistplot(controller*sys, w, Ms_circles=[2.7], Mt_circles=[3], xlims=(-2, 2), ylims=(-1, 3))

hinfnorm2(input_comp_sensitivity(sys, controller))

gangoffourplot(sys, controller, w, xlabel="", sigma=false, titlefont=8)

Kgmf, γ, info = glover_mcfarlane(sys, 1.05; W1=controller)
@assert isstable(Kgmf)

f1 = bodeplot([controller*sys, Kgmf*sys], w, plot_title="Input Loop transfers", lab=["Pole placement" "" "GMF" ""]); vline!([2*4.85], sp=1, lab="Fundamental limitation", l=(:dash, :black))
f2 = nyquistplot([controller*sys, Kgmf*sys], xlims=(-4, 4), ylims=(-1, 5), Ms_circles=[2.7], Mt_circles=[3], lab=["Pole placement" "GMF"])
f3 = bodeplot(controller, w, lab="Pole placement")
bodeplot!(Kgmf, w, plot_title="Controllers", lab="GMF", legend=:bottomleft)
f4 = sigmaplot([
    input_sensitivity(sys, controller),
    input_sensitivity(sys, Kgmf)
    ], w, title="Input S", lab=["Pole placement" "GMF"], legend=:bottomright)
plot(f1,f2,f3,f4, size=(1000,1000))

gangoffourplot(sys, [controller, Kgmf], w, xlabel="", sigma=false, titlefontsize=8)

sigmaplot(sensitivity.(Ref(sys), [controller, Kgmf]), w, lab=["Pole placement" "GMF"], legend=:bottomright)

dmf1 = plot(diskmargin(sys*controller), title="Simultaneous Output diskmargin", lab="Pole placement")
dmf2 = plot(diskmargin(controller*sys), title="Input diskmargin", lab="Pole placement")
plot!(dmf1, diskmargin(sys*Kgmf), title="Simultaneous Output diskmargin", lab="GMF")
plot!(dmf2, diskmargin(Kgmf*sys), title="Input diskmargin", lab="GMF")
plot(dmf1, dmf2)

rampsim(sys) = lsim(sys[:, 1], (x,t)->[min(t, 1)], 0:0.01:3.5, method=:zoh) # helper function to simulate ramp input
plot([
    rampsim(feedback(sys*controller)),
    rampsim(feedback(sys*Kgmf)),
], ylab=["Pos" "Angle"], plot_title="Position command step response", lab=["Pole placement" "" "GMF" ""], legend=:bottomright)

plot([
    impulse(feedback(sys, controller), 8, method=:zoh),
    impulse(feedback(sys, Kgmf), 8, method=:zoh),
], ylab=["Pos" "Angle"], plot_title="Disturbance step response", lab=["Pole placement" "" "GMF" ""], legend=:bottomright)

sim_pp  = rampsim(feedback(sys*controller))
sim_gmf = rampsim(feedback(sys*Kgmf))
@gif for i in 1:3:length(sim_pp.t)
    p, a = sim_pp.y[:, i]
    plot([p; p - 0.5sin(a)], [0; 0.5cos(a)], lw = 1, markershape = :square, markersize=[8, 1], lab = "PP", dpi = 200, size=(700, 200))

    p, a = sim_gmf.y[:, i]
    plot!([p; p - 0.5sin(a)], [0; 0.5cos(a)], lw = 1, markershape = :square, markersize=[8, 1], lab = "GMF", xlims = (-1, 2.2),
        ylims = (-0.2, 0.6), title = "Inverted pendulum control",
        dpi = 200, aspect_ratio = 1)
end

Cp = pid(-1, 5; state_space=true) # Position controller
Gtotal = feedback(Gecl, Cp, Y1=[1]) # Indicate that the outer controller can only see the cart position (y1)
plot(step.([Gecl, Gtotal], 20), lab=["Only angle feedback" "" "Cascade control" ""], legend=:bottomright, ylims=[(-5, 0.2) (-Inf, Inf)])