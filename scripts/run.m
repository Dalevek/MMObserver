%ctls.configureLogger()

s = Simulation("Simple example",-0.25,-0.25,-0.25);

% # run simulation without second layer (single model)
ra = runObserverMulti(s, "none","integral-finite");
% # run simulation with second layer (multi model)
rMMa = runObserverMulti(s, "RLS", "integral-finite");

time = ra('time');

x1  = ra('x1')
x2  = ra('x2')
x3  = ra('x3') 

xe1S = ra('xe1')
xe2S = ra('xe2')
xe3S = ra('xe3')

xe1M = rMMa('xe1')
xe2M = rMMa('xe2')
xe3M = rMMa('xe3')

e1 = ra('e1')
e2 = ra('e2')
e3 = ra('e3')

e1MM = rMMa('e1')
e2MM = rMMa('e2')
e3MM = rMMa('e3')

con = 10000.0

figure
plot(time,x1,time,xe1S,time,xe1M,'LineWidth',3)
xlabel('time $$[s]$$','Interpreter','Latex')
ylabel('motor speed $$[rad/s]$$','Interpreter','Latex')
legend({'motor speed $$\\omega(t)$$','estimated motor speed (SO) $$\hat{\omega}(t)$$','estimated motor speed (MO) $$\hat{\omega}(t)$$'},'Interpreter','Latex')

figure
plot(time,x2,time,xe2S,time,xe2M,'LineWidth',3)
xlabel('time $$[s]$$','Interpreter','Latex')
ylabel('motor current $[A]$','Interpreter','Latex')
legend({'motor current $i(t)$','estimated motor current (SO) $\hat{i}(t)$','estimated motor current (MO) $\\hat{i}(t)$'},'Interpreter','Latex')

figure
plot(time,x3,time,xe3S,time,xe3M,'LineWidth',3)
xlabel('time $$[s]$$','Interpreter','Latex')
ylabel('load torque $[Nm]$','Interpreter','Latex')
legend({'load torque $T_L(t)$','estimated load torque (SO) $\hat{T}_L(t)$','estimated load torque (MO) $\hat{T}_L(t)$'},'Interpreter','Latex')


figure
plot(time,e1,time,e1MM,'LineWidth',3)
xlabel('time $$[s]$$','Interpreter','Latex')
ylabel('observation error $e_1$','Interpreter','Latex')
legend({'single layer observer $e_1(t)$','multi layer observer $e_1(t)$'},'Interpreter','Latex')

figure
plot(time,e2,time,e2MM,'LineWidth',3)
xlabel('time $$[s]$$','Interpreter','Latex')
ylabel('observation error $e_2$','Interpreter','Latex')
legend({'single layer observer $e_2(t)$','multi layer observer $e_2(t)$'},'Interpreter','Latex')

figure
plot(time,e3,time,e3MM,'LineWidth',3)
xlabel('time $$[s]$$','Interpreter','Latex')
ylabel('observation error $e_3$','Interpreter','Latex')
legend({'single layer observer $e_3(t)$','multi layer observer $e_3(t)$'},'Interpreter','Latex')

