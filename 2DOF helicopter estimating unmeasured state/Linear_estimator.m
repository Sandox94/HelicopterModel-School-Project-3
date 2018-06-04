%%Parameters of the system part 1
lcm = 0.015;    %[m]            % Distance between the pivot point and center of mass
mheli = 0.479;  %[kg]           % Total moving mass of heli
Jeqp = 0.0172; Jeqy= 0.0210;    %[kg/m2] % Moment of inertia pitch and yaw axis
g = 9.81;        %[m/s2]         % Acceleration due to gravity

%%Parameter of the system part 2
Kpp = 0.0556;    %[Nm/V]         % Torque constant on pitch for pitch motor
Kyy = 0.21084;   %[Nm/V]         % Torque constant on yaw for yaw motor
Kpy = 0.005;     %[Nm/V]         % Torque constant on pitch from yaw motor
Kyp = 0.15;      %[Nm/V]         % Torque constan on pitch from pithc motor
Bp = 0.01;       %[N/V]          % Damping friction factor about pitch
By = 0.08;       %[N/V]          % Damping friction fator about yaw

%%Operation point
theta_op = (-10*pi)/180;
psi_op = pi/2;
w_theta_op = 0;
w_psi_op = 0;
Vmp_op = (Kyy*mheli*g*cos(theta_op)*lcm)/(Kyy*Kpp-Kyp*Kpy);
Vmy_op = (Kyp*Vmp_op)/Kyy;
Uop = [Vmp_op;Vmy_op];

%%Simplify some complex operation
Co_1 = Jeqp+mheli*lcm^2;
Co_2 = Jeqy+mheli*cos(theta_op)^2*lcm^2;
Co_3 = (Kyp*Vmp_op-Kyy*Vmy_op)*2*mheli*lcm^2*sin(theta_op)*cos(theta_op);

%%Matrix operation
Ac = [  0,0,1,0;
        0,0,0,1;
        (mheli*g*lcm*sin(theta_op))/Co_1,0,-Bp/Co_1,0;
        Co_3/Co_2^2,0,0,-By/Co_2];

Bc = [  0,0;
        0,0;
        Kpp/Co_1, -Kpy/Co_1;
        Kyp/Co_2, -Kyy/Co_2];
Cc = [  1,0,0,0;
        0,1,0,0];
G = eye(4);
V = diag([0.1,0.1]);
W = diag([60,70,80,70]);

R = diag([7,10]);
Q = diag([5,1,5,1]);

%%Kalman filter
plant = ss(Ac,[Bc,G],Cc,0);
[kalmf,L,P,M] = kalman(plant,W,V);

%%LQR Grain matrix
[K] = lqr(Ac,Bc,Q,R);
