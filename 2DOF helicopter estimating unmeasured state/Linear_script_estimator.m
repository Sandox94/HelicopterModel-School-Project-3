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
theta_op = -10*pi/180;              % Pitch operation point
psi_op = pi/2;                    % Yaw operation point
w_theta_op = 0; w_psi_op = 0;           % Pitch and Yaw-acceleration operation point
Vmp_op = (Kyy*mheli*g*cos(theta_op)*lcm)/(Kyy*Kpp-Kyp*Kpy); % Volt-motor-pitch operation point
Vmy_op = (Kyp*Vmp_op)/Kyy;      % Volt-motor-yaw operation point

%%Simplify some operation
J1= Jeqp+mheli*(lcm.^2);
J2= Jeqy+mheli*(cos(theta_op).^2)*(lcm.^2);
K = (Kyp*Vmp_op-Kyy*Vmy_op)*2*mheli*(lcm.^2)*sin(theta_op)*cos(theta_op);

%%Matrix operation
Ac1= [0;0;(mheli*g*lcm*sin(theta_op))/J1;K/(J2.^2)];
Ac2= [0,1,0;0,0,1;0,-Bp/J1,0;0,0,-By/J2];
Ac= [Ac1,Ac2];
Bc= [0,0;0,0;Kpp/J1,-Kpy/J1;Kyp/J2,-Kyy/J2];
Cc= [1,0,0,0;0,1,0,0];
G = eye(4);
Q=diag([10,20,40,30]);
R= diag([0.1,0.1]);

%%Kalman filter
plant = ss(Ac,[Bc,G],Cc,0);         %DC = 0 since it is not present
[kalmf,L,P,M] = kalman(plant,Q,R);




