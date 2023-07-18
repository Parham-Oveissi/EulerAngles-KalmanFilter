function d_theta = Euler_Angular_Velocity_ODE(t, theta, omega, Process_STD)

Phi = theta(1);
Theta = theta(2);
Psi = theta(3);

S = [1 sind(Phi)*tand(Theta) cosd(Phi)*tand(Theta); ...
    0 cosd(Phi) -sind(Phi); ... 
    0 sind(Phi)*secd(Theta) cosd(Phi)*secd(Theta)];


% omega = [2*sind(0.01*t) ; -3*cosd(0.02*t); 0+sind(0.03*t)]*5;
% omega = [2*cosd(0.01*t) ; -2*cosd(0.02*t); 2*cosd(0.03*t)]*10;

% omega = [3; 3; 0];

d_theta = S*omega + Process_STD*randn(3,1);


