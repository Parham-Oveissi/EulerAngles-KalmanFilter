clc; clear; close all

dt = 0.1;
end_time = 20;
time = 0:dt:end_time;

IC = [0 0 0];

omega = [0; 0; 10];

Process_STD = 0.1;
Process_Var = Process_STD^2;

Measurement_noise_std = 0.5;
Measurement_noise_var = Measurement_noise_std^2;

[Time, Disturbed_Theta] = ode45(@(Time, Theta) Euler_Angular_Velocity_ODE(Time, Theta, omega, Process_STD), time, IC);

Noisy_Measurement = Disturbed_Theta + Measurement_noise_std*randn(size(Disturbed_Theta));


%% Kalman Filter Implementation
Estimate = [0; 0; 0]; % Initial state Estimate

P = 25*eye(length(Estimate)); % Initial state Covariance
Estimates_Vec(1,:) = Estimate;

Q = Process_Var*eye(size(Disturbed_Theta,2));
R = Measurement_noise_var*eye(size(Noisy_Measurement,2));


for kk = 2:length(Disturbed_Theta)
    Estimate = Estimate + dt*f(Estimate,omega);
    

    [A, B] = local_linearizer(Estimate, omega, dt);
    
    P = A*P*A' + Q;
    C = eye(size(Noisy_Measurement,2));
    K = (P*C')/(C*P*C'+ R);
    Estimate = Estimate + K*(transpose(Noisy_Measurement(kk,:)) - C*Estimate);
    P = P - K*C*P;

    Estimates_Vec(kk,:) = Estimate;

end

%% Plots
figure(1)
plot_EulerAngles(Time, Noisy_Measurement)
figure(2)
plot_EulerAngles(Time, Estimates_Vec)
pause(1)

figure(3)
Animate_attitude(Time, Noisy_Measurement, 'deg', 1e-10)
figure(4)
Animate_attitude(Time, Estimates_Vec, 'deg', 1e-10)

