clc; clear; close all

syms Phi Theta Psi omega1 omega2 omega3

dt = 0.1;
end_time = 20;
time = 0:dt:end_time;

IC = [2 1 3];

omega = [0; 0; 10];

Process_STD = 0.1;
Process_Var = Process_STD^2;

Measurement_noise_std = 0.01;
Measurement_noise_var = Measurement_noise_std^2;

[Time, Disturbed_Theta] = ode45(@(Time, Theta) Euler_Angular_Velocity_ODE(Time, Theta, omega, Process_STD), time, IC);
[~, Not_Disturbed_Theta] = ode45(@(Time, Theta) Euler_Angular_Velocity_ODE(Time, Theta, omega, 0), time, IC);


Noisy_Measurement = Disturbed_Theta + Measurement_noise_std*randn(size(Disturbed_Theta));


%% Kalman Filter Implementation
Estimate = [0; 0; 0]; % Initial state Estimate

P = 400*eye(length(Estimate)); % Initial state Covariance
Estimates_Vec(1,:) = Estimate;
% P_vec(1,:) = diag(P);

Q = Process_Var*eye(size(Disturbed_Theta,2));
R = Measurement_noise_var*eye(size(Noisy_Measurement,2));
C = eye(size(Noisy_Measurement,2));


[A_sym, B_sym] = local_linearizer(dt);

for kk = 2:length(Disturbed_Theta)
    A = double(subs(A_sym,[Phi Theta Psi omega1 omega2 omega3], ...
    [Estimate; omega]'));

    B = double(subs(B_sym,[Phi Theta Psi omega1 omega2 omega3], ...
    [Estimate; omega]'));
    Estimate = Estimate + dt*f(Estimate,omega);
    % [A, B] = local_linearizer(Estimate, omega, dt);
    
    
    P = A*P*A' + Q;

    K = (P*C')/(C*P*C'+ R);
    Estimate = Estimate + K*(transpose(Noisy_Measurement(kk,:)) - C*Estimate);
    P = P - K*C*P;

    Estimates_Vec(kk,:) = Estimate;
    % P_vec(kk,:) = diag(P);

end

Estimation_error = abs(Estimates_Vec - Disturbed_Theta);
Measurement_error = abs(Noisy_Measurement - Disturbed_Theta);
Model_error = abs(Not_Disturbed_Theta - Disturbed_Theta);
%% Plots
figure(1)
plot_EulerAngles(Time, Noisy_Measurement,'Attitude Measurement')
figure(2)
plot_EulerAngles(Time, Estimates_Vec,'Attitude Estimation')
pause(1)

figure(3)
Animate_attitude(Time, Noisy_Measurement, 'deg', 1e-10)
figure(4)
Animate_attitude(Time, Estimates_Vec, 'deg', 1e-10)

figure(5)
semilogy(Time, Estimation_error,'LineWidth',3)
xlabel('Time')
ylabel('Angle (deg)')
axis tight
grid on
legend('Roll (\phi) error', 'Pitch (\Theta) error', 'Yaw (\Psi) error')
title('Estimation Error')

% figure(6)
% semilogy(Time,P_vec, 'LineWidth',3)
% xlabel('Time')
% ylabel('Variance (deg^{2})')
% axis tight
% grid on
% legend('Roll (\phi) var', 'Pitch (\Theta) var', 'Yaw (\Psi) var')
% title('Estimate Variance')

figure(6)
semilogy(Time, Measurement_error,'LineWidth',3)
xlabel('Time')
ylabel('Angle (deg)')
axis tight
grid on
legend('Roll (\phi) error', 'Pitch (\Theta) error', 'Yaw (\Psi) error')
title('Measurement Error')

figure(7)
semilogy(Time, Model_error,'LineWidth',3)
xlabel('Time')
ylabel('Angle (deg)')
axis tight
grid on
legend('Roll (\phi) error', 'Pitch (\Theta) error', 'Yaw (\Psi) error')
title('Model Error')


figure(8);subplot(3,1,1);
plot(Time, Disturbed_Theta(:,1),'r', 'LineWidth',3)
hold on
plot(Time, Estimates_Vec(:,1),'r--', 'LineWidth',3)
ylabel('Roll (\phi)')
legend('Disturbed System', 'Estimates')

subplot(3,1,2);
plot(Time, Disturbed_Theta(:,2),'g', 'LineWidth',3)
hold on
plot(Time, Estimates_Vec(:,2),'g--', 'LineWidth',3)
ylabel('Pitch (\Theta)')
legend('Disturbed System', 'Estimates')

subplot(3,1,3);
plot(Time, Disturbed_Theta(:,3),'b', 'LineWidth',3)
hold on
plot(Time, Estimates_Vec(:,3),'b--', 'LineWidth',3)
ylabel('Yaw (\Psi)')
xlabel('Time')
legend('Disturbed System', 'Estimates')
