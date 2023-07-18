function plot_EulerAngles(Time, Theta)
plot(Time, Theta(:,1),'r', 'LineWidth',3)
hold on
plot(Time, Theta(:,2),'g', 'LineWidth',3)
plot(Time, Theta(:,3),'b', 'LineWidth',3)

xlabel('Time')
ylabel('Angle (deg)')
axis tight
grid on
legend('Roll (\phi)', 'Pitch (\Theta)', 'Yaw (\Psi)')
title('Attitude')