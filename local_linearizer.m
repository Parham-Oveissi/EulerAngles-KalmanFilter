function [A_sym, B_sym] = local_linearizer(dt)
% function [A_sym, B_sym] = local_linearizer(theta, omega, dt)

syms Phi Theta Psi omega1 omega2 omega3


S = [1 sind(Phi)*tand(Theta) cosd(Phi)*tand(Theta); ...
    0 cosd(Phi) -sind(Phi); ... 
    0 sind(Phi)*secd(Theta) cosd(Phi)*secd(Theta)];

s_w = S*[omega1; omega2; omega3];

F = [Phi; Theta; Psi] + dt*s_w;

v1 = [Phi; Theta; Psi];
v2 = [omega1; omega2; omega3];

A_sym = jacobian(F,v1);
% A = double(subs(A_sym,[Phi Theta Psi omega1 omega2 omega3], ...
%     [theta; omega]'));

B_sym = jacobian(F,v2);
% B = double(subs(B_sym,[Phi Theta Psi omega1 omega2 omega3], ...
%     [theta; omega]'));
