function s_w = f(theta, omega)


Phi = theta(1);
Theta = theta(2);
Psi = theta(3);

S = [1 sind(Phi)*tand(Theta) cosd(Phi)*tand(Theta); ...
    0 cosd(Phi) -sind(Phi); ... 
    0 sind(Phi)*secd(Theta) cosd(Phi)*secd(Theta)];

s_w = S*omega;
end