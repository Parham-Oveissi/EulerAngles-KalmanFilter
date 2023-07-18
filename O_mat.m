function [ Otheta ] = O_mat( theta, n)
%   Constructs O_n(theta). 
%   n is 1,2, or 3. 
%   theta is in degrees. 


if n==3
    Otheta      = blkdiag(cosd(theta), cosd(theta), 1);
    Otheta(1,2) = sind(theta);
    Otheta(2,1) = -sind(theta);
elseif n==2
    Otheta      = blkdiag(cosd(theta), 1, cosd(theta));
    Otheta(1,3) = -sind(theta);
    Otheta(3,1) = sind(theta);
elseif n==1
    Otheta      = blkdiag(1, cosd(theta), cosd(theta));
    Otheta(2,3) = sind(theta);
    Otheta(3,2) = -sind(theta);
end

end

