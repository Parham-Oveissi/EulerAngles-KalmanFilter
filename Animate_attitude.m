function Animate_attitude(t, y, ang, pause_time)


if strcmp(ang,'rad')
    y = y*180/pi; % convert rads to degrees
end


for ii = 1:length(t)
    O_DA = O_mat(y(ii,1),1)*O_mat(y(ii,2),2)*O_mat(y(ii,3),3);
    
    % quiver3(0,0,0, O_DA(1,1), O_DA(2,1),O_DA(3,1),'r')
    % hold on
    % quiver3(0,0,0, O_DA(1,2), O_DA(2,2),O_DA(3,2),'g')
    % quiver3(0,0,0, O_DA(1,3), O_DA(2,3),O_DA(3,3),'b')

    quiver3(0,0,0, O_DA(1,1), O_DA(1,2),O_DA(1,3),'r')
    hold on
    quiver3(0,0,0, O_DA(2,1), O_DA(2,2),O_DA(2,3),'g')
    quiver3(0,0,0, O_DA(3,1), O_DA(3,2),O_DA(3,3),'b')

    view(135,45)
    hold off
    axis([-1 1 -1 1 -1 1])
    axis square
    pause(pause_time)

end


