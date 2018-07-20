% figure
% hold on
% for b=0:0.1:2*pi
%     rotated = rotateVector( [10 0 0]', [1 1 1]', b );
%     scatter3( rotated(1), rotated(2), rotated(3) )
% end
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid on

function [ newVector ] = rotateVector( vector, direction, angle )
    
    v = vector;
    u = direction/norm(direction);
    b = angle;

    T = eye(3) + (1-cos(b))*(u*u.'-eye(3)) + sin(b)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
    newVector = T*v;

end