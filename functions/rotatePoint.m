% figure
% hold on
% for b=0:0.1:2*pi
%     rotated = rotatePoint( [10 0 0]', [0 1 0]', b, [9 0 0]' );
%     scatter3( rotated(1), rotated(2), rotated(3) )
% end
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid on

function [ newPoint ] = rotatePoint( point, direction, angle, center )

    p = point;
    u = direction/norm(direction);
    b = angle; 
    c = repmat(center,1,size(p,2));
  
    T = eye(3) + (1-cos(b))*(u*u.'-eye(3)) + sin(b)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];

    relativ = p - c;
    newPoint = T*relativ + c;

end
