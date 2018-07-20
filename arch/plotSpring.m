function [ vec ] = plotSpring( p1, p2, radius, nCoils )


    offset = 0.033;
    
    springSize = norm(p2-p1)-offset*2;
    teta = 0:(pi/30):2*pi*nCoils;
    x = radius*sin(teta);
    y = radius*cos(teta);
    z = teta/(2*pi*nCoils)*springSize + offset;
    
    x = [0 0 x 0 0];
    y = [0 0 y 0 0];
    z = [0 offset z offset+springSize 2*offset+springSize];
    
%     figure
%     plot3(x,y,z)
%     xlabel('x');
%     ylabel('y');
%     hold on
%     axis equal
%     grid on

    
    
    % Rotate
    direction = p2-p1;
    uInitial  = [0 0 1]';
    rotation_vector = cross(uInitial,direction);
    rotation_angle  = acos(dot(uInitial,direction)/(norm(uInitial)*norm(direction)));

    newPoint = rotatePoint([x(:),y(:),z(:)]', rotation_vector, rotation_angle, [0 0 0]');
    x(:) = newPoint(1,:); y(:) = newPoint(2,:); z(:) = newPoint(3,:);

    % Translate
    x = x + p1(1);
    y = y + p1(2);
    z = z + p1(3);

    vec = [x;y;z];
end

