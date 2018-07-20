function [ vec ] = plotTyre( uFWH, cFWH, tyr_UnloadedRadius, tyr_AspectRatio, tyr_width)


N_TETA = 30;
N_PHI = 30;


uInitial = [0 0 1]';
cInitial = [0 0 0]';
    
% Calculate elipe
a = tyr_width/2;
b = tyr_AspectRatio/100*tyr_width /2;
c = tyr_UnloadedRadius-b;

theta  = linspace(-pi, pi, N_TETA)   ; % Poloidal angle
phi    = linspace(0, 2*pi, N_PHI) ; % Toroidal angle


[t, p] = meshgrid(phi, theta);

x = (c + b.*cos(p)) .* cos(t);
y = (c + b.*cos(p)) .* sin(t);
z = a.*sin(p);


% Rotate
rotation_vector = cross(uInitial,uFWH);
rotation_angle  = acos(dot(uInitial,uFWH)/(norm(uInitial)*norm(uFWH)));

newPoint = rotatePoint([x(:),y(:),z(:)]', rotation_vector, rotation_angle, cInitial);
x(:) = newPoint(1,:); y(:) = newPoint(2,:); z(:) = newPoint(3,:);


% Translate
x = x + cFWH(1);
y = y + cFWH(2);
z = z + cFWH(3);

vec = cat(3,x,y,z);

end

