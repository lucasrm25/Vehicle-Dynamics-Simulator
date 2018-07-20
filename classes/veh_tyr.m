classdef veh_tyr < matlab.mixin.Copyable  
    
    properties
        unloadedRadius    % [m]
        aspectRatio       % -
        width             % [m]
        compressionLength % [m]
        stiffness         % [N/m] 
    end
    
    methods
        function obj = veh_tyr(unloadedRadius, aspectRatio, width, compressionLength, stiffness)
            obj.unloadedRadius      = unloadedRadius;
            obj.aspectRatio         = aspectRatio;
            obj.width               = width;
            obj.compressionLength   = compressionLength;
            obj.stiffness           = stiffness;
        end
        
        function vec = plot(obj, uWH, cWH)
            N_TETA = 30;
            N_PHI = 30;


            uInitial = [0 0 1]';
            cInitial = [0 0 0]';

            % Calculate elipe
            a = obj.width/2;
            b = obj.aspectRatio/100*obj.width /2;
            c = obj.unloadedRadius-b;

            theta  = linspace(-pi, pi, N_TETA)   ; % Poloidal angle
            phi    = linspace(0, 2*pi, N_PHI) ; % Toroidal angle


            [t, p] = meshgrid(phi, theta);

            x = (c + b.*cos(p)) .* cos(t);
            y = (c + b.*cos(p)) .* sin(t);
            z = a.*sin(p);


            % Rotate
            rotation_vector = cross(uInitial,uWH);
            rotation_angle  = acos(dot(uInitial,uWH)/(norm(uInitial)*norm(uWH)));

            newPoint = rotatePoint([x(:),y(:),z(:)]', rotation_vector, rotation_angle, cInitial);
            x(:) = newPoint(1,:); y(:) = newPoint(2,:); z(:) = newPoint(3,:);


            % Translate
            x = x + cWH(1);
            y = y + cWH(2);
            z = z + cWH(3);

            vec = cat(3,x,y,z);
        end
    end
    
end

