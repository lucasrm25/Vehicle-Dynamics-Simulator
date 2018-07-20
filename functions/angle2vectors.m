function ang_rad = angle2vectors( v1, v2 , vn)
%     vc = cross(v1,v2);
%     acos( dot(v1,v2) / norm(v1) / norm(v2) ) * sign(dot(vc,vn))   
    ang_rad = atan((cross(v1,v2).'*vn/norm(vn)) / (v1.'*v2));
% angle2vectors([0 1 0],[0.1 1 0],[0 0 1])
end

