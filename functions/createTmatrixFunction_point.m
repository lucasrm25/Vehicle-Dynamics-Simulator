function [ mFunction ] = createTmatrixFunction_point( p1, p2, p3, p10, p20, p30, ptarget )

    v1 = p10-p30;
    v2 = p20-p30;
    v3 = cross(v1,v2);

    v1n = p1-p3;
    v2n = p2-p3;
    v3n = cross(v1n,v2n); 

    Tmatrix_c = [v1 v2 v3]\(ptarget-p30);

    mFunction = vpa( [v1n v2n v3n]*Tmatrix_c + p3 );
end

