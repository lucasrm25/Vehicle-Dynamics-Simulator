% Derivate matrix 'f_in' in respect to each cell element of 'derivator'
%
% 'derivator' is a cell array containing all the derivators
% 'f_out' is a cell array with the same size of 'derivator'

function [f_out] = diff2(f_in,derivator)

    f_out = cell(size(derivator));
    syms substitute;
    for i=1:numel(derivator)
        f_aux       = subs(f_in,derivator(i),substitute);
        f_aux_diff  = diff(f_aux,substitute);
        f_out{i}    = subs(f_aux_diff,substitute,derivator(i));
    end
    if length(derivator) == 1
        f_out = f_out{1};
    end
end