
% Example:
%               syms  = 2x1  sym
%               array = 2x50 double 
%
% result:       structArray 1x50 struct array with field from syms


function [ structArray ] = syms2structArray( syms, array, dim )

    cell_array = num2cell(array);
    fieldNames = cellfun(@char, sym2cell(syms),'UniformOutput', false);
    structArray = cell2struct(cell_array,fieldNames, dim);
end

