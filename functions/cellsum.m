function csum = cellsum(cellarray)
    csum = 0;
    for it = 1:numel(cellarray)
        csum = csum + cellarray{it};
    end
end

