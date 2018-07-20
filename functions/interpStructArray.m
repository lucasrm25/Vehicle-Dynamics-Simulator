function [ Vq ] = interpStructArray(structArray,XFieldName ,Xq, YFieldName, Yq, fieldNames )

    X = structArray.(XFieldName);
    Y = structArray.(YFieldName);
    for i = 1:numel(fieldNames)
        V{i}  = structArray.(fieldNames{i});
        if min(size(V{i})) == 1
            Vq(i) = interp1(X,V{i},Xq,'pchip',0);
        else
            Vq(i) = interp2(X,Y,V{i},Xq,Yq,'cubic',0);
        end
    end
    Vq = Vq';
end

