function [  ] = nice_Pie( pieHandle, pie_X, show_label_pie_bigger_than )

%   color_out = {'b','r','y','g','k','c','m','b','r','y','g','k','c','m'}; % Farben für Einflussplot

%         color_grad = colormap(jet(500));
%         rgb_index = floor( 1+ rand(numel(pie_X),1)*size(color_grad,1));                               %RANDOM
%         rgb_index = ceil( 1+ (pie_X-min(pie_X))* (size(color_grad,1)-1)/(max(pie_X)-min(pie_X)) );      %JET

    color_grad = lines(numel(pie_X));
    rgb_index = 1:size(color_grad,1);

    hPatch = findobj(pieHandle,'Type','patch');
    for iColor= 1:numel(pie_X)
        color_out = color_grad(rgb_index(iColor),:);
        set(hPatch(iColor),{'FaceColor'},{color_out});
    end

    hText  = findobj(pieHandle,'Type','text');
    for iColor= 1:numel(pie_X)
        color_out = color_grad(rgb_index(iColor),:);
        set(hText(iColor),{'Color'},{color_out});
        set(hText(iColor),{'BackgroundColor'},{'w'});
        if pie_X(iColor) < show_label_pie_bigger_than
            set(hText(iColor),{'String'},{''});
        end
    end        

    Scale_text = 1.05;
    textPositions_cell = get(hText,{'Position'}); % cell array
    textPositions = cell2mat(textPositions_cell); % numeric array
    textPositions = textPositions * Scale_text; % scale position
    set(hText,{'Position'},num2cell(textPositions,2)); % set new position
end

