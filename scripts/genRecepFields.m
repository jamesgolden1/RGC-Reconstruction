function [ rf, pixRF ] = genRecepFields(sta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
rf = cell(length(sta),1);

for icell = 1:length(sta)
    stacell = sta{icell}(:,26);
    meanpix = mean(stacell);
    stdpix = std(stacell);
    rf{icell} = find(stacell < meanpix - stdpix/2);
end
pixRF = cell(length(rf),1);

for i = 1:size(sta{1},1)
    listadd = [];
    for j = 1:length(rf)
        if ismember(i, rf{j})
            listadd = [listadd, j];
        end
    end
    pixRF{i} = listadd;
end

end
