function [features_All_together,all_descriptions] = extract_graph_interplay_features(coords) % coords as input 

numGroups = length(coords);
features_All_together = [];
all_descriptions = [];

for i = 1:numGroups
    for j = i+1:numGroups
            for k = j+1:numGroups
                if (isempty(coords{i})== false && isempty(coords{j})== false && isempty(coords{k})== false)
                    
                    [features_All,description] = getGraphInterplayFeatures(coords{i},coords{j},coords{k});% Graph_interplay(marker1,marker2,marker3);
                    features_All_together=[features_All_together,features_All];
                    all_descriptions=[all_descriptions,description];
                
                end
            end
    end
end


% For testing: 
% marker1=coords{i};
% marker2=coords{j};
% marker3=coords{k};