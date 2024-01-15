function [boundsout] = bwnuclei2bounds(bwnuclei)

[ERegionalMinima, LRegionalMinima] = bwboundaries(bwnuclei, 'noholes'); %ERegionalMinima may be cell boundaries

PRegionalMinima = regionProperties(bwnuclei, LRegionalMinima);
if ~isempty(PRegionalMinima)
    nuclei=ERegionalMinima;
    for i=1:numel(PRegionalMinima)
        centroids_xy(i,:)=PRegionalMinima(i).Centroid;  
    end
    centroids_rc=[centroids_xy(:,2),centroids_xy(:,1)];
    boundsout=Veta2struct(nuclei, centroids_rc);
else
    boundsout(1).r=[];
    boundsout(1).c=[];
    boundsout(1).centroid_r=[];
    boundsout(1).centroid_c=[];
end
end
    

