function [features_All,Complete_feature_list] = getGraphInterplayFeatures(marker1,marker2,marker3)

% Here we group the markers into target graphs
% we considered 3 groups


Graph1_cells_coordinates = marker1; % Associated with Low

Graph2_cells_coordinates = marker2; % Associated with High

Graph3_cells_coordinates = marker3; % Middle coordinates


lineWidth=2.5;


%% First Feature Set
x = [Graph1_cells_coordinates(:,1);Graph2_cells_coordinates(:,1);Graph3_cells_coordinates(:,1)]; % all nuclei centroids (x coordinate)
y = [Graph1_cells_coordinates(:,2);Graph2_cells_coordinates(:,2);Graph3_cells_coordinates(:,2)]; % all nuclei centroids (y coordinate)

all_centroids_coordinates=[Graph1_cells_coordinates;Graph2_cells_coordinates;Graph3_cells_coordinates];


Graph1_cells_indices = (1:size(Graph1_cells_coordinates,1))';

Graph2_cells_indices = (size(Graph1_cells_coordinates,1)+1:size(Graph1_cells_coordinates,1)+size(Graph2_cells_coordinates,1))';

Graph3_cells_indices = (size(Graph1_cells_coordinates,1)+size(Graph2_cells_coordinates,1)+1:length(all_centroids_coordinates))';

% Get the delaunay triangulation...
del = delaunay(x,y);    

% Returns a set of triangles such that no data points are contained in any 
% triangle's circumcircle. Each row of the numt-by-3 matrix TRI defines 
% one such triangle and contains indices into the vectors X and Y. When 
% the triangles cannot be computed (such as when the original data is 
% collinear, or X is empty), an empty matrix is returned.
    

% Plotting Delaunay
%triplot(del,x,y,'Color','y');

%Building Adjacency matrix based on delaunay
for i = 1:size(del,1)
    t = [x(del(i,:)),y(del(i,:))];
    
    
    AdjMat1(del(i,1), del(i,2)) = true;
    AdjMat1(del(i,2), del(i,3)) = true;
    AdjMat1(del(i,1), del(i,3)) = true;

    % making the symetric matrix
    AdjMat1(del(i,2), del(i,1)) = true;
    AdjMat1(del(i,3), del(i,2)) = true;
    AdjMat1(del(i,3), del(i,1)) = true;
    
    % Adjacency matrix containing distances
    AdjMat2(del(i,1), del(i,2)) = sqrt( ( t(1,1) - t(2,1) )^2 + (t(1,2) - t(2,2))^2 );
    AdjMat2(del(i,2), del(i,1)) = AdjMat2(del(i,1), del(i,2)); 
    AdjMat2(del(i,2), del(i,3)) = sqrt( ( t(2,1) - t(3,1) )^2 + (t(2,2) - t(3,2))^2 );
    AdjMat2(del(i,3), del(i,2)) = AdjMat2(del(i,2), del(i,3));
    AdjMat2(del(i,1), del(i,3)) = sqrt( ( t(1,1) - t(3,1) )^2 + (t(1,2) - t(3,2))^2 );
    AdjMat2(del(i,3), del(i,1)) = AdjMat2(del(i,1), del(i,3));
end
% AdjMat1 = AdjMat1 | AdjMat1'; % to make it symmetric (old way)

% Reading through adjacency matrix to remove long edges
AdjMat3=AdjMat1;
AdjMat4=AdjMat2;



threshold= 100;% must be set manually

for j = 1:size(AdjMat2,1)
    for k = 1:size(AdjMat2,2)
        if AdjMat2(j,k) > threshold
            AdjMat3(j,k)=0;
            AdjMat4(j,k)=0;
        end
    end
end

HH=0;
HH_dist=[];
LL=0;
LL_dist=[];
MM=0;
MM_dist=[];

HL=0;
HL_dist=[];
HM=0;
HM_dist=[];
LM=0;
LM_dist=[];


% Counting relationships
for j = 1:size(AdjMat3,1)
    for k = 1:size(AdjMat3,2)
        if AdjMat4(j,k) ~=0
            if ismember(j,Graph2_cells_indices) && ismember(k,Graph2_cells_indices)
                HH=HH+1;
                HH_dist = [HH_dist;AdjMat4(j,k)];
            end
            
            if ismember(j,Graph1_cells_indices) && ismember(k,Graph1_cells_indices)
                LL=LL+1;
                LL_dist = [LL_dist;AdjMat4(j,k)];
            end
            
            if ismember(j,Graph3_cells_indices) && ismember(k,Graph3_cells_indices)
                MM=MM+1;
                MM_dist = [MM_dist;AdjMat4(j,k)];
            end     
            
            if (ismember(j,Graph2_cells_indices) && ismember(k,Graph1_cells_indices)) || (ismember(j,Graph1_cells_indices) && ismember(k,Graph2_cells_indices))
                HL=HL+1;
                HL_dist = [HL_dist;AdjMat4(j,k)];
            end
            
            if (ismember(j,Graph2_cells_indices) && ismember(k,Graph3_cells_indices)) || (ismember(j,Graph3_cells_indices) && ismember(k,Graph2_cells_indices))
                HM=HM+1;
                HM_dist = [HM_dist;AdjMat4(j,k)];
            end
            
            if (ismember(j,Graph1_cells_indices) && ismember(k,Graph3_cells_indices)) || (ismember(j,Graph3_cells_indices) && ismember(k,Graph1_cells_indices))
                LM=LM+1;
                LM_dist = [LM_dist;AdjMat4(j,k)];
            end
        end
    end
end

vfeature(1) = HH;
vfeature(2) = LL; 
vfeature(3) = MM;
vfeature(4) = HL;
vfeature(5) = HM;
vfeature(6) = LM;

vfeature(7) = mean(HH_dist);
vfeature(8) = mean(LL_dist);
vfeature(9) = mean(MM_dist);
vfeature(10) = mean(HL_dist);
vfeature(11) = mean(HM_dist);
vfeature(12) = mean(LM_dist);

vfeature(13) = std(HH_dist);
vfeature(14) = std(LL_dist);
vfeature(15) = std(MM_dist);
vfeature(16) = std(HL_dist);
vfeature(17) = std(HM_dist);
vfeature(18) = std(LM_dist);

if length(HH_dist) ~=0
    vfeature(19) = min(HH_dist);
    vfeature(25) = max(HH_dist);
    vfeature(31) = median(HH_dist);
    vfeature(37) = mode(HH_dist);
    vfeature(43) = kurtosis(HH_dist);
    vfeature(49) = skewness(HH_dist);
else 
    vfeature(19) = 0;
    vfeature(25) = 0;
    vfeature(31) = 0;
    vfeature(37) = 0;
    vfeature(43) = 0;
    vfeature(49) = 0;
end

if length(LL_dist) ~=0
    vfeature(20) = min(LL_dist);
    vfeature(26) = max(LL_dist);
    vfeature(32) = median(LL_dist);
    vfeature(38) = mode(LL_dist);
    vfeature(44) = kurtosis(LL_dist);
    vfeature(50) = skewness(LL_dist);
else 
    vfeature(20) = 0;
    vfeature(26) = 0;
    vfeature(32) = 0;
    vfeature(38) = 0;
    vfeature(44) = 0;
    vfeature(50) = 0;
end

if length(MM_dist) ~=0
    vfeature(21) = min(MM_dist);
    vfeature(27) = max(MM_dist);
    vfeature(33) = median(MM_dist);
    vfeature(39) = mode(MM_dist);
    vfeature(45) = kurtosis(MM_dist);
    vfeature(51) = skewness(MM_dist);
else 
    vfeature(21) = 0;
    vfeature(27) = 0;
    vfeature(33) = 0;
    vfeature(39) = 0;
    vfeature(45) = 0;
    vfeature(51) = 0;
end

if length(HL_dist) ~=0
    vfeature(22) = min(HL_dist);
    vfeature(28) = max(HL_dist);
    vfeature(34) = median(HL_dist);
    vfeature(40) = mode(HL_dist);
    vfeature(46) = kurtosis(HL_dist);
    vfeature(52) = skewness(HL_dist);
else
    vfeature(22) = 0;
    vfeature(28) = 0;
    vfeature(34) = 0;
    vfeature(40) = 0;
    vfeature(46) = 0;
    vfeature(52) = 0;
end
    
if length(HM_dist) ~=0
    vfeature(23) = min(HM_dist);
    vfeature(29) = max(HM_dist);
    vfeature(35) = median(HM_dist);
    vfeature(41) = mode(HM_dist);
    vfeature(47) = kurtosis(HM_dist);
    vfeature(53) = skewness(HM_dist);
else
    vfeature(23) = 0;
    vfeature(29) = 0;
    vfeature(35) = 0;
    vfeature(41) = 0;
    vfeature(47) = 0;
    vfeature(53) = 0;
end

if length(LM_dist) ~=0
    vfeature(24) = min(LM_dist);
    vfeature(30) = max(LM_dist);
    vfeature(36) = median(LM_dist);
    vfeature(42) = mode(LM_dist);
    vfeature(48) = kurtosis(LM_dist);
    vfeature(54) = skewness(LM_dist);
else
    vfeature(24) = 0;
    vfeature(30) = 0;
    vfeature(36) = 0;
    vfeature(42) = 0;
    vfeature(48) = 0;
    vfeature(54) = 0;
end


feature_list(1) = {'HH Count'}; % High group interconnection
feature_list(2) = {'LL Count'}; % Low group interconnection
feature_list(3) = {'MM Count'}; % Middle group interconnection
feature_list(4) = {'HL Count'};
feature_list(5) = {'HM Count'};
feature_list(6) = {'LM Count'};

feature_list(7) = {'Average Distance HH'};
feature_list(8) = {'Average Distance LL'};
feature_list(9) = {'Average Distance MM'};
feature_list(10) = {'Average Distance HL'};
feature_list(11) = {'Average Distance HM'};
feature_list(12) = {'Average Distance LM'};

feature_list(13) = {'STD Distance HH'};
feature_list(14) = {'STD Distance LL'};
feature_list(15) = {'STD Distance MM'};
feature_list(16) = {'STD Distance HL'};
feature_list(17) = {'STD Distance HM'};
feature_list(18) = {'STD Distance LM'};

feature_list(19) = {'Min Distance HH'};
feature_list(20) = {'Min Distance LL'};
feature_list(21) = {'Min Distance MM'};
feature_list(22) = {'Min Distance HL'};
feature_list(23) = {'Min Distance HM'};
feature_list(24) = {'Min Distance LM'};

feature_list(25) = {'Max Distance HH'};
feature_list(26) = {'Max Distance LL'};
feature_list(27) = {'Max Distance MM'};
feature_list(28) = {'Max Distance HL'};
feature_list(29) = {'Max Distance HM'};
feature_list(30) = {'Max Distance LM'};

feature_list(31) = {'Median Distance HH'};
feature_list(32) = {'Median Distance LL'};
feature_list(33) = {'Median Distance MM'};
feature_list(34) = {'Median Distance HL'};
feature_list(35) = {'Median Distance HM'};
feature_list(36) = {'Median Distance LM'};

feature_list(37) = {'Mode Distance HH'};
feature_list(38) = {'Mode Distance LL'};
feature_list(39) = {'Mode Distance MM'};
feature_list(40) = {'Mode Distance HL'};
feature_list(41) = {'Mode Distance HM'};
feature_list(42) = {'Mode Distance LM'};

feature_list(43) = {'Kurtosis Distance HH'};
feature_list(44) = {'Kurtosis Distance LL'};
feature_list(45) = {'Kurtosis Distance MM'};
feature_list(46) = {'Kurtosis Distance HL'};
feature_list(47) = {'Kurtosis Distance HM'};
feature_list(48) = {'Kurtosis Distance LM'};

feature_list(49) = {'Skewness Distance HH'};
feature_list(50) = {'Skewness Distance LL'};
feature_list(51) = {'Skewness Distance MM'};
feature_list(52) = {'Skewness Distance HL'};
feature_list(53) = {'Skewness Distance HM'};
feature_list(54) = {'Skewness Distance LM'};

clear AdjMat1 AdjMat3 AdjMat4


%% Comparing subgraphs set of features (2nd Feature set)



thresh=300; % a threshold for distance between centroids (edge threshold)
edges = constructGlobalGraph(all_centroids_coordinates,thresh); 


edges= edges | edges';


% Removing Family 1
edges1 = edges;

ind = Graph1_cells_indices;
edges1(ind,:) = 0;
edges1(:,ind) = 0;


[CCGfeatsM1,feature_list1] = cluster_graph_features_optimized(edges1,all_centroids_coordinates); 


% Removing Family 2
edges2 = edges;

ind = Graph2_cells_indices;
edges2(ind,:) = 0;
edges2(:,ind) = 0;

[CCGfeatsM2,feature_list2] = cluster_graph_features_optimized(edges2,all_centroids_coordinates); 


% Removing Family 3
edges3 = edges;

ind = Graph3_cells_indices;
edges3(ind,:) = 0;
edges3(:,ind) = 0;
    

[CCGfeatsM3,feature_list3] = cluster_graph_features_optimized(edges3,all_centroids_coordinates); 


ccM1= cell2mat(CCGfeatsM1);
ccfeatM1=full(ccM1);

ccM2= cell2mat(CCGfeatsM2);
ccfeatM2=full(ccM2);

ccM3= cell2mat(CCGfeatsM3);
ccfeatM3=full(ccM3);


%% Measuring statistics on distances in a delauny global graph (3rd Feature set)

% We already built a global delauny for the first set of features: AdjMat2 
% AdjMat2 contains the distances while AdjMat1 only shows connections

OneDArray = reshape(AdjMat2.',1,[]);
OneDArray2 = (nonzeros(OneDArray)'); % to get an array of the distances where they exist and remove zeros

MeanDist= mean(OneDArray2); %mean(mean(AdjMat2));
feats(1)= MeanDist;
feature_list4(1) = {'Average Distance'};

StdDist= std(OneDArray2);
feats(2)= StdDist;
feature_list4(2) = {'Distance STD'};

SkewDist= skewness(OneDArray2);
feats(3)= SkewDist;
feature_list4(3) = {'Distance Skewness'};

KurtDist= kurtosis(OneDArray2);
feats(4)= KurtDist;
feature_list4(4) = {'Distance Kurtosis'};

MedDist= median(OneDArray2);
feats(5)= MedDist;
feature_list4(5) = {'Median Distance'};

ModeDist= mode(OneDArray2);
feats(6)= ModeDist;
feature_list4(6) = {'Distance Mode'};

MinDist= min(OneDArray2);
feats(7)= MinDist;
feature_list4(7) = {'Min Distance'};

MaxDist= max(OneDArray2);
feats(8)= MaxDist;
feature_list4(8) = {'Max Distance'};

LenDist= size(OneDArray2,2);% Count of nonezero values in adjacency matrix
feats(9)= LenDist;
feature_list4(9) = {'No. Nonezero Distances'};

LenVertDist= size(AdjMat2,2);% Count of the vertices of the global graph
feats(10)= LenVertDist;
feature_list4(10) = {'No. vertices'};


%% 4th Feature Set (trianglur relationships)
% Get the delaunay triangulation...
Tri_count = 0;
some_edge_keep = [];
MaxEdge = 0;
MaxEdge_keep = [];
MinEdge = 0;
MinEdge_keep = [];
AvgEdge = 0;
AvgEdge_keep = [];
area_keep = [];
perimeter_keep = [];
Tri_ind = [];

for i = 1:size(del,1)
    t2 = [x(del(i,:)),y(del(i,:))];
    
    if ((ismember(del(i,1),Graph1_cells_indices) && ismember(del(i,2),Graph2_cells_indices)&& ismember(del(i,3),Graph3_cells_indices))...
            || (ismember(del(i,1),Graph3_cells_indices) && ismember(del(i,2),Graph2_cells_indices)&& ismember(del(i,3),Graph1_cells_indices))...
            || (ismember(del(i,1),Graph2_cells_indices) && ismember(del(i,2),Graph1_cells_indices)&& ismember(del(i,3),Graph3_cells_indices))...
            || (ismember(del(i,1),Graph2_cells_indices) && ismember(del(i,2),Graph3_cells_indices)&& ismember(del(i,3),Graph1_cells_indices))...
            || (ismember(del(i,1),Graph3_cells_indices) && ismember(del(i,2),Graph1_cells_indices)&& ismember(del(i,3),Graph2_cells_indices))...
            || (ismember(del(i,1),Graph1_cells_indices) && ismember(del(i,2),Graph3_cells_indices)&& ismember(del(i,3),Graph2_cells_indices)))
         
                Tri_count = Tri_count+1;
                
                Tri_ind = [Tri_ind;del(i,:)];
                
                some3 = AdjMat2(del(i,1), del(i,2)) + AdjMat2(del(i,2), del(i,3)) + AdjMat2(del(i,3), del(i,1));
               
                AvgEdge = some3/3;
                AvgEdge_keep = [AvgEdge_keep;AvgEdge];
                
                Area3 = polyarea(x(del(i,:))',y(del(i,:))');
                area_keep = [area_keep;Area3];

                polyin = polyshape(t2);
                P = perimeter(polyin); 
                perimeter_keep = [perimeter_keep;P];
                                
                MaxEdge = max( [AdjMat2(del(i,1), del(i,2)) , AdjMat2(del(i,2), del(i,3)) , AdjMat2(del(i,3), del(i,1))]);
                MaxEdge_keep = [MaxEdge_keep;MaxEdge];
                                
                MinEdge = min( [AdjMat2(del(i,1), del(i,2)) , AdjMat2(del(i,2), del(i,3)) , AdjMat2(del(i,3), del(i,1))]);
                MinEdge_keep = [MinEdge_keep;MinEdge];

    end
end


%disp('Triangles Found')

if (isempty(Tri_ind) == false)
    TriFeats(1)= Tri_count;
    feature_list5(1) = {'Triangle Count'};

    MeanArea= mean(area_keep); 
    TriFeats(2)= MeanArea;
    feature_list5(2) = {'Average Triangle Area'};

    StdArea= std(area_keep); 
    TriFeats(3)= StdArea;
    feature_list5(3) = {'STD Triangle Area'};

    MinArea= min(area_keep); 
    TriFeats(4)= MinArea;
    feature_list5(4) = {'Minimum Triangle Area'};

    MaxArea= max(area_keep); 
    TriFeats(5)= MaxArea;
    feature_list5(5) = {'Maximum Triangle Area'};

    MedianArea= median(area_keep); 
    TriFeats(6)= MedianArea;
    feature_list5(6) = {'Median Triangle Area'};

    ModeArea= mode(area_keep); 
    TriFeats(7)= ModeArea;
    feature_list5(7) = {'Mode Triangle Area'};

    KurtArea= kurtosis(area_keep); 
    TriFeats(8)= KurtArea;
    feature_list5(8) = {'Kurtosis Triangle Area'};

    SkwArea= skewness(area_keep); 
    TriFeats(9)= SkwArea;
    feature_list5(9) = {'Skewness Triangle Area'};


    % Perimeter Related:
    MeanPerim= mean(perimeter_keep);
    TriFeats(10)= MeanPerim;
    feature_list5(10) = {'Average Triangle Perimeter'};

    StdPerim= std(perimeter_keep); 
    TriFeats(11)= StdPerim;
    feature_list5(11) = {'STD Triangle Perimeter'};

    MinPerim= min(perimeter_keep); 
    TriFeats(12)= MinPerim;
    feature_list5(12) = {'Minimum Triangle Perimeter'};

    MaxPerim= max(perimeter_keep); 
    TriFeats(13)= MaxPerim;
    feature_list5(13) = {'Maximum Triangle Perimeter'};

    MedianPerim= median(perimeter_keep); 
    TriFeats(14)= MedianPerim;
    feature_list5(14) = {'Median Triangle Perimeter'};

    ModePerim= mode(perimeter_keep); 
    TriFeats(15)= ModePerim;
    feature_list5(15) = {'Mode Triangle Perimeter'};

    KurtPerim= kurtosis(perimeter_keep); 
    TriFeats(16)= KurtPerim;
    feature_list5(16) = {'Kurtosis Triangle Perimeter'};

    SkwPerim= skewness(perimeter_keep); 
    TriFeats(17)= SkwPerim;
    feature_list5(17) = {'Skewness Triangle Perimeter'};

    
    % MaxEdge Related:
    MeanMxEdg= mean(MaxEdge_keep);
    TriFeats(18)= MeanMxEdg;
    feature_list5(18) = {'Average Triangle MaxEdge'};

    StdMxEdg= std(MaxEdge_keep); 
    TriFeats(19)= StdMxEdg;
    feature_list5(19) = {'STD Triangle MaxEdge'};

    MinMxEdg= min(MaxEdge_keep); 
    TriFeats(20)= MinMxEdg;
    feature_list5(20) = {'Minimum Triangle MaxEdge'};

    MaxMxEdg= max(MaxEdge_keep); 
    TriFeats(21)= MaxMxEdg;
    feature_list5(21) = {'Maximum Triangle MaxEdge'};

    MedianMxEdg= median(MaxEdge_keep); 
    TriFeats(22)= MedianMxEdg;
    feature_list5(22) = {'Median Triangle MaxEdge'};

    ModeMxEdg= mode(MaxEdge_keep); 
    TriFeats(23)= ModeMxEdg;
    feature_list5(23) = {'Mode Triangle MaxEdge'};

    KurtMxEdg= kurtosis(MaxEdge_keep); 
    TriFeats(24)= KurtMxEdg;
    feature_list5(24) = {'Kurtosis Triangle MaxEdge'};

    SkwMxEdg= skewness(MaxEdge_keep); 
    TriFeats(25)= SkwMxEdg;
    feature_list5(25) = {'Skewness Triangle MaxEdge'};

    
    % MinEdge Related:
    MeanMinEdg= mean(MinEdge_keep);
    TriFeats(26)= MeanMinEdg;
    feature_list5(26) = {'Average Triangle MinEdge'};

    StdMinEdg= std(MinEdge_keep); 
    TriFeats(27)= StdMinEdg;
    feature_list5(27) = {'STD Triangle MinEdge'};

    MinMinEdg= min(MinEdge_keep); 
    TriFeats(28)= MinMinEdg;
    feature_list5(28) = {'Minimum Triangle MinEdge'};

    MaxMinEdg= max(MinEdge_keep); 
    TriFeats(29)= MaxMinEdg;
    feature_list5(29) = {'Maximum Triangle MinEdge'};

    MedianMinEdg= median(MinEdge_keep); 
    TriFeats(30)= MedianMinEdg;
    feature_list5(30) = {'Median Triangle MinEdge'};

    ModeMinEdg= mode(MinEdge_keep); 
    TriFeats(31)= ModeMinEdg;
    feature_list5(31) = {'Mode Triangle MinEdge'};

    KurtMinEdg= kurtosis(MinEdge_keep); 
    TriFeats(32)= KurtMinEdg;
    feature_list5(32) = {'Kurtosis Triangle MinEdge'};

    SkwMinEdg= skewness(MinEdge_keep); 
    TriFeats(33)= SkwMinEdg;
    feature_list5(33) = {'Skewness Triangle MinEdge'};
    
else 
    
    TriFeats(1)= 0;
    feature_list5(1) = {'Triangle Count'};

    TriFeats(2)= 0;
    feature_list5(2) = {'Average Triangle Area'};

    TriFeats(3)= 0;
    feature_list5(3) = {'STD Triangle Area'};

    TriFeats(4)= 0;
    feature_list5(4) = {'Minimum Triangle Area'};

    TriFeats(5)= 0;
    feature_list5(5) = {'Maximum Triangle Area'};

    TriFeats(6)= 0;
    feature_list5(6) = {'Median Triangle Area'};

    TriFeats(7)= 0;
    feature_list5(7) = {'Mode Triangle Area'};

    TriFeats(8)= 0;
    feature_list5(8) = {'Kurtosis Triangle Area'};

    TriFeats(9)= 0;
    feature_list5(9) = {'Skewness Triangle Area'};


    TriFeats(10)= 0;
    feature_list5(10) = {'Average Triangle Perimeter'};

    TriFeats(11)= 0;
    feature_list5(11) = {'STD Triangle Perimeter'};

    TriFeats(12)= 0;
    feature_list5(12) = {'Minimum Triangle Perimeter'};

    TriFeats(13)= 0;
    feature_list5(13) = {'Maximum Triangle Perimeter'};

    TriFeats(14)= 0;
    feature_list5(14) = {'Median Triangle Perimeter'};

    TriFeats(15)= 0;
    feature_list5(15) = {'Mode Triangle Perimeter'};

    TriFeats(16)= 0;
    feature_list5(16) = {'Kurtosis Triangle Perimeter'};

    TriFeats(17)= 0;
    feature_list5(17) = {'Skewness Triangle Perimeter'};


    TriFeats(18)= 0;
    feature_list5(18) = {'Average Triangle MaxEdge'};

    TriFeats(19)= 0;
    feature_list5(19) = {'STD Triangle MaxEdge'};

    TriFeats(20)= 0;
    feature_list5(20) = {'Minimum Triangle MaxEdge'};

    TriFeats(21)= 0;
    feature_list5(21) = {'Maximum Triangle MaxEdge'};

    TriFeats(22)= 0;
    feature_list5(22) = {'Median Triangle MaxEdge'};

    TriFeats(23)= 0;
    feature_list5(23) = {'Mode Triangle MaxEdge'};

    TriFeats(24)= 0;
    feature_list5(24) = {'Kurtosis Triangle MaxEdge'};

    TriFeats(25)= 0;
    feature_list5(25) = {'Skewness Triangle MaxEdge'};


    TriFeats(26)= 0;
    feature_list5(26) = {'Average Triangle MinEdge'};

    TriFeats(27)= 0;
    feature_list5(27) = {'STD Triangle MinEdge'};

    TriFeats(28)= 0;
    feature_list5(28) = {'Minimum Triangle MinEdge'};

    TriFeats(29)= 0;
    feature_list5(29) = {'Maximum Triangle MinEdge'};

    TriFeats(30)= 0;
    feature_list5(30) = {'Median Triangle MinEdge'};

    TriFeats(31)= 0;
    feature_list5(31) = {'Mode Triangle MinEdge'};

    TriFeats(32)= 0;
    feature_list5(32) = {'Kurtosis Triangle MinEdge'};

    TriFeats(33)= 0;
    feature_list5(33) = {'Skewness Triangle MinEdge'};
end
%% Combining all features names in a cell
Complete_feature_list= [feature_list, feature_list1, feature_list2, feature_list3, feature_list4, feature_list5];

% Combining all GraphInterplay features into one array
features = [vfeature';ccfeatM1';ccfeatM2';ccfeatM3';feats';TriFeats'];
features_All = features';
