function [feats,feature_list] = cluster_graph_features_optimized(edge,bounds)


%1) Number of Nodes
N = length(edge);
feats{1} = N;
feature_list(1) = {'Number of Nodes'};


%2) Number of edge
E = nnz(edge);
feats{2} = E;
feature_list(2) = {'Number of edge'};


%3) Average Degree
feats{3} = E/N;
feature_list(3) = {'Average Degree'};


%%% Eccentricity calculation
% generate distance matrix
X = [bounds(:,1),bounds(:,2)];
D = squareform(pdist(X,'euclidean')) + eye(length(X));

% create sparse distance-weighted edge matrix
edge = triu(edge); % force edge to be upper triangular
sym_edge = sparse(edge | edge'); % full symmetric matrix
weighted = sym_edge .* D; 


% Calculating the longest edge in the graph
D_second = squareform(pdist(X,'euclidean')); % This is the Distance matrix with zeros on the main diagonal of the matrix
longest = max(max(D_second));
feats{4} = longest;
feature_list(4) = {'Longest Edge'};

% Calculating the sum of longest edges between every node and all the other nodes in the graph
SumOfLongest = sum(max(D_second));
feats{5} = SumOfLongest;
feature_list(5) = {'Sum of Longest Edges'};

% Calculating the median of longest edges between every node and all the other nodes in the graph
MedianOfLongest = median(max(D_second));
feats{6} = MedianOfLongest;
feature_list(6) = {'Median of Longest Edges'};

% Calculating the shortest edge in the graph
A=D_second;
A(eye(size(A))&A == 0) = inf;
shortest = min(min(A));
feats{7} = shortest;
feature_list(7) = {'Shortest Edge'};

% Calculating the sum of shortest edges between every node and all the other nodes in the graph
SumOfShortest =  sum(min(A));
feats{8} = SumOfShortest;
feature_list(8) = {'Sum of Shortest Edges'};

% Calculating the median of shortest edges between every node and all the other nodes in the graph
MedianOfShortest = median(min(A));
feats{9} = MedianOfShortest;
feature_list(9) = {'Median of Shortest Edges'};

% Calculating the average of shortest edges between every node and all the other nodes in the graph
MeanOfShortest = mean(min(A));
feats{10} = MeanOfShortest;
feature_list(10) = {'Average of Shortest Edges'};

% Calculating the average of longest edges between every node and all the other nodes in the graph
MeanOfLongest = mean(max(D_second));
feats{11} = MeanOfLongest;
feature_list(11) = {'Average of Longest Edges'};

% Calculating the Column-wise sum on the distance matrix
% Calculating the sum of distances from every node to all the other
% neigboring nodes and then get the statistics on the resulting vector

col_wise_sum = sum(D_second,1);

Max_col_wise_sum = max(col_wise_sum);
feats{12} = Max_col_wise_sum;
feature_list(12) = {'Max of Sum of Node Distances'};

Min_col_wise_sum = min(col_wise_sum);
feats{13} = Min_col_wise_sum;
feature_list(13) = {'Min of Sum of Node Distances'};

Median_col_wise_sum = median(col_wise_sum);
feats{14} = Median_col_wise_sum;
feature_list(14) = {'Median of Sum of Node Distances'};

Skewness_col_wise_sum = skewness(col_wise_sum);
feats{15} = Skewness_col_wise_sum;
feature_list(15) = {'Skewness of Sum of Node Distances'};

Kurtosis_col_wise_sum = kurtosis(col_wise_sum);
feats{16} = Kurtosis_col_wise_sum;
feature_list(16) = {'Kurtosis of Sum of Node Distances'};

Mode_col_wise_sum = mode(col_wise_sum);
feats{17} = Mode_col_wise_sum;
feature_list(17) = {'Mode of Sum of Node Distances'};

Mean_col_wise_sum = mean(col_wise_sum);
feats{18} = Mean_col_wise_sum;
feature_list(18) = {'Average of Sum of Node Distances'};


%% clustering coefficients
%[~,network] = graphconncomp(sym_edge);
G = graph(sym_edge)
network= conncomp(G)

for n=1:N
    nodes = find(network==n);
    En(nodes) = sum(sum(edge(nodes, nodes)));
    kn(nodes) = length(nodes);
end

%11) Clustering Coefficient C
% ratio beteween A: the number of edge between neighbors of node n and B:
% the number of possible edge between the neighbors of node n

Cn = 2*En ./ ( kn .* (kn-1) );
Cn( isnan(Cn) ) = 0; % account for divide by zero
feats{19} = sum(Cn)/N;
feature_list(19) = {'Clustering Coefficient C'};


% Clustering Coefficient D
Dn = 2*(kn + En) ./ ( kn .* (kn+1) );
Dn( isnan(Dn) ) = 0;
feats{20} = sum(Dn)/N;
feature_list(20) = {'Clustering Coefficient D'};


% Clustering Coefficient E
% count isolated nodes
iso_nodes = sum(kn == 1);
feats{21} = sum( Cn(kn > 1) ) / (N - iso_nodes);
feature_list(21) = {'Clustering Coefficient E'};

% Number of connected components
feats{22} = length(kn(kn>1));
feature_list(22) = {'Number of connected components'};


% Giant connected component ratio
feats{23} = max(kn) / N;
feature_list(23) = {'giant connected component ratio'};

% Average Connected Component Size
feats{24} = mean(kn(kn>1));
feature_list(24) = {'average connected component size'};


% Number / Percentage of Isolated Nodes
feats{25} = iso_nodes;
feature_list(25) = {'number isolated nodes'};

feats{26} = iso_nodes/N;
feature_list(26) = {'percentage isolated nodes'};


% Number/Percentage of End points
feats{27} = sum(kn==2);
feature_list(27) = {'number end nodes'};
feats{28} = sum(kn==2)/N;
feature_list(28) = {'percentage end nodes'};

% Edge length statistics
edge_lengths = weighted(:);
edge_lengths(edge_lengths==0) = []; % remove zero edge lengths

feats{29} = sum(edge_lengths)/length(edge_lengths); % mean edge-length
feature_list(29) = {'mean edge length'};
feats{30} = std(edge_lengths); % standard deviation
feature_list(30) = {'standard deviation edge length'};
feats{31} = skewness(edge_lengths); % skewness
feature_list(31) = {'skewness edge length'};
feats{32} = kurtosis(edge_lengths); % kurtosis
feature_list(32) = {'kurtosis edge length'};


n=length(edge);
Dedge=double(edge);


% Eigenvector Centrality- Unweighted
G=graph(Dedge,'upper');
C=centrality(G,'eigenvector');

feats{33} = (sum(C))/n; % mean Eigenvector Centrality
feature_list(33) = {'Average Eigenvector Centrality_Unweighted'};
feats{34} = max(C);
feature_list(34) = {'Max Eigenvector Centrality_Unweighted'};
feats{35} = min(C); 
feature_list(35) = {'Min Eigenvector Centrality_Unweighted'};
feats{36} = std(C); 
feature_list(36) = {'STD Eigenvector Centrality_Unweighted'};


%% To visualize the Cluster Graphs based on Centrality values
% I first need to complete the graph with the coordinates of the neuclei
% centroids:
A2= [bounds(:,1),bounds(:,2)];

%% Calculating the area and perimeter under the graph
[mm,nn]=size(bounds);
if mm>=3

 
    Ax=A2(:,1);
    Ay=A2(:,2);
    j=boundary(Ax,Ay);
    
    Xs1=Ax(j);
    Ys1=Ay(j);%return the original x,y coordinate corresponding to each index k
    Area1 = polyarea(Xs1,Ys1);%calculate the area of boundary
    feats{37}=Area1;
    feature_list(37) = {'Area of Boundry_0.5_default'};

    j2=boundary(Ax,Ay,1);

    Xs2=Ax(j2);
    Ys2=Ay(j2);%return the original x,y coordinate corresponding to each index k
    Area2 = polyarea(Xs2,Ys2);%calculate the area of boundary
    feats{38}=Area2;
    feature_list(38) = {'Area of Boundry_1_MaxShrink'};

    j3=boundary(Ax,Ay,0);
   
    Xs3=Ax(j3);
    Ys3=Ay(j3);%return the original x,y coordinate corresponding to each index k
    Area3 = polyarea(Xs3,Ys3);%calculate the area of boundary
    feats{39}=Area3;
    feature_list(39) = {'Area of Boundry_0_ConvexHull'};

    % The perimeter ratio here is the perimeter of the maxshrink over the
    % perimeter of the convex hull

    XYs2=[Xs2 Ys2];
    polyin2 = polyshape(XYs2);
    P2 = perimeter(polyin2); % Perimeter of the MaxShrink
    XYs3=[Xs3 Ys3];
    polyin3 = polyshape(XYs3);
    P3 = perimeter(polyin3); % Perimeter of the convex hull
    feats{40}=P2/P3;
    feature_list(40) = {'Perimeter_Ratio'};
else 
    feats{37}=0;
    feature_list(37) = {'Area of Boundry_0.5_default'};
    feats{38}=0;
    feature_list(38) = {'Area of Boundry_1_MaxShrink'};
    feats{39}=0;
    feature_list(39) = {'Area of Boundry_0_ConvexHull'};
    feats{40}=NaN;
    feature_list(40) = {'Perimeter_Ratio'};
end


feats{41} = sum(bounds(:,1))/N; % mean row of Centroid
feature_list(41) = {'Centroid_r'};
feats{42} = sum(bounds(:,2))/N;% mean column of Centroid
feature_list(42) = {'Centroid_c'};
