function edges = constructGlobalGraph(bounds,t)
% t is a threshold for distance between centroids (edge threshold)
X = [bounds(:,1),bounds(:,2)];

%distance matrix
D = pdist(X,'euclidean');

%define edges
edges = triu(true(length(bounds)), 1) & squareform(t > D);% edges shorter than the threshold is accepted only (very long edges will be removed)

