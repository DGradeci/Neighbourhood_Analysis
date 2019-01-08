% calculate the distance distribution of all cells in frame
% returns a histogram of the data
function [n_dist, mean_dist, std_dist] = calculate_neighbour_distances(P, bins)

DT_0 = delaunayTriangulation(P(:,1:2));
edges = DT_0.edges;
edges_xy = [P(edges(:,1),1:2) P(edges(:,2),1:2)];
d = sqrt((edges_xy(:,1)-edges_xy(:,3)).^2 + (edges_xy(:,2)-edges_xy(:,4)).^2); 
[n_dist, ~] = hist(d, linspace(0,800,bins));

mean_dist = mean(d);
std_dist = std(d);

return