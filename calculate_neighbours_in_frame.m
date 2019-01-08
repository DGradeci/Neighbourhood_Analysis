function [Neigh_Analysis] = calculate_neighbours_in_frame(neighbors, P, frame_number)

% total number of neighbours for each cell
S2 = sum(neighbors,1);

num_cells = size(neighbors,1);

% building a table analising a single cell over time. First row corresponds
% to the true_ID, second is total n.of neighbors, third is relative number
% of green neigh, fourth is number of red neigh.

% single frame calculations
Neigh_Analysis = zeros(num_cells,8);
Neigh_Analysis(:,1) = P(:,3);         % cell ID
Neigh_Analysis(:,2) = P(:,4);         % cell type
Neigh_Analysis(:,3) = S2';            % total number of neighbours
Neigh_Analysis(:,6) = frame_number;   % original movie frame number


% SLOW VERSION
% loop through each cell, find the neighbours and calculate the number of
% each type
for i=1:size(neighbors,1)
    d_centroid = 0; n_neigh = 0;
    for j=1:size(neighbors,2)
        if neighbors(i,j)==1
            if P(j,4)==0 % cell type
              Neigh_Analysis(i,4)= Neigh_Analysis(i,4)+1;
            else
              Neigh_Analysis(i,5)=Neigh_Analysis(i,5)+1;
            end
            d_centroid = d_centroid + sqrt((P(i,1)-P(j,1)).^2 + (P(i,2)-P(j,2)).^2);
            n_neigh = n_neigh + 1;
        end
    end
    
    if n_neigh > 0
        Neigh_Analysis(i,7) = d_centroid / n_neigh;
    end
end



% calculate denisty based on sum area of triangles enclosing cell
DT = delaunayTriangulation(P(:,1:2));

% from this calculate the attached triangles for each vertex
tri = vertexAttachments(DT);

for i=1:size(tri,1)
    
    num_neighbours = length(tri{i});
    
    sum_tri = 0.;
    for t=1:num_neighbours
        [this_tri_vert] = DT.ConnectivityList(tri{i}(t),:);
        [~,V] = convhull(DT.Points(this_tri_vert,:));
        sum_tri = sum_tri + V;
    end
    
    density = (num_neighbours/sum_tri);
    Neigh_Analysis(i,8) = density;
end
    
    
    
    

% % % FAST VERSION
% % % find the non-zero entries
% [i,j] = find(neighbors>0);
% GFP_neighbors = P(j,4) == 0;
% RFP_neighbors = P(j,4) == 1;
% % now sum up the neighbours
% for g=1:length(GFP_neighbors)
%     Neigh_Analysis(i(GFP_neighbors(g)),4) = Neigh_Analysis(i(GFP_neighbors(g)),4) + 1;
% end
% for r=1:length(RFP_neighbors)
%     Neigh_Analysis(i(RFP_neighbors(r)),5) = Neigh_Analysis(i(RFP_neighbors(r)),5) + 1;
% end


return