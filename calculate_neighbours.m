% calculate the neighbouring cells given the voronoi tesselation and cull
% those that have centroids greater than a certain threshold distance apart
%
% the idea is to only calculate neighbours that are a 'sensible' distance
% apart
%
% arguments:
% P                  - matrix of the centroids and cell type data
% distance_threshold - a threshold distance to cull non-interacting cells
%
% returns:
% n - an array (num_tracks x num_tracks), where each index is either 0 (non
% contacting) or 1 (contacting)

function [n] = calculate_neighbours(P, distance_threshold)

% calculate the voronoi tesselation based on the points (and an internal
% Delaunay triangulation)


% check to see whether there are going to be any non-unique points
if any( diff (sortrows(P(:,1),1)) == 0 )
    disp('Warning! Non-unique points found');
end


try
    [~,CC] = voronoin(P(:,1:2)); 
catch
    if (size(P,1) < 3)
        disp('Warning, need at least 3 points for triangulation');
        n = [];
        return
    end
end
% make a sparse array for the data to save memory
n = uint8(zeros(size(CC,1),size(CC,1)));


% for ref_cell = 1:size(CC,1)
%     for cmp_cell = ref_cell+1:size(CC,1)
%         
%         % get the common vertices
%         common_vertices = intersect(CC{ref_cell,1},CC{cmp_cell,1});
%         
%         % if there are common vertices, proceed
%         if ~isempty(common_vertices) && length(common_vertices) > 1
%             
%             % calculate the distance between the ref and cmp cells
%             inter_cell_distance = sqrt((P(ref_cell,1)-P(cmp_cell,1)).^2+(P(ref_cell,2)-P(cmp_cell,2)).^2);
%             
%             % if it's closer than threshold distance, store it
%             if inter_cell_distance <= distance_threshold
%                 n(ref_cell,cmp_cell) = 1;
%                 n(cmp_cell,ref_cell) = 1;
%             end
%         end
%         
%     end
% end


for ref_cell = 1:size(CC,1)
    
    % now find cells close enough to compare
    cmp_cells = ref_cell+1:size(CC,1);
    inter_cell_distances = sqrt((P(ref_cell,1)-P(cmp_cells,1)).^2+(P(ref_cell,2)-P(cmp_cells,2)).^2);
    cmp_cells_cull = cmp_cells(inter_cell_distances <= distance_threshold);
    
    for cmp_cell = cmp_cells_cull
        
        % get the common vertices
        common_vertices = intersect(CC{ref_cell,1},CC{cmp_cell,1});
        
        % if there are common vertices, proceed
        if ~isempty(common_vertices) && length(common_vertices) > 1
                      
            % if it's closer than threshold distance, store it
            %if inter_cell_distance <= distance_threshold
            n(ref_cell,cmp_cell) = 1;
            n(cmp_cell,ref_cell) = 1;
            %end
        end
        
    end
end


return