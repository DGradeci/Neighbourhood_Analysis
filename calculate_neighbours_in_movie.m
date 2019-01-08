% calculate all the properties of the tracks given as input
% 
% arguments:
% tracks_GFP - any tracks for an GFP movie (optional)
% tracks_RFP - any tracks for an RFP movie (optional)
%
% usage
% [neighbours] = calculate_neighbours_in_movie(GFP, RFP)
% [neighbours] = calculate_neighbours_in_movie(GFP, []) - GFP only
% [neighbours] = calculate_neighbours_in_movie([], RFP) - RFP only


function [all_neighbours] = calculate_neighbours_in_movie(tracks_GFP, tracks_RFP)

% data_structure_answer =  questdlg('what experiment are you analysing?','experiment',...
%     '100% MDCK WT','100% MDCK Scrb_{kd}','90% MDCK wt - 10% MDCK Scrb_{kd}','50% MDCK wt - 50% MDCK Scrb_{kd}','10% MDCK wt - 90% MDCK Scrb_{kd}');
%     
%     if data_structure_answer == '100% MDCK WT'
%         distance_threshold = 100.0;
%     elseif data_structure_answer == '100% MDCK Scrb_{kd}'
%         distance_threshold = 133.0;
%     elseif data_structure_answer == '90% MDCK wt - 10% MDCK Scrb_{kd}'
%         distance_threshold = 100.0;    
%     elseif data_structure_answer == '50% MDCK wt - 50% MDCK Scrb_{kd}'
%         distance_threshold = 133.0;    
%     elseif data_structure_answer == '10% MDCK wt - 90% MDCK Scrb_{kd}'
%         distance_threshold = 133.0;    
%     end
experiment_type = {'100% MDCK WT','100% MDCK Scrb kd ','90% MDCK wt - 10% MDCK Scrb kd','50% MDCK wt - 50% MDCK Scrb kd','10% MDCK wt - 90% MDCK Scrb kd'};
[indx,kk] = listdlg('ListString',experiment_type,'SelectionMode','single','ListSize',[300,300],'Name','Experiment Seeding');
if indx==1 || indx==3
    distance_threshold = 100.0;
else
    distance_threshold = 133.0;
end

if size(tracks_RFP,2)==7
   tracks_RFP(:,8)=1;
end
if size(tracks_GFP,2)==7
    tracks_GFP(:,8)=0;
end

if isempty(tracks_GFP)
    disp('tracks_GFP is empty')
elseif isempty(tracks_RFP)
    disp('tracks_RFP is empty')
end

if isempty(tracks_GFP) && isempty(tracks_RFP)
    disp('ERROR! No tracks have been found.')
    return
end
    
% make sure that the correct cell type flags have been set
if ~isempty(tracks_GFP) && max(tracks_GFP(:,8)) ~= 0
    tracks_GFP(:,8) = 0;
end
if ~isempty(tracks_RFP) && max(tracks_RFP(:,8)) ~= 0
    tracks_RFP(:,8) = 1;
end

% ========================================================================
% OPTIONS
% ========================================================================
% threshold distance above which we do not consider cells to be touching
% distance_threshold = 100.;
% distance_threshold = 100.;

% put all of the data together
all_data = cat(1,tracks_GFP,tracks_RFP);

% now get the x,y, cellID and celltype
x_pos = all_data(:,1);
y_pos = all_data(:,2);
cell_ID = all_data(:,4);
cell_type = all_data(:,8); %correction for new data structure
frame_no = all_data(:,3);
max_frame_no = max(frame_no); % maximum frame number from the movie

% make some space to save all of the neighbour data
all_neighbours = zeros(size(tracks_GFP,1)+size(tracks_RFP,1),8);

% make some space for the distance histogram
dist_hist_bins = 32;
dist_hist = zeros(dist_hist_bins, max_frame_no);
dist_mean = zeros(max_frame_no,2);

% a counter for the size of neighbours array
c = 1;

% now loop through each frame and calculate the neighbours
for i=1:max_frame_no

    % give the user an update on progress
    if mod(i,10) == 0
        fprintf('Completed %d of %d frames... \n',i,max(frame_no));
    end
    
    % get the cell positions for this frame
    frm_idx = frame_no(:) == (i);   %logical array
    P = [];
    P(:,1) = x_pos(frm_idx);    %x frame_n
    P(:,2) = y_pos(frm_idx);    %y frame_n
    P(:,3) = cell_ID(frm_idx);   % this is the cell ID number from the tracking
    P(:,4) = cell_type(frm_idx); % this is the cell type (GFP or RFP?)

    % remove duplicates
    [P] = remove_duplicates(P);

    % check that we have enough neighbours for the analysis
    if (size(P,1) < 3)
        fprintf('Warning, less than 3 cells found in frame %d...\n',i);
        continue;
    end
    
    % calculate the neighbor matrix
    [n_matrix] = calculate_neighbours(P, distance_threshold);
    
    fprintf('%d, P: (%d, %d) n_matrix: (%d,%d)\n',i, size(P,1), size(P,2), size(n_matrix,1), size(n_matrix,2));
    
    % get the neighbours in the current frame
    [neighbours] = calculate_neighbours_in_frame(n_matrix, P, i);
    
    
    % append this to a large table for all of the data, format is:
    % [cellID, cell type, #neighbours, #GFP_neighbours, #RFP_neighbours, frame_number]
    all_neighbours(c:c+size(neighbours,1)-1,:) = neighbours;
    c = c + size(neighbours,1);

    % calculate the distances between cells
    [n_dist_matrix, mean_dist, std_dist] = calculate_neighbour_distances(P, dist_hist_bins);
    dist_hist(:,i+1) = n_dist_matrix;
    dist_mean(i+1,:) = [mean_dist std_dist];

end 

% get rid of empty space at the end
all_neighbours(c+1:end,:) = [];


% display the distance histogram
figure
hold on
imagesc(dist_hist,'YData',linspace(0,800,dist_hist_bins));
plot(dist_mean(:,1),'w-', 'LineWidth',2);
hold off
xlabel('Frame number');
ylabel('Neighbour distance (pixels)');
% xlim([1,max_frame_no]);
ylim([0,800]);
colorbar();

% % use a for loop to see the content of the cell array CC
% [V,CC] = voronoin(P);
% for i=1:length(CC) 
%     disp(CC{i})
% end


return