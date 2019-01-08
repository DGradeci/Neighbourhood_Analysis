
% Function Taking in the full tracks from the software and returning a data
% structure where we have the Neighbouhood matrix, the gfp stucture and rfp
% structure.

function [data_structure]=run_neighbourhood_data_structure(gfp,rfp)
if isempty(gfp) && isempty(rfp)
    disp('ERROR! No tracks have been found.')
    return
end
if isempty(gfp) && ~isempty(rfp)
    [all_neighbours]=calculate_neighbours_in_movie([], rfp.tracks);
end

if ~isempty(gfp) && isempty(rfp)
    [all_neighbours]=calculate_neighbours_in_movie(gfp.tracks,[]);
end

if ~isempty(gfp) && ~isempty(rfp)
    [all_neighbours]=calculate_neighbours_in_movie(gfp.tracks, rfp.tracks);
end
data_structure.N=all_neighbours;
data_structure.gfp=gfp;
data_structure.rfp=rfp;
end

%  [data_structure1]=run_neighbourhood_data_structure(gfp0,rfp0);
%  [data_structure2]=run_neighbourhood_data_structure(gfp1,rfp1);
%  [data_structure3]=run_neighbourhood_data_structure(gfp2,rfp2);
%  [data_structure4]=run_neighbourhood_data_structure(gfp3,rfp3);
%  [data_structure5]=run_neighbourhood_data_structure(gfp00,rfp00);
%  [data_structure6]=run_neighbourhood_data_structure(gfp01,rfp01);
%  [data_structure7]=run_neighbourhood_data_structure(gfp02,rfp02);
%  [data_structure8]=run_neighbourhood_data_structure(gfp03,rfp03);