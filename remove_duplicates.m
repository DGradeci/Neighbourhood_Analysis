function [ P_clean ] = remove_duplicates( P )
% Find and remove duplicates in the set of points

[sorted_P] = sortrows(P,1);
[duplicates_X] = diff(sorted_P(:,1)) == 0;
[duplicates_Y] = diff(sorted_P(:,1)) == 0;
[duplicates] = duplicates_X & duplicates_Y;


if any( duplicates )
    fprintf('Warning! Non-unique points found: ');
    [duplicate_P] = find(duplicates);
    for i = 1:length(duplicate_P)
        fprintf('%d, ',duplicate_P(i));
    end
    fprintf('\n');
    sorted_P(duplicate_P,:) = [];
    P_clean = sorted_P;
    return
end

P_clean = P;

return

