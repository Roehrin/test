function [ROIs_elec]= get_ROI_centroid(nifti_atlas, tbl_idx2label)

I = find(nifti_atlas.anatomy);

% convert the indices into voxel coordinates
[x,y,z] = ind2sub(nifti_atlas.dim,I);
% convert voxel coordinates into coord in mm in real space
coord_mm = nifti_atlas.transform *[x,y,z,ones(length(x),1)]';

% get the centroids of the ROIs
centroid_coord = zeros(height(tbl_idx2label),3);
for id_abrv = 1:height(tbl_idx2label)
    
    bln = nifti_atlas.anatomy(I) == tbl_idx2label.index(id_abrv);
    centroid_coord(id_abrv,:) = mean(coord_mm(1:3,bln),2)';
end

ROIs_elec = struct('elecpos', centroid_coord,...
    'chanpos',centroid_coord, 'label',{tbl_idx2label.abbreviation});