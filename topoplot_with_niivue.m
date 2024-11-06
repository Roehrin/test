%% write topography on electrode mesh for display using Niivue
addpath('E:\roehri\MATLAB\Matlab_External_Toolboxes\fieldtrip\')
output_dir = '\\129.195.132.250\users\common_resources\template_parcelation\derivatives\source_modelling\sub-template\eeg';
leadfield_name = fullfile(output_dir, 'sub-template_leadfield.mat');
load(leadfield_name, 'leadfield');
ft_defaults

%% 3D electrode mesh
elec = ft_read_sens('\\129.195.132.250\users\common_resources\template_parcelation\sub-template\eeg\sub-template_electrodes.tsv');
elec.label(chan_idx) = [];
elec.elecpos(chan_idx,:) = [];
elec.chanpos(chan_idx,:) = [];
elec.chantype(chan_idx,:) = [];
elec.chanunit(chan_idx,:) = [];
% fig = test_leadfield(leadfield, headmodel,elec)
tri = projecttri(elec.elecpos, 'delaunay');
elec_mesh = struct('pos', elec.elecpos, 'tri', tri, 'coordsys', 'ras');
%% decimate leadfield 
leadfield_mat = cat(2, leadfield.leadfield{:});
leadfield_mat = reshape(leadfield_mat, 208, 3, []);
leadfield_mat = leadfield_mat(:,:,1:5:end);
n_src = size(leadfield_mat,3);

leadfield_mat = leadfield_mat/max(abs(leadfield_mat),[],'all');
leadfield_mat = 1e2*reshape(leadfield_mat, 208, []);
% save unnormalised leadfield (just scaled up)
ft_write_headshape('.\src\leadfield.gifti', elec_mesh, 'data', leadfield_mat,...
    'format', 'gifti')

leadfield_mat = reshape(leadfield_mat, 208, 3, []);
% fro norm of the source
% for id = 1:n_src
%     leadfield_mat(:,3*(id-1)+(1:3)) = leadfield_mat(:,3*(id-1)+(1:3))./...
%         norm(leadfield_mat(:,3*(id-1)+(1:3)), 'fro');
% end
% leadfield_mat = leadfield_mat/max(abs(leadfield_mat),[],'all');
% abs(max) norm of the source
%%save normalised leadfield (/max(src|ori))
for id = 1:3*n_src
    leadfield_mat(:,id) = leadfield_mat(:,id)./...
        max(abs(leadfield_mat(:,id)));
end
% same unnormalised leadfield (just scaled up)
leadfield_mat = 9*1e1*reshape(leadfield_mat, 208, []);
ft_write_headshape('.\src\norm_leadfield.gifti', elec_mesh, 'data', leadfield_mat,...
    'format', 'gifti')

%% 2D electrode mesh
load('GSN-HydroCel-257_withoutE.mat', 'layout')
EGI_2elec = struct('label', {layout.label(4:end-2)}, ...
    'elecpos', 100*[layout.pos(4:end-2,:), zeros(257,1)]);
[cheek_channels, chan_idx] = get_cheek_channels(EGI_2elec.label);
EGI_2elec.label(chan_idx) = [];
EGI_2elec.elecpos(chan_idx,:) = [];
% fig = test_leadfield(leadfield, headmodel,elec)
tri = projecttri(EGI_2elec.elecpos, 'delaunay');
EGI_2elec = struct('pos', EGI_2elec.elecpos, 'tri', tri, 'coordsys', 'ras');
ft_write_headshape('.\src\elecMesh2D.gifti', EGI_2elec,...
    'format', 'gifti')

%% headmesh



