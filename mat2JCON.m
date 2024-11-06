%% script to convert graph with node positions in a 3D spaces
addpath('E:\roehri\Matlab\Matlab_External_Toolboxes\BrainConnectivityToolbox\2019_03_03_BCT')
addpath(genpath('E:\roehri\Matlab\preprocessing_toolbox'))
load('\\129.195.132.250\Users\roehri\Ambizione_datasets_analysis\Milan_dataset\derivatives\iEEG_coregistration_on_MNI\projected_shafts_milan_gmavg.mat')
load('\\129.195.132.250\Users\roehri\from_Louise\Leakage_Matlab_codes\ROI_info.mat',...
    'connectome')
elec = all_sens_positions{2};
addpath('E:\roehri\MATLAB\Matlab_External_Toolboxes\fieldtrip\')
ft_defaults
elec2JCON('jcon_test.jcon', elec);

Tbl_idx2label = readtable(fullfile('E:\roehri\Matlab\preprocessing_toolbox\C_ROI_time-series',...
    'res-Lausanne2_dseg.tsv'), 'FileType', 'text', 'Delimiter', '\t');
ROIidx2remove = [58:62 122:126 129]; % for Lausanne 2 only
Tbl_idx2label(ROIidx2remove,:) = [];

atlas_filename = fullfile('src', 'sub-template_atlas-L2018_res-scale2_dseg.nii.gz');
if endsWith(atlas_filename, '.nii.gz')
    [atlas_filename, temp_folder] = decompress_mri(atlas_filename);
end
atlas = ft_read_mri(atlas_filename);
atlas.coordsys = 'ras'; % we suppose that the coord sys is ras
mask = ismember(atlas.anatomy, ROIidx2remove);
atlas.anatomy(mask) = 0;
if exist(temp_folder, 'dir')
    delete(atlas_filename)
    rmdir(temp_folder)
    clear temp_folder
end
[ROIs_elec]= get_ROI_centroid(atlas, Tbl_idx2label);

JCON_init = struct("name", "ROI", 'test', 5464,...
"nodeColormap", "red", 'nodeScale', 6);
jcon_struct = elec2JCON('./src/Lausanne2_ROI.jcon', ROIs_elec, [], [], JCON_init);

n_channels = length(ROIs_elec.label);
R = double(rand(n_channels*(n_channels-1)/2,1) > 0.95);
R(R > 0) = rand(1, sum(R));
conn_mat = eye(n_channels);
conn_mat(triu(true(n_channels),1)) = 2*R;
jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_edges.jcon', ROIs_elec, conn_mat, [], JCON_init);

logConn = connectome;
logConn(connectome > 0) = log10(connectome(connectome >0));

jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_SC.jcon', ROIs_elec, connectome, [], JCON_init);

jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_SC_CC.jcon', ROIs_elec, connectome,...
    'clustcoeff', JCON_init);

jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_SC_NE.jcon', ROIs_elec, connectome,...
    'nodaleff', JCON_init);

jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_SClogfiber.jcon', ROIs_elec, logConn, [], JCON_init);

jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_SClogfiber_CC.jcon', ROIs_elec, logConn,...
    'clustcoeff', JCON_init);

jcon_struct = elec2JCON('./src/Lausanne2_ROI_n_SClogfiber_NE.jcon', ROIs_elec, logConn,...
    'nodaleff', JCON_init);

figure;
edge_struct = jcon_struct.edges;
dgrp = digraph(reshape(edge_struct, n_channels, n_channels).*~eye(n_channels));
plt_dgrph = plot(dgrp);
plt_dgrph.XData = jcon_struct.nodes.X;
plt_dgrph.YData = jcon_struct.nodes.Y;
plt_dgrph.ZData = jcon_struct.nodes.Z;

%% make EGI electrodes

load('GSN-HydroCel-257_withoutE.mat', 'layout')
load("\\129.195.132.250\Users\common_resources\template_parcelation\derivatives\mri_preprocessing\sub-template\anat\sub-template_desc-head_headmodel.mat", 'headmodel')
load("\\129.195.132.250\Users\common_resources\template_parcelation\derivatives\mri_preprocessing\sub-template\anat\sub-template_desc-head_sourcemodel.mat", 'sourcemodel')

EGI_2elec = struct('label', {layout.label(4:end-2)}, ...
    'elecpos', 1e2*[layout.pos(4:end-2,:), zeros(257,1)]);
mesh = headmodel.bnd(1);
[cheek_channels, chan_idx] = get_cheek_channels(EGI_2elec.label);
EGI_2elec.label(chan_idx) = [];
EGI_2elec.elecpos(chan_idx,:) = [];

jcon_struct = elec2JCON('./src/EGI257_2Dlayout.jcon', EGI_2elec);

if 0
    elec = ft_read_sens('GSN-HydroCel-257.sfp');

    elec.label(1:3) = [];
    elec.chanpos(1:3,:) = [];
    elec.elecpos(1:3,:) = [];
    elec.elecpos = 10*elec.elecpos;
    elec.chanpos = 10*elec.chanpos;
    % apply ICP
    P = elec.elecpos';
    Q = mesh.pos';
    [Ricp, Ticp] = icp(Q,P, 15, 'Matching', 'kDtree', 'Extrapolation', true);
    Dicp = Ricp * P + repmat(Ticp, 1, length(P));
    elec.elecpos = Dicp';
    elec.chanpos = Dicp';
else
    elec = ft_read_sens('\\129.195.132.250\users\common_resources\template_parcelation\sub-template\eeg\sub-template_electrodes.tsv');
    elec.label(chan_idx) = [];
    elec.elecpos(chan_idx,:) = [];
    elec.chanpos(chan_idx,:) = [];
    elec.chantype(chan_idx,:) = [];
    elec.chanunit(chan_idx,:) = [];
end
figure
ft_plot_mesh(headmodel.bnd(3), 'facecolor',[0.9 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);
ft_plot_mesh(headmodel.bnd(1),'edgecolor','none','facecolor',0.4*[1 1 1],'facealpha', 0.3);
ft_plot_sens(elec, 'style', '.b', 'label', 'label')

jcon_struct = elec2JCON('./src/EGI257_electrodes.jcon', elec);

if 0
    elec_tbl = table(elec.label, elec.elecpos(:,1),...
        elec.elecpos(:,2),elec.elecpos(:,3),...
        'VariableNames', {'name', 'x', 'y', 'z'});
    elec_tsv_fnm = 'sub-template_electrodes.tsv';
    writetable(elec_tbl, fullfile(elec_tsv_fnm), ...
        'FileType', 'text', 'Delimiter', '\t');
end
% source model
src_elec = struct('label', {compose('src%04i', 1:sum(sourcemodel.inside))},...
    'elecpos', sourcemodel.pos(sourcemodel.inside,:));
jcon_struct = elec2JCON('./src/source_points.jcon', src_elec);

figure;
n_channels = length(jcon_struct.nodes.names);
edge_struct = jcon_struct.edges;
dgrp = digraph(reshape(edge_struct, n_channels, n_channels).*~eye(n_channels));
plt_dgrph = plot(dgrp);
plt_dgrph.XData = jcon_struct.nodes.X;
plt_dgrph.YData = jcon_struct.nodes.Y;
plt_dgrph.ZData = jcon_struct.nodes.Z;
axis vis3d