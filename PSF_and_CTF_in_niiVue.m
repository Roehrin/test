%% write topography on electrode mesh for display using Niivue
addpath('E:\roehri\MATLAB\Matlab_External_Toolboxes\fieldtrip\')
output_dir = '\\129.195.132.250\users\common_resources\template_parcelation\derivatives\source_modelling\sub-template\eeg';
leadfield_name = fullfile(output_dir, 'sub-template_leadfield.mat');
load(leadfield_name, 'leadfield');
ft_defaults

method_list = {'eloreta', 'sloreta', 'mne', 'wmne', 'wmne2'};
lfd_mat = cat(2, leadfield.leadfield{:});
%% Compute Resolution matrix and PSF
for id_meth = 1:length(method_list)
    invop_name = fullfile(output_dir, ...
        sprintf('sub-template_desc-%s_inversefilters.mat', method_list{id_meth}));
    load(invop_name, 'inverse_filters');

    inv_op = cat(1, inverse_filters.filter{:});
    n_src = length(inv_op)/3;
    R = inv_op*lfd_mat;

    R_small = zeros(n_src);

    for id1 = 1:n_src
        curr_id1 = 3*(id1-1)+(1:3);
        for id2 = 1:n_src
            curr_id2 = 3*(id2-1)+(1:3);
            R_small(id1,id2) = norm(R(curr_id1,curr_id2),'fro');
        end
    end

    R_downs = R_small(:,1:50:end);
    R_downs = 1e2*R_downs./(max(R_downs,[],1));
    leadfield.PSF = zeros(length(leadfield.inside), size(R_downs,2));
    leadfield.PSF(leadfield.inside,:) = R_downs;

    cfg = [];
    cfg.filename  = fullfile('.','src',sprintf('PSF_%s',method_list{id_meth}));
    cfg.filetype  ='nifti';
    cfg.parameter = 'PSF';
    cfg.precision = 'single';
    ft_sourcewrite(cfg, leadfield);
end


%% Compute Resolution matrix and PSF for different lambdas
preproc_toolbox_path = 'E:\roehri\Matlab\preprocessing_toolbox';
addpath(genpath(preproc_toolbox_path))
lambda = [0 1e-5 1e-2 1e-1 1e0];
method_list = {'eloreta', 'sloreta', 'mne', 'wmne'};
noiseC = eye(length(leadfield.label));
inputs = struct('lambda',0.05,'method','',...
    'prewhiten','no', 'norm_lf', 'no', 'scalesourcecov','yes');
n_src = length(lfd_mat)/3;
id_src = 3*(0:250:n_src-1)+(1:3)';
id_src = id_src(:);
n_src2 = length(id_src)/3;
for id_meth = 1:length(method_list)
    leadfield.PSF = NaN(length(leadfield.inside), n_src2, length(lambda));
    for id_lmbd = 1:length(lambda)
        inputs.method = method_list{id_meth};
        inputs.lambda = lambda(id_lmbd);
        inverse_filters = create_inverse_solution(leadfield, noiseC, ...
                noiseC, inputs);

        inv_op = cat(1, inverse_filters.filter{:});
        R = inv_op*lfd_mat(:,id_src);
       
        R_small = zeros(n_src,n_src2);

        for id1 = 1:n_src
            curr_id1 = 3*(id1-1)+(1:3);
            for id2 = 1:n_src2
                curr_id2 = 3*(id2-1)+(1:3);
                R_small(id1,id2) = norm(R(curr_id1,curr_id2),'fro');
            end
        end

        R_small = 1e2*R_small./(max(R_small,[],1));
        leadfield.PSF(leadfield.inside,:,id_lmbd) = R_small;
                                               
    end
    leadfield.PSF = reshape(leadfield.PSF, length(leadfield.inside),[]);
    cfg = [];
    cfg.filename  = fullfile('.','src',sprintf('PSF_%s_varLambda',method_list{id_meth}));
    cfg.filetype  ='nifti';
    cfg.parameter = 'PSF';
    cfg.precision = 'single';
    ft_sourcewrite(cfg, leadfield);
end
