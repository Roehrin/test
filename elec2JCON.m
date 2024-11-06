function JCON_struct = elec2JCON(filename, elec, conn_mat, node_type, JCON_init)
% elec2JCON is a function to convert Fieldtrip electrode/sensor structure
% into a JCON (a JSON file for storing connectivity informations) used in
% NiiVue
%
% Syntax:  elec2JCON(elec, conn_mat)
%
% Inputs:
%   filename      	- string, fullfile name of the JCON file to write
%   elec            - Fieldtrip structure of electrodes/sensors
%   conn_mat     	- NxN matrix, containing a connectivity matrix
%   node_type       - string, measure which should be encoded in the node
%                   size and colour (e.g. strength, clustering coeff., see
%                   BCT toolbox)
%   JCON_init       - structure, containing JSON fields that the user wants
%                   to modify from the default one (e.g., 'name', or 
%                   'edgeScale', 'nodeColormap').
% Outputs:
%
% Example: 
%
%      
%   OR
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Toolbox required: Brain connectivity toolbox
%
% See also:
% Author: 
%               Nicolas Roehri, Ph.D., neuroscience
%
%
% email address:    nicolas.roehri@unige.ch, roehri.nicolas@gmail.com, 
%
% Oct 2024; Last revision: 22-Oct-2024
if ~exist('filename', 'var') || ~(isstring(filename) || ischar(filename))
    error('Incorrect filename variable')
end

folder_name = fileparts(filename);
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
    
if ~exist('elec', 'var') || ~isstruct(elec) || ~isfield(elec, 'label') || ...
        (~isfield(elec, 'elecpos') && ~isfield(elec, 'chanpos'))
    error('Incorrect elec variable.')
end
n_channels = length(elec.label);

if isfield(elec, 'elecpos')
    elecpos = elec.elecpos;
else
    elecpos = elec.chanpos;
end

if ~exist('conn_mat', 'var') || isempty(conn_mat)
    conn_mat = zeros(n_channels);
elseif ~ismatrix(conn_mat) || size(conn_mat,1) ~= size(conn_mat,2)
    error('conn_mat must be a square matrix.')
end

if ~exist('JCON_init', 'var')
    JCON_init = struct();
elseif ~isstruct(JCON_init)
    error('JCON_init must be a structure.')
end
% amke it positive and symmetric
conn_mat = abs(conn_mat).*~eye(n_channels);
conn_mat = .5*(conn_mat+conn_mat');

measure_types = {'strength', 'clustcoeff', 'centrality', 'btwcentrality',...
    'localEff', 'shortestpath', 'nodaleff'};
if ~exist('node_type', 'var') || isempty(node_type)
    node_type = 'strength';
elseif ~ismember(node_type, measure_types)
    error('node_type must be an option of this list: %s.', ...
        strjoin(measure_types, ', '))
end

% create main JCON struct
JCON_struct = struct("name", "connectome",...
	"nodeColormap", "inferno",...
	"nodeColormapNegative", "inferno",...
	"nodeMinColor", 0,...
	"nodeMaxColor", 10,...
	"nodeScale", 3,...
	"edgeColormap", "warm",...
	"edgeColormapNegative", "winter",...
	"edgeMin", 0,...
	"edgeMax", 2,...
	"edgeScale", 2,"nodes", struct(), "edges", struct());

init_fldname = fieldnames(JCON_init);
if ~isempty(init_fldname)
    JCON_fldnms = fieldnames(JCON_struct);
    [bln] = ismember(init_fldname, JCON_fldnms);
    for id = find(bln)'
        JCON_struct.(init_fldname{id}) = JCON_init.(init_fldname{id});
    end
end

% compute network measure for node size
switch node_type
    case measure_types{1}
        ntwrk_meas = sum(conn_mat,2);
    case measure_types{2}
        ntwrk_meas = clustering_coef_wu(conn_mat);
    case measure_types{3}
        ntwrk_meas = eigenvector_centrality_und(conn_mat);
    case measure_types{4}
        W = 1./conn_mat;
        ntwrk_meas = betweenness_wei(W);
    case measure_types{5}
        ntwrk_meas = efficiency_wei(conn_mat,2);
    case measure_types{6}
        D = 1./conn_mat;
        D(isinf(D)) = 0;
        ntwrk_meas = sum(distance_wei(D),2)/(n_channels-1);
    case measure_types{7}
        D = 1./conn_mat;
        D(isinf(D)) = 0;
        ntwrk_meas = distance_wei(D);
        ntwrk_meas = 1./ntwrk_meas;
        ntwrk_meas(isinf(ntwrk_meas)) = 0;
        ntwrk_meas = sum(ntwrk_meas,2)/(n_channels-1);
    otherwise
        error('unknown metric')
end
% normalise ntwrk_meas to avoid oversized nodes (should be scalable in GUI)
M = max(ntwrk_meas);
if M == 0
    ntwrk_meas = ones(1,n_channels);
else
    ntwrk_meas = ntwrk_meas./max(ntwrk_meas);
end

nodes_struct = struct("names", {elec.label'}, "X", elecpos(:,1)',...
		"Y",elecpos(:,2)', "Z", elecpos(:,3)',...
		"Color", ntwrk_meas, "Size", ntwrk_meas);

edge_struct = conn_mat.*tril(ones(n_channels),-1);
if M ~= 0 % if not null mat
    edge_struct = edge_struct./(max(edge_struct(:)));
end


JCON_struct.nodes = nodes_struct;
JCON_struct.nodeMinColor = 0;
JCON_struct.nodeMaxColor = 1;
if all(edge_struct == 0)
    JCON_struct.edges = reshape(eye(n_channels),1,[]);
    JCON_struct.edgeMin = 0;
    JCON_struct.edgeMax = 0;
else
    JCON_struct.edges = reshape(edge_struct+eye(n_channels),1,[]);
    JCON_struct.edgeMin = 0;
    JCON_struct.edgeMax = 1;
end

json_txt = jsonencode(JCON_struct);
try
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', json_txt);
catch EM
    warning(EM.message)
end
fclose(fid);

