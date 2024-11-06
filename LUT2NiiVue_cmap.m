function LUT2NiiVue_cmap(input_fname, output_fname)
% LUT2NiiVue_cmap is a function to convert Freesurfer LUT colormap into a
% JSON file for storing colormap informations used in NiiVue.
%
% Syntax:  LUT2NiiVue_cmap(input_fname, output_fname)
%
% Inputs:
%   input_fname     - string, fullfile name of the LUT file to read from
%   output_fname    - string, fullfile name of the JSON file to be written
%
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
% Toolbox required: none
%
% See also:
% Author:
%               Nicolas Roehri, Ph.D., neuroscience
%
%
% email address:    nicolas.roehri@unige.ch, roehri.nicolas@gmail.com,
%
% Oct 2024; Last revision: 22-Oct-2024


if ~exist('input_fname', 'var') || ~exist(input_fname, 'file')
    error('Incorrect input filename variable.')
end

if ~exist('output_fname', 'var') || ~(isstring(output_fname) || ischar(output_fname))
    error('Incorrect output filename variable.')
end

folder_name = fileparts(output_fname);
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end


[code, name, rgbv] = read_LUT(input_fname);

% Lausanne atlas LUT start at 1 not zero so we have to add one element
if code(1) ~= 0
    code = [0; code];
    name = ['Unknown'; name];
    rgbv = [zeros(1,4); 15*ones(1,4); rgbv(2:end,:)];
end

% force the alpha so that 0 is transparent others aren't
A = [0; 255*ones(length(name)-1,1)];

JSON_cmap = struct('R', rgbv(:,1), 'G', rgbv(:,2), 'B', rgbv(:,3),...
    'I', code, 'A', A, 'labels', {name});


json_txt = jsonencode(JSON_cmap);
try
    fid = fopen(output_fname, 'w');
    fprintf(fid, '%s', json_txt);
catch EM
    fclose(fid);
    rethrow(EM)
end
fclose(fid);

end


function [code, name, rgbv, tt] = read_LUT(fname)
% adapted from Freesurfer read_fscolorlut(fname)
code = [];
try
    fp = fopen(fname,'r');
    if(fp == -1)
        fprintf('ERROR: could not open %s\n',fname);
        return;
    end
    
    tt = [];
    name = '';
    nthitem = 1;
    while(1)
        
        % scroll through any blank lines or comments %
        while(1)
            tline = fgetl(fp);
            if(~isempty(tline) && tline(1) == -1); break; end
            if(~isempty(deblank(tline)) && tline(1) ~= '#'); break; end
        end
        if(tline(1) == -1); break; end
        
        c = sscanf(tline,'%d',1);
        n = sscanf(tline,'%*d %s',1);
        r = sscanf(tline,'%*d %*s %d',1);
        g = sscanf(tline,'%*d %*s %*d %d',1);
        b = sscanf(tline,'%*d %*s %*d %*d %d',1);
        v = sscanf(tline,'%*d %*s %*d %*d %*d %d',1);
        t = sscanf(tline,'%*d %*s %*d %*d %*d %*d %d',1);
        code(nthitem,1) = c;
        name = strvcat(name,n');
        rgbv(nthitem,:) = [r g b v];
        if(~isempty(t)) tt(nthitem,1) = t; end
        
        nthitem = nthitem + 1;
    end
    name = cellstr(name);
catch EM
    fclose(fp);
    rethrow(EM)
end
fclose(fp);
end