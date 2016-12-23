% load_modal_params.m: Load estimates of parameters from the posterior
%     distribution (mean or mode). Else, the standard procedure gets the prior
%     mean values.
%

% Check to see if all necessary variables are in the workspace, and if not,
% immediately run initializePrograms.
necessary_vars = {'spath'};
for v = 1:length(necessary_vars)
  if ~exist(necessary_vars{v},'var')
    initializePrograms;
    break;
  end
end

if exist('spath_overwrite','var')
    spath = spath_overwrite;
    disp(['spath: ',spath]);
end

% Set infile
indataType = {'mhparaMode'};
if ~exist('infile0', 'var')
  infile.mhparaMode = [spath, 'mode_in'];
else
  infile.mhparaMode = infile0;
end

disp(['infile: ',infile.mhparaMode]);

fid.mhparaMode = fopen(infile.mhparaMode,'r');
params = fread(fid.mhparaMode,[npara,1],'single');
params = params.*(1-para_mask)+para_fix.*para_mask;
fclose(fid.mhparaMode);

para = params;
getPara_script;
