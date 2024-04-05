%--------------------------------------------------------------------------
% DTQP_standardizedinputs2.m
% Standardized input function
%--------------------------------------------------------------------------
% NOTE: only works if you use the standardized input format
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------

% check if external inputs are defined (e.g., with DTQP_sensmethods)
if exist('externalInput','var')
    n = length(externalInput);
else
    n = 0;
end

% first external input will define auxdata
if n >= 1
    auxdata = externalInput{1};
else
    auxdata = [];
end

% second external input will define opts
if n >= 2
    opts = externalInput{2};
end

% potentially set current file name and path
if ~isfield(opts,'general')
   opts.general = [];
end
if ~isfield(opts,'mname')
    stack = dbstack('-completenames');
    mname = stack(2).name; % get calling function
    % mname(end-4:end) = [];
	opts.general.mname = mname;
end
mname = opts.general.mname;
if ~isfield(opts.general,'mpath')
    mpath = mfoldername(mname,'');
    opts.general.mpath = mpath;
end