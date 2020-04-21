%--------------------------------------------------------------------------
% DTQP_standardizedinputs.m
% Standardized input function
%--------------------------------------------------------------------------
% NOTE: only works if you use the standardized input format
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = DTQP_standardizedinputs(f,in)

% initialize
p = [];

% number of inputs
n = length(in);

% get user options for this example
opts = feval(f);

% potentially add inputs to p or opts
if n >= 1
    p = in{1};
end
if n >= 2
    opts = in{2};
end
if n > 3
    warning('too many input arguments...');
end

% potentially set current file name and path
if ~isfield(opts,'general')
   opts.general = [];
end
if ~isfield(opts,'mname')
    mname = char(f);
    mname(end-4:end) = [];
	opts.general.mname = mname;
end
mname = opts.general.mname;
if ~isfield(opts.general,'mpath')
    mpath = mfoldername(mname,'');
    opts.general.mpath = mpath;
end

end