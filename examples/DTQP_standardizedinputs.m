%--------------------------------------------------------------------------
% DTQP_standardizedinputs.m
% Standardized input function
%--------------------------------------------------------------------------
% NOTE: only works if you use the standardized input format
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = DTQP_standardizedinputs(f,in)

% number of inputs
n = length(in);

% get user options for this example
[p,opts] = feval(f);

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
[mpath,mname] = fileparts(mfilename('fullpath'));
if ~isfield(opts,'mpath')
	opts.mpath = mpath;
end
if ~isfield(opts,'mname')
	opts.mname = mname;
end

end