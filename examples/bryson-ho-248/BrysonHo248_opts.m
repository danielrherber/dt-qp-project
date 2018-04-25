%--------------------------------------------------------------------------
% BrysonHo248_opts.m
% User options function for BrysonHo248 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = BrysonHo248_opts

% test number
num = 1;

switch num

case 1
    opts.Quadmethod = 'CEF';
    opts.Defectmethod = 'ZO';
    opts.NType = 'ED';
    p.nt = 100;

case 2
    opts.Quadmethod = 'G';
    opts.Defectmethod = 'PS';
    opts.NType = 'LGL';
    p.nt = 100;
    
case 3      
    opts.Quadmethod = 'CQHS';
    opts.Defectmethod = 'HS';
    opts.NType = 'ED';
    p.nt = 100;
    
end

end