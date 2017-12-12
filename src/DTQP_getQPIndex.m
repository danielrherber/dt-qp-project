%--------------------------------------------------------------------------
% DTQP_getQPIndex.m
% Optimization variable index generating functions
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function I = DTQP_getQPIndex(x,xtype,Flag,p)
    % check if the inputs are allowed
    if (Flag == 0)
        if (xtype == 1) || (xtype == 2)
            error('getQPIndex::not allowed with Flag = 0')
        end
    end

    % get indices for the variable type 
    switch xtype
        case 0 % singleton dimension
            I = 1;
        case {1,2} % control or states
            I = ((x-1)*p.nt+1):(x*p.nt);
        case 3 % parameters
            I = (x-p.nu-p.ns) + (p.nu+p.ns)*p.nt;
        case 4 % initial states
            I = (x-1)*p.nt+1;
        case 5 % final states
            I = x*p.nt;
    end

    % determine if continuous option is desired
    if (Flag == 1)
        if (xtype ~= 1) && (xtype ~= 2)
            I = I*ones(1,p.nt); % continuous option, output Nt elements
        end
    end
end