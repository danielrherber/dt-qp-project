%--------------------------------------------------------------------------
% DTQP_getQPIndex.m
% Optimization variable index generating functions
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function I = DTQP_getQPIndex(x,xtype,Flag,nt,nu,ny)
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
        case {1,2} % controls or states
            I = ((x-1)*nt+1):(x*nt);
        case 3 % parameters
            I = (x-nu-ny) + (nu+ny)*nt;
        case {4,6} % initial controls or states 
            I = (x-1)*nt+1;
        case {5,7} % final controls or states
            I = x*nt;
    end

    % determine if continuous option is desired
    if (Flag == 1)
        if (xtype ~= 1) && (xtype ~= 2)
            I = I*ones(1,nt); % continuous option, output Nt elements
        end
    end
    
    % make column vector
    I = I(:);
end