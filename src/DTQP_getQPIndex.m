%--------------------------------------------------------------------------
% DTQP_getQPIndex.m
% Optimization variable index generating functions
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function I = DTQP_getQPIndex(x,xtype,Flag,nt,I_stored)

% get indices for the variable type
if (xtype == 1) || (xtype == 2) % controls or states
    I = I_stored(:,x);
    return
elseif (xtype == 5) || (xtype == 7) % final states or controls
    I = I_stored(end,x);
elseif (xtype == 0) % singleton dimension
    I = 1;
else % parameters, initial states or controls
    I = I_stored(1,x);
end

% determine if continuous option is desired
if Flag && (xtype ~= 1) && (xtype ~= 2)
    I = repmat(I,nt,1); % continuous option, output nt elements
end

end