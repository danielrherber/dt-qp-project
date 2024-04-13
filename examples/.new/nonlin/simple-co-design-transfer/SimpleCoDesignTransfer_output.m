%--------------------------------------------------------------------------
% SimpleCoDesignTransfer_output.m
% Output function for SimpleCoDesignTransfer example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = SimpleCoDesignTransfer_output(T,U,Y,P,F,in,opts)

% solution
sol = []; % no exact solution

% outputs
O(1).value = F;
O(1).label = 'v';
O(2).value = opts.timer.sym;
O(2).label = 'Tsym';
O(3).value = opts.timer.create-opts.timer.sym;
O(3).label = 'Tint';
O(4).value = opts.timer.qpsolver;
O(4).label = 'Topt';
O(5).value = O(3).value + O(4).value;
O(5).label = 'T';
O(6).value = in.output.iterations;
O(6).label = 'Iter';
O(7).value = in.output.funcCount;
O(7).label = 'funcCount';

end