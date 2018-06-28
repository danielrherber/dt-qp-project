%--------------------------------------------------------------------------
% LQRScalarTransfer_output.m
% Output function for LQRScalarTransfer example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = LQRScalarTransfer_output(T,U,Y,P,F,in,opts)

% extract parameter structure
p = in.p;

% solution on T
sol(1).T = T;
sol(1).U = real(LQRScalarTransfer_U(p.a,p.b,p.c,p.d,in.tf,p.q,p.r,T));
sol(1).Y = real(LQRScalarTransfer_Y(p.a,p.b,p.c,p.d,in.tf,p.q,p.r,T));
sol(1).F = real(LQRScalarTransfer_F(p.a,p.b,p.c,p.d,in.tf,p.q,p.r));

% flag to decide if the symbolic solution is needed
symflag = 0;

% check if we need the symbolic solution
if any(isnan(sol(1).U(:))) || any(isnan(sol(1).U(:))) || isnan(sol(1).F)
    [sol(1).U,sol(1).Y,sol(1).F] = ...
        LQRScalarTransfer_SymFun(p.a,p.b,p.c,p.d,in.tf,p.q,p.r,T);    
    symflag = 1;
end

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = (DTQP_nodes_CGL(999)+1)*in.tf/2;
    sol(2).F = sol(1).F;
    if symflag
        [sol(2).U,sol(2).Y,~] = ...
            LQRScalarTransfer_SymFun(p.a,p.b,p.c,p.d,in.tf,p.q,p.r,sol(2).T);
    else
        sol(2).U = real(LQRScalarTransfer_U(p.a,p.b,p.c,p.d,in.tf,p.q,p.r,sol(2).T));
        sol(2).Y = real(LQRScalarTransfer_Y(p.a,p.b,p.c,p.d,in.tf,p.q,p.r,sol(2).T));
    end
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY);
O(1).label = 'Ymax';
O(2).value = max(errorU);
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';