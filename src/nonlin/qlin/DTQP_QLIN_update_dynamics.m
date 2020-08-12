%--------------------------------------------------------------------------
% DTQP_QLIN_update_dynamics.m
% Update the dynamics in the quasilinearization process
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function setup = DTQP_QLIN_update_dynamics(setup,A,B,G,d,T,X,param)

% convert to DTQP compatible functions
if ~isempty(A)
    setup.A = DTQP_QLIN_update_tmatrix(A,T,X,param);
end

if ~isempty(B)
    setup.B = DTQP_QLIN_update_tmatrix(B,T,X,param);
end

if ~isempty(G)
    setup.G = DTQP_QLIN_update_tmatrix(G,T,X,param);
end

if ~isempty(d)
    setup.d = DTQP_QLIN_update_tmatrix(d,T,X,param);
end

end