%--------------------------------------------------------------------------
% DTQP_NLP_parse_objective.m
% Parse element.lagrange to obtain nonlinear objective function and
% linear/quadratic matrices (as needed)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Project Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,c,lqdoflag,in,opts] = DTQP_NLP_parse_objective(setup,in,opts)

% check if there is a X field
if isfield(setup.element,'lagrange') && ~isempty(setup.element.lagrange)

    % extract
    L = setup.element.lagrange;

    % not considered an LQDO problem anymore
    lqdoflag = false;

else % only linear-quadratic objective terms

    % construct objective function matrices
    H = DTQP_createH(setup.L,setup.M,in,opts); % create Hessian
    f = DTQP_createf(setup.l,setup.m,in,opts); % create gradient
    c = DTQP_createc(setup.cL,setup.cM,in,opts); % determine constants

    % fmincon expects H directly
    H = H/2;

    % default field value
    in.obj = [];

    % still can be considered an LQDO problem
    lqdoflag = true;

    return

end

% false because cannot detect quadratic terms symbolically yet in DTQP_NLP_symb.m
linflagOb = false;

% check internal information is already available for the objective
if isfield(in.internalinfo,'obj')

    % extract
    obj = in.internalinfo.obj;

else

    % calculate derivatives
    [obj,opts] = DTQP_NLP_symb(L,in,linflagOb,false,opts);

    % store to internal information structure
    in.internalinfo.obj = obj;

end

% check if there are OLQ elements
if isfield(obj,'Ilin')

    % check if there are any quadratic objective function terms
    if ~isempty(obj.Ilin)

        % assign
        in.IDlin = obj.Ilin;
        in.IDnon = obj.Inon;

        % TODO: construct LQDO objective terms

    end
end

% assign
in.obj = obj;

% construct objective function matrices
H = DTQP_createH(setup.L,setup.M,in,opts); % create Hessian
f = DTQP_createf(setup.l,setup.m,in,opts); % create gradient
c = DTQP_createc(setup.cL,setup.cM,in,opts); % determine constants

% fmincon expects H directly
H = H/2;

end