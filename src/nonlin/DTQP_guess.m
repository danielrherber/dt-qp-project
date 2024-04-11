%--------------------------------------------------------------------------
% DTQP_guess.m
% Construct a solution guess on the current grid
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function in = DTQP_guess(setup,in)

% check if a guess was provided
if isfield(setup,'guess') && ~isempty(setup.guess)

    % extract
    guess = setup.guess;

    % check different guess cases
    if size(guess.X,1) == 2

        % linear interpolation based on guess at end points
        X0 = interp1([in.t(1) in.t(end)],guess.X,in.t,'linear');

    else % arbitrary mesh provided

        % check if an interpolation method was provided
        if isfield(guess,'method')
            method = guess.method;
        else
            method = 'spline'; % default
        end

        % spline interpolation based on guess at end points
        X0 = interp1(guess.T,guess.X,in.t,method);

    end

    % extract
    X0uy = X0(:,1:(in.nu+in.ny));
    X0p = X0(1,end-in.np+1:end);

    % assign
    in.X0 = [X0uy(:);X0p(:)];

    return
end

% default guess
in.X0 = ones(in.nx,1);

end