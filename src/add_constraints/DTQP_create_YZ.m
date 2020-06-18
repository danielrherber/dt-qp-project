%--------------------------------------------------------------------------
% DTQP_create_YZ.m
% Create matrices for general linear constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_create_YZ(YZ,in)

% total number of constraints
nYZ = length(YZ);

% initialize storage arrays
AeqIsav = cell(nYZ,1); AeqJsav = AeqIsav; AeqVsav = AeqIsav;
beqIsav = AeqIsav; beqVsav = AeqIsav;

% initialize the number of mixed linear constraints
N = 0;

% go through each constraint structure
for k = 1:nYZ

    % initialize as a boundary constraint
    YZflag = 0;

    % go through each substructure
    for j = 1:length(YZ(k).linear)

        % check if any variables types are the controls or states
        if YZ(k).linear(j).right < 3
            YZflag = 1; % needs to be a path constraint
        end

        % check if there are any time varying matrices
        if iscell(YZ(k).linear(j).matrix)
            YZflag = 1; % needs to be a path constraint
        end

    end % end for j

    % check if the scalar term is time varying
    if iscell(YZ(k).b)
        YZflag = 1; % needs to be a path constraint
    end

    % generate the sequences
    if YZflag
        % path constraint
        [AI,AJ,AV,b] = DTQP_path(YZ(k),in);

        % combine
        AeqIsav{k} = AI+N;
        AeqJsav{k} = AJ;
        AeqVsav{k} = AV;
        beqIsav{k} = N+(1:in.nt)';
        beqVsav{k} = b;

        % update the current number of mixed linear constraints
        N = N + in.nt;
    else
        % boundary constraint
        [AJ,AV,b] = DTQP_boundary(YZ(k),in);

        % combine
        AeqIsav{k} = (N+1)*ones(length(AV),1);
        AeqJsav{k} = AJ;
        AeqVsav{k} = AV;
        beqIsav{k} = N+1;
        beqVsav{k} = b;

        % update the current number of mixed linear constraints
        N = N + 1;
    end

end

% combine
AeqI = vertcat(AeqIsav{:});
AeqJ = vertcat(AeqJsav{:});
AeqV = vertcat(AeqVsav{:});
beqI = vertcat(beqIsav{:});
beqV = vertcat(beqVsav{:});

% create sparse matrices
Aeq = sparse(AeqI,AeqJ,AeqV,N,in.nx);
beq = sparse(beqI,1,beqV,N,1);

end % end function