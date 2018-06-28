%--------------------------------------------------------------------------
% DTQP_defects_ZO.m
% Create matrices for the zero-order hold method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_ZO(A,B,G,d,in,opts)

    % extract some of the variables in p   
    nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
    nt = in.nt; h = in.h; 

    % matrix form of I in the formulas
    K = kron(eye(ny),ones(nt-1,1));

    % initialize storage arrays
    Isav = {}; Jsav = {}; Vsav = {};

    %----------------------------------------------------------------------
    % calculate matrices and sequencing vectors
    %----------------------------------------------------------------------   
    % find time dependent matrices
    if isa(A,'double') % not time varying
        % initialize
        if strcmpi(opts.dt.mesh,'ED')
            Aexpm(1,:,:) = expm( A*h(1) );
            At = repmat(Aexpm,nt-1,1,1);
        else
            At = zeros(nt-1,ns,ns);
            for i = 1:nt-1
                At(i,:,:) = expm( A*h(i) );
            end
        end
    else
       error('A matrix cannot be time varying with ZOH defect method') 
    end   
    Bt = DTQP_convolution(A,B,in,opts);
    Gt = DTQP_convolution(A,G,in,opts);
    dt = DTQP_convolution(A,d,in,opts);

    Jy = DTQP_indexcolumns(nt,ny,nu); % optimization variable (column) locations
    Jys = [Jy;Jy+1]; % combine to create paired optimization variable locations

    if nu > 0
        Ju = DTQP_indexcolumns(nt,nu,0); % optimization variable (column) locations
    end

    if np > 0
        Jp = kron(nt*(nu+ny)+(1:np)', ones(nt-1,1)); % optimization variable (column) locations
    end
    %----------------------------------------------------------------------

    % defect constraint of row continuous constraints
    for i = 1:ny
        % current defect constraint row indices
        DefectIndices = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

        %------------------------------------------------------------------
        % controls
        %------------------------------------------------------------------
        if nu > 0
            % defect constraint (row) locations
            Iu = repmat(DefectIndices,nu,1);

            % extract matrices
            Bv = reshape(Bt(:,i,:),[],1);

            % theta values
            V3 = -Bv; % theta 3

            % combine
            Is = Iu; Js = Ju; Vs = V3;

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % states
        %------------------------------------------------------------------
        % if ny > 0 % there always is at least one state
            % defect constraint (row) locations
            Iy = repmat(DefectIndices,ny,1);

            % extract matrices
            Av = reshape(At(:,i,:),[],1);

            % theta values
            V1 = -Av; % theta 1
            V2 = K(:,i); % theta 2

            % combine
            Is = [Iy;Iy]; Js = Jys; Vs = [V1;V2];

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        % end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % parameters
        %------------------------------------------------------------------
        if np > 0
            % defect constraint (row) locations
            Is = repmat(DefectIndices,np,1);

            % extract matrices
            Gv = reshape(Gt(:,i,:),[],1);

            % theta values
            Vs = -Gv; % theta 5  

            % combine
            Js = Jp;

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
        %------------------------------------------------------------------
    end

    % combine
    If = vertcat(Isav{:});
    Jf = vertcat(Jsav{:});
    Vf = vertcat(Vsav{:});

	% output sparse matrix   
    Aeq = sparse(If,Jf,Vf,ny*(nt-1),nx);

    %------------------------------------------------------------------
	% disturbance
    %------------------------------------------------------------------
    if nd > 0
        % initialize storage arrays
        Isav = {}; Vsav = {};

        % defect constraint of row continuous constraints
        for i = 1:ny
            % defect constraint (row) locations
            Is = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

            % extract matrices
            dv = reshape(dt(:,i,:),[],1);

            % nu values
            Vs = dv; % nu

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Vsav{end+1} = Vs;

        end

        % combine
        If = vertcat(Isav{:});
        Vf = vertcat(Vsav{:});

        % output sparse matrix
        beq = sparse(If,1,Vf,ny*(nt-1),1); 
    else
        % output sparse matrix
        beq = sparse([],[],[],ny*(nt-1),1);
    end
    %------------------------------------------------------------------
end