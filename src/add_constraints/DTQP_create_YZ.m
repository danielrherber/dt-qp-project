%--------------------------------------------------------------------------
% DTQP_create_YZ.m
% Create matrices for general linear constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_create_YZ(YZ,p)
    % total number of constraints
    nYZ = length(YZ);

    % initialize sequences
    AeqI = []; AeqJ = []; AeqV = []; beqI = []; beqV = [];

    % go through each constraint structure
    for i = 1:nYZ
        
        % initialize as a boundary constraint        
        YZflag = 0;
        
        % go through each substructure
        for j = 1:length(YZ(i).linear)
            
            % check if any variables types are the controls or states
            if YZ(i).linear(j).right < 3
                YZflag = 1; % needs to be a path constraint
            end
            
            % check if there are any time varying matrices
            if iscell(YZ(i).linear(j).matrix)
                YZflag = 1; % needs to be a path constraint
            end
            
        end % end for j
        
        % check if the scalar term is time varying
        if iscell(YZ(i).b)
            YZflag = 1; % needs to be a path constraint
        end
        
        % current number of mixed equality constraints + 1
        N = length(beqI);
        
        % generate the sequences
        if YZflag
            % path constraint
            [AI,AJ,AV,b] = DTQP_path(YZ(i),p);
            
            % combine
            AeqI = [AeqI,AI+N];
            AeqJ = [AeqJ,AJ]; AeqV = [AeqV,AV];
            beqI = [beqI,N+(1:p.nt)]; beqV = [beqV,b];
            
        else
            % boundary constraint
            [AJ,AV,b] = DTQP_boundary(YZ(i),p);
            
            % combine
            AeqI = [AeqI,(N+1)*ones(1,length(AV))];
            AeqJ = [AeqJ,AJ]; AeqV = [AeqV,AV];
            beqI = [beqI,N+1]; beqV = [beqV,b];
            
        end

    end

    % create sparse matrices
    Aeq = sparse(AeqI,AeqJ,AeqV,length(beqI),p.nx);
    beq = sparse(beqI,1,beqV,length(beqI),1);
    
end % end function