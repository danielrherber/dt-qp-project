%--------------------------------------------------------------------------
% DTQP_linkage.m
% Create sequences for linkage constraint terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [LA,Lb,RA,Rb] = DTQP_linkage(left,right,in)

% transcribe the linkage constraints for left
if ~isempty(left)

    % get the constraint structure for the left phase
    for idx = 1:length(left)
        LEFT(idx) = left(idx).left;
    end

    % get the constant value for the constraint
    for idx = 1:length(left)
        LEFT(idx).b = left(idx).b;
    end

    % create the sequences
    [LA,Lb] = DTQP_create_YZ(LEFT,in);

else

    % set as empty matrices
    LA = []; Lb = []; % probably need correct sizes

end

% transcribe the linkage constraints for right
if ~isempty(right)

    % get the constraint structure for the right phase
    for idx = 1:length(right)
        RIGHT(idx) = right(idx).right;
    end

    % get the constant value for the constraint
    for idx = 1:length(right)
        RIGHT(idx).b = 0; % always zero for right
    end

    % create the sequences
    [RA,Rb] = DTQP_create_YZ(RIGHT,in);

else

    % set as empty matrices
    RA = []; Rb = []; % probably need correct sizes

end

end % end function