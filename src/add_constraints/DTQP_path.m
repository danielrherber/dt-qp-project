%--------------------------------------------------------------------------
% DTQP_path.m
% Create sequences for path constraint terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,J,V,b] = DTQP_path(YZ,in)

% extract some of the variables
nt = in.nt; t = in.t; ini = in.i; I_stored = in.I_stored; auxdata = in.auxdata;

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

% row locations
Ir = (1:nt)';

% go through each substructure
for j = 1:length(YZ.linear) % loop through the extended variables

    % find time dependent matrix
    YZt = DTQP_tmultiprod(YZ.linear(j).matrix,auxdata,t);

    % variable locations for the variable type
    C = ini{YZ.linear(j).right};

    % go through each variable of the current type
    for i = 1:length(C)

        % row locations
        Is = Ir;

        % column locations
        Js = DTQP_getQPIndex(C(i),YZ.linear(j).right,1,nt,I_stored);

        % nt values assigned
        if length(size(YZt))==3
            Vs = YZt(:,:,i);
        else
            Vs = YZt(:,i,:);
        end

        % remove zeros
        ZeroIndex = find(~Vs);
        Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

        % combine
        Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

    end % end for i

end % end for j

% combine
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
V = vertcat(Vsav{:});

% find time dependent vector
b = DTQP_tmultiprod(YZ.b,auxdata,t);

end % end function