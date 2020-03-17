%--------------------------------------------------------------------------
% DTQP_qlin_update4tmatrix.m
% Convert the linearized/quadracized function to DTQP_tmatrix form
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function B = DTQP_qlin_update4tmatrix(A,T,X,param)

try
    B = cell2mat(A); % constant matrix

    % check if it is numeric
    if ~isnumeric(B)
        B = convert4tmatrix(A,T,X,param);
    end
catch
    B = convert4tmatrix(A,T,X,param);
end

end

function B = convert4tmatrix(A,T,X,param)

% initialize
B = cell(size(A));

for ix = 1:size(A,1)
    for jx = 1:size(A,2)

        % extract entry
        a = A{ix,jx};

        % if function
        if isa(a,'function_handle')
            B{ix,jx} = @(t) a(t,T,X,param); % function
        else
            B{ix,jx} = a; % value
        end

    end
end

end