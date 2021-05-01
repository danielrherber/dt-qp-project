%--------------------------------------------------------------------------
% DTQP_checkCellNumeric.m
% Check if the cell matrix can be converted to a numeric matrix
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function flag = DTQP_checkCellNumeric(opt,varargin)

switch opt
    %----------------------------------------------------------------------
    case 1

    % try to convert each input
    try

        % initialize all inputs as numeric
        flag = true;

        % go through each input
        for k = 1:nargin

            %
            A = cell2mat(varargin{k});

            %
            if ~isnumeric(A)

                % some inputs are not numeric
                flag = false;
                break

            end
        end

    catch

        % some inputs are not numeric
        flag = false;

    end
    %----------------------------------------------------------------------
    case 'element'

    % extract
    A = varargin{1};

    % initialize all inputs as numeric
    flag = true(size(A));

    % go through each element in the matrix
    for k = 1:numel(A)

        % check if the element is numeric
        flag(k) = isnumeric(A{k});

    end
    %----------------------------------------------------------------------
end

end