%--------------------------------------------------------------------------
% DTQP_bnds.m
% Create sequences for simple upper and lower bound terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,V] = DTQP_bnds(bnd,in)

    % extract some of the variables
    nt = in.nt; nu = in.nu; ny = in.ny;

    % find time dependent matrix (column vector)
    Bndt = DTQP_tmatrix(bnd.matrix,in.p,in.t);

    % variable locations for the variable type
    C = in.i{bnd.right};

    % initialize storage arrays
    Isav = cell(length(C),1); Vsav = Isav;

    % go through each variable of the current type
    for k = 1:length(C)

        switch bnd.right

            % control or states
            case {1,2} 
                % rows in lb/ub
                Isav{k} = DTQP_getQPIndex(C(k),bnd.right,1,nt,nu,ny);

                % nt values assigned
                if length(size(Bndt))==3
                    Vs = Bndt(:,:,k);
                else
                    Vs = Bndt(:,k,:);
                end

                % combine
                Vsav{k} = Vs;

            % parameters, initial states, or final states
            case {3,4,5,6,7}
                % row in lb/ub
                Isav{k} = DTQP_getQPIndex(C(k),bnd.right,0,nt,nu,ny);

                % single value assigned
                if length(size(Bndt))==3
                    Vs = Bndt(1,:,k);
                else
                    Vs = Bndt(1,k,:);
                end

                % combine
                Vsav{k} = Vs;

        end % end switch

    end % end for i
    
    % combine
    I = vertcat(Isav{:});
    V = vertcat(Vsav{:});

end % end function