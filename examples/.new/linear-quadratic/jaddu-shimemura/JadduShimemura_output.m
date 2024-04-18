%--------------------------------------------------------------------------
% JadduShimemura_output.m
% Output function for JadduShimemura example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = JadduShimemura_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% outputs for the different cases
switch auxdata.examplenum
	%----------------------------------------------------------------------
    case 0
        % solution on T
        sol(1).T = T;
        sol(1).U = JadduShimemura0_U(T);
        sol(1).Y = JadduShimemura0_Y(T);
        sol(1).F = JadduShimemura0_F();

        % solution on high resolution T
        if opts.general.plotflag
            sol(2).T = linspace(in.t0,in.tf,1e4)';
            sol(2).U = JadduShimemura0_U(sol(2).T);
            sol(2).Y = JadduShimemura0_Y(sol(2).T);
            sol(2).F = sol(1).F;
        end

        % errors
        errorU = abs(U-sol(1).U);
        errorY = abs(Y-sol(1).Y);
        errorF = abs(F-sol(1).F);

        % outputs
        O(1).value = max(errorY(:,1));
        O(1).label = 'Ymax';
        O(2).value = max(errorU(:,1));
        O(2).label = 'Umax';
        O(3).value = max(errorF);
        O(3).label = 'F';
        O(4).value = max(in.QPcreatetime);
        O(4).label = 'QPcreatetime';
        O(5).value = max(in.QPsolvetime);
        O(5).label = 'QPsolvetime';
	%----------------------------------------------------------------------
    case 1
        sol = [];

        % optimal parameters
        params = JadduShimemura1_parameters;

        % extract
        t1 = params(1); t2 = params(2); y1_t1 = params(3);

        % solution on T
        T1 = T(T <= t1); T2 = T(T > t1 & T <= t2); T3 = T(T > t2);

        % controls
        U_t1 = JadduShimemura1_U_t1(T1,t1,y1_t1);
        U_t2 = JadduShimemura1_U_t2(T2);
        U_t3 = JadduShimemura1_U_t3(T3,t1,t2,y1_t1);
        Uactual = [U_t1;U_t2;U_t3];

        % states
        Y_t1 = JadduShimemura1_Y_t1(T1,t1,y1_t1);
        Y_t2 = JadduShimemura1_Y_t2(T2,t1,y1_t1);
        Y_t3 = JadduShimemura1_Y_t3(T3,t1,t2,y1_t1);
        Yactual = [Y_t1;Y_t2;Y_t3];

        % objective
        Factual = JadduShimemura1_objective(y1_t1);

        % assign
        sol(1).T = T;
        sol(1).U = Uactual;
        sol(1).Y = Yactual;
        sol(1).F = Factual;

        % errors
        errorU = abs(U-sol(1).U);
        errorY = abs(Y-sol(1).Y);
        errorF = abs(F-sol(1).F);

        %
        TT = linspace(in.t0,in.tf,1e4)';

        % solution on T
        T1 = TT(TT <= t1); T2 = TT(TT > t1 & TT <= t2); T3 = TT(TT > t2);

        % controls
        U_t1 = JadduShimemura1_U_t1(T1,t1,y1_t1);
        U_t2 = JadduShimemura1_U_t2(T2);
        U_t3 = JadduShimemura1_U_t3(T3,t1,t2,y1_t1);
        U = [U_t1;U_t2;U_t3];

        % states
        Y_t1 = JadduShimemura1_Y_t1(T1,t1,y1_t1);
        Y_t2 = JadduShimemura1_Y_t2(T2,t1,y1_t1);
        Y_t3 = JadduShimemura1_Y_t3(T3,t1,t2,y1_t1);
        Y = [Y_t1;Y_t2;Y_t3];

        % objective
        F = JadduShimemura1_objective(y1_t1);

        % assign
        sol(2).T = TT;
        sol(2).U = U;
        sol(2).Y = Y;
        sol(2).F = F;

        % outputs
        O(1).value = max(errorY(:,1));
        O(1).label = 'Ymax';
        O(2).value = max(errorU(:,1));
        O(2).label = 'Umax';
        O(3).value = max(errorF);
        O(3).label = 'F';
        O(4).value = max(in.QPcreatetime);
        O(4).label = 'QPcreatetime';
        O(5).value = max(in.QPsolvetime);
        O(5).label = 'QPsolvetime';
	%----------------------------------------------------------------------
    case 2
        sol = [];

        % optimal parameters
        params = JadduShimemura2_parameters;

        % extract
        t1 = params(1); y2_t1 = params(2);

        % solution on T
        T1 = T(T <= t1); T2 = T(T > t1);

        % controls
        U_t1 = JadduShimemura2_U_t1(T1,t1,y2_t1);
        U_t2 = JadduShimemura2_U_t2(T2,t1,y2_t1);
        Uactual = [U_t1;U_t2];

        % states
        Y_t1 = JadduShimemura2_Y_t1(T1,t1,y2_t1);
        Y_t2 = JadduShimemura2_Y_t2(T2,t1,y2_t1);
        Yactual = [Y_t1;Y_t2];

        % objective
        Factual = JadduShimemura2_objective(t1);

        % assign
        sol(1).T = T;
        sol(1).U = Uactual;
        sol(1).Y = Yactual;
        sol(1).F = Factual;

        % errors
        errorU = abs(U-sol(1).U);
        errorY = abs(Y-sol(1).Y);
        errorF = abs(F-sol(1).F);

        %
        TT = linspace(in.t0,in.tf,1e4)';

        % solution on T
        T1 = TT(TT <= t1); T2 = TT(TT > t1);

        % controls
        U_t1 = JadduShimemura2_U_t1(T1,t1,y2_t1);
        U_t2 = JadduShimemura2_U_t2(T2,t1,y2_t1);
        Uactual = [U_t1;U_t2];

        % states
        Y_t1 = JadduShimemura2_Y_t1(T1,t1,y2_t1);
        Y_t2 = JadduShimemura2_Y_t2(T2,t1,y2_t1);
        Yactual = [Y_t1;Y_t2];

        % objective
        Factual = JadduShimemura2_objective(t1);

        % assign
        sol(2).T = TT;
        sol(2).U = Uactual;
        sol(2).Y = Yactual;
        sol(2).F = Factual;

        % outputs
        O(1).value = max(errorY(:,1));
        O(1).label = 'Ymax';
        O(2).value = max(errorU(:,1));
        O(2).label = 'Umax';
        O(3).value = max(errorF);
        O(3).label = 'F';
        O(4).value = max(in.QPcreatetime);
        O(4).label = 'QPcreatetime';
        O(5).value = max(in.QPsolvetime);
        O(5).label = 'QPsolvetime';
	%----------------------------------------------------------------------
    case 3
        sol = [];

        % optimal parameters
        params = JadduShimemura3_parameters;

        % extract
        y2_t1 = params(1);

        % solution on T
        T1 = T(T <= 0.5); T2 = T(T > 0.5);

        % controls
        U_t1 = JadduShimemura3_U_t1(T1,y2_t1);
        U_t2 = JadduShimemura3_U_t2(T2,y2_t1);
        Uactual = [U_t1;U_t2];

        % states
        Y_t1 = JadduShimemura3_Y_t1(T1,y2_t1);
        Y_t2 = JadduShimemura3_Y_t2(T2,y2_t1);
        Yactual = [Y_t1;Y_t2];

        % objective
        Factual = JadduShimemura3_objective(y2_t1);

        % assign
        sol(1).T = T;
        sol(1).U = Uactual;
        sol(1).Y = Yactual;
        sol(1).F = Factual;

        % errors
        errorU = abs(U-sol(1).U);
        errorY = abs(Y-sol(1).Y);
        errorF = abs(F-sol(1).F);

        %
        TT = linspace(in(1).t0,in(2).tf,1e4)';

        % solution on T
        T1 = TT(TT <= 0.5); T2 = TT(TT > 0.5);

        % controls
        U_t1 = JadduShimemura3_U_t1(T1,y2_t1);
        U_t2 = JadduShimemura3_U_t2(T2,y2_t1);
        Uactual = [U_t1;U_t2];

        % states
        Y_t1 = JadduShimemura3_Y_t1(T1,y2_t1);
        Y_t2 = JadduShimemura3_Y_t2(T2,y2_t1);
        Yactual = [Y_t1;Y_t2];

        % objective
        Factual = JadduShimemura3_objective(y2_t1);

        % assign
        sol(2).T = TT;
        sol(2).U = Uactual;
        sol(2).Y = Yactual;
        sol(2).F = Factual;

        % outputs
        O(1).value = max(errorY(:,1));
        O(1).label = 'Ymax';
        O(2).value = max(errorU(:,1));
        O(2).label = 'Umax';
        O(3).value = max(errorF);
        O(3).label = 'F';
        O(4).value = max(in.QPcreatetime);
        O(4).label = 'QPcreatetime';
        O(5).value = max(in.QPsolvetime);
        O(5).label = 'QPsolvetime';
end

end