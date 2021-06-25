%--------------------------------------------------------------------------
% DTQP_MESH_pts.m
% Generate the time mesh (vector of discrete time values). Also,
% potentially generate the quadrature weights and differentiation matrix
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [t,w,D] = DTQP_MESH_pts(in,dt)

% extract
if isfield(in,'t0')
    t0 = in.t0;
else
    t0 = 0;
end
tf = in.tf;

% initialize
w = []; % empty
D = []; % empty

switch upper(dt.mesh)
    %----------------------------------------------------------------------
    case 'ED' % equidistant node points

    t = linspace(t0,tf,dt.nt)';

    %----------------------------------------------------------------------
    case 'LGL' % Legendre-Gauss-Lobatto nodes

    % scaled nodes and quadrature weights
    if strcmpi(dt.quadrature,'G') % Gaussian quadrature
        [tau,w] = lobpts_stored(dt.nt);
    else % other
        tau = lobpts(dt.nt); % using chebfun
    end

    % differentiation matrix
    if strcmpi(dt.defects,'PS')
        D = legslbdiff_stored(dt.nt,tau);
    else
        D = []; % empty
    end

    % unscale mesh
    t = (tau + (tf+t0)/(tf-t0))*(tf-t0)*0.5;

    %----------------------------------------------------------------------
    case 'CGL' % Chebyshev-Gauss-Lobatto nodes

    % scaled nodes and quadrature weights
    if strcmpi(dt.quadrature,'CC') % Clenshaw-Curtis quadrature
        [tau,w] = chebpts(dt.nt); % using chebfun
    else % other
        tau = chebpts(dt.nt); % using chebfun
    end

    % differentiation matrix
    if strcmpi(dt.defects,'PS')
        D = diffmat(dt.nt,'chebkind2'); % using chebfun
    else
        D = []; % empty
    end

    % unscale mesh
    t = (tau + (tf+t0)/(tf-t0))*(tf-t0)*0.5;

    %----------------------------------------------------------------------
    case 'USER' % user-defined nodes

    if ~isfield(dt,'t')
        error('ERROR: opts.dt.t does not exist with USER option specified')
    else
        t = dt.t(:);
    end

    if strcmpi(dt.defects,'PS')
        error('PS option cannot handle USER mesh')
    end

    %----------------------------------------------------------------------
end

end

% determine Legendre-Gauss-Lobatto scaled nodes and quadrature weights
% (with look up for small nt)
function [tau,w] = lobpts_stored(nt)

% check if the values are precomputed below
if nt > 10
    [tau,w] = lobpts(nt); % using chebfun
else

    % NOTE: these were manually calibrated to ensure zero error between
    % lobpts and the decimal representations
    switch nt
        case 2
            tau = [-1.0; 1.0];
            w = [1.0, 1.0];
        case 3
            tau = [-1.0; 0; 1.0];
            w = [0.33333333333333333333333333333333, 1.3333333333333333333333333333333, 0.33333333333333333333333333333333];
        case 4
            tau = [-1.0; -0.44721359549995793928183473374626; 0.44721359549995793928183473374626; 1.0];
            w = [0.16666666666666666666666666666667, 0.8333333333333333, 0.8333333333333333, 0.16666666666666666666666666666667];
        case 5
            tau = [-1.0; -0.65465367070797714379829245624686; 0; 0.65465367070797714379829245624686; 1.0];
            w = [0.1, 0.5444444444444445, 0.711111111111111111111111111111, 0.5444444444444445, 0.1];
        case 6
            tau = [-1.0; -0.76505532392946462572069776797434;-0.28523151648064509755542417224206; 0.28523151648064509755542417224206; 0.76505532392946462572069776797434; 1.0];
            w = [0.066666666666666666666666666666667, 0.37847495629784705384324183796707, 0.55485837703548657184882131332415, 0.55485837703548657184882131332415, 0.37847495629784705384324183796707, 0.066666666666666666666666666666667];
        case 7
            tau = [-1.0; -0.83022389627856685301310335489688; -0.46884879347071423127957245924335; 0; 0.46884879347071423127957245924335; 0.83022389627856685301310335489688; 1.0];
            w = [0.047619047619047619047619047619048, 0.27682604736156646296763028658461, 0.43174538120986305500537127954885, 0.4876190476190479, 0.43174538120986305500537127954885, 0.27682604736156646296763028658461, 0.047619047619047619047619047619048];
        case 8
            tau = [-1.0; -0.87174014850960657163625455723377; -0.59170018143314229153162386865006; -0.20929921790247887902758350264776; 0.20929921790247887902758350264776; 0.59170018143314229153162386865006; 0.87174014850960657163625455723377; 1.0];
            w = [0.035714285714285714285714285714286, 0.21070422714350597881427518132114, 0.34112269248350451933404769988556, 0.41245879465870388669301860318228, 0.41245879465870388669301860318228, 0.34112269248350451933404769988556, 0.21070422714350597881427518132114, 0.035714285714285714285714285714286];
        case 9
            tau = [-1.0; -0.8997579954114600653269917529542; -0.67718627951073773196810634544818; -0.36311746382617815509519232364255; 0; 0.36311746382617815509519232364255; 0.67718627951073773196810634544818; 0.8997579954114600653269917529542; 1.0];
            w = [0.027777777777777777777777777777778, 0.165495361560805853695654832336, 0.27453871250016137484095679610618, 0.34642851097304622198791435039311, 0.3715192743764169, 0.34642851097304622198791435039311, 0.27453871250016137484095679610618, 0.165495361560805853695654832336, 0.027777777777777777777777777777778];
        case 10
            tau = [-1.0; -0.91953390816645885763591650174931; -0.7387738651055050231875043209584; -0.47792494981044447710516465122055; -0.16527895766638703300976942500711; 0.16527895766638703300976942500711; 0.47792494981044447710516465122055; 0.7387738651055050231875043209584; 0.91953390816645885763591650174931; 1.0];
            w = [0.022222222222222222222222222222222, 0.13330599085106978329839932939649, 0.22488934206312646835179691606754, 0.29204268367968377884125175114605, 0.32753976118389743765746402459627, 0.32753976118389743765746402459627, 0.29204268367968377884125175114605, 0.22488934206312646835179691606754, 0.13330599085106978329839932939649, 0.022222222222222222222222222222222];
    end

end

end

% determine Legendre-Gauss-Lobatto differentiation matrix
% (with look up for small nt)
function D = legslbdiff_stored(nt,tau)

% check if the values are precomputed below
if nt > 10
    D = legslbdiff(nt,tau); % using chebfun
else

    % NOTE: these were manually calibrated to ensure zero error between
    % lobpts and the decimal representations
    switch nt
        case 2
            D = [-0.5, Inf; ...
                Inf, 0.5];
        case 3
            D = [-1.5, 2.0, -0.5; ...
                -0.5, 0, 0.5; ...
                0.5, -2.0, 1.5];
        case 4
            D = [-3.0,  4.04508497187473726, -1.54508497187473703, 0.5; ...
                -0.80901699437494734, 0, 1.11803398874989485, -0.309016994374947396; ...
                0.309016994374947396, -1.11803398874989485, 0, 0.80901699437494734; ...
                -0.5, 1.54508497187473703, -4.04508497187473726, 3.0];
        case 5
            D = [-5.0,  6.75650248872424175, -2.66666666666666667, 1.41016417794242654, -0.5; ...
                -1.24099025303098354, 0, 1.7457431218879387, -0.76376261582597332, 0.259009746969017129; ...
                0.375, -1.3365845776954531, 0, 1.3365845776954531, -0.375; ...
                -0.259009746969017129, 0.76376261582597332, -1.7457431218879387, 0,1.24099025303098354; ...
                0.5, -1.41016417794242654, 2.66666666666666667, -6.75650248872424175, 5.0];
        case 6
            D = [-7.5, 10.1414159363196674, -4.03618727030534785, 2.24468464817616686, -1.34991331419048799, 0.5; ...
                -1.78636494833909443, 0, 2.52342677742945654, -1.15282815853592968, 0.653547507429800278, -0.237781177984231373; ...
                0.484951047853569128, -1.72125695283023306, 0, 1.75296196636786594, -0.786356672223240571, 0.269700610832038945; ...
                -0.269700610832038945, 0.786356672223240571, -1.75296196636786594, 0, 1.72125695283023306, -0.484951047853569128; ...
                0.237781177984231373, -0.653547507429800278, 1.15282815853592968, -2.52342677742945654, 0, 1.78636494833909443; ...
                -0.5, 1.34991331419048799, -2.24468464817616686, 4.03618727030534785, -10.1414159363196674, 7.5];
        case 7
            D = [ -10.5, 14.2015766029198112, -5.66898522554550688, 3.2, -2.0499648130767425, 1.31737343570243448, -0.5; ...
                -2.44292601424428923, 0.000000000000000222044604925031308, 3.45582821429428533, -1.59860668809836692, 0.961339797288711773, -0.602247179635785779, 0.226611870395445336; ...
                0.625256665515342092, -2.2158042831699718, 0, 2.26669808708599918, -1.06644190400637462, 0.61639083551757945, -0.226099400942574691; ...
                -0.3125, 0.907544471268820874, -2.00696924058875315, 0, 2.00696924058875315, -0.907544471268820874, 0.3125; ...
                0.226099400942574691, -0.61639083551757945, 1.06644190400637462, -2.26669808708599918, 0, 2.2158042831699718, -0.625256665515342092; ...
                -0.226611870395445336, 0.602247179635785779, -0.961339797288711773, 1.59860668809836692, -3.45582821429428533, 0, 2.44292601424428923; ...
                0.5, -1.31737343570243448, 2.0499648130767425, -3.2, 5.66898522554550688, -14.2015766029198112, 10.5];
        case 8
            D = [-14.0, 18.9375986071173656, -7.5692898193484881, 4.29790816426517441, -2.81018898925794858, 1.94165942554412219, -1.29768738832023223, 0.5; ...
                -3.20991570300298612, 0, 4.54358506456656563, -2.11206121431454186, 1.29423205091350124, -0.869448098331492947, 0.57356541494026414, -0.219957514771304291; ...
                0.792476681320514298, -2.80647579473643427, 0, 2.87551740597250483, -1.3727858318060282, 0.845022556506510369, -0.53703958615766112, 0.20328456890059271; ...
                -0.372150435728594853, 1.07894468879045324, -2.37818723351550609, 0, 2.38892435915823897, -1.13535801688111171, 0.661157350900311469, -0.243330712723791032; ...
                0.243330712723791032, -0.661157350900311469, 1.13535801688111171, -2.38892435915823897, 0, 2.37818723351550609, -1.07894468879045324, 0.372150435728594853; ...
                -0.20328456890059271, 0.53703958615766112, -0.845022556506510369, 1.3727858318060282, -2.87551740597250483, 0.000000000000000222044604925031308, 2.80647579473643427, -0.792476681320514298; ...
                0.219957514771304291, -0.57356541494026414, 0.869448098331492947, -1.29423205091350124, 2.11206121431454186, -4.54358506456656563, 0, 3.20991570300298612; ...
                -0.5, 1.29768738832023223,  -1.94165942554412219, 2.81018898925794858, -4.29790816426517441, 7.5692898193484881, -18.9375986071173656, 14.0];
        case 9
            D = [ -18.0, 24.349745171593046, -9.73870165721154457, 5.5449639069493788, -3.65714285714285714, 2.59074567655935484, -1.87444087344698285, 1.28483063269958864, -0.5; ...
                -4.08701370203367009, 0, 5.78680581663731175, -2.69606544031405626, 1.66522164500538494, -1.14565373845513219, 0.816756381741385651, -0.555704981283716815, 0.215654018702498951; ...
                0.985360090074507089, -3.48835875343445689, 0, 3.57668094012561522, -1.71783215719506299, 1.0798038112826307, -0.738349277190386122, 0.492350938315507636, -0.189655591978356464; ...
                -0.444613449281090645, 1.28796075006390698, -2.83445891207942058, 0, 2.85191596846289519, -1.37696489376051212, 0.855726185092675284, -0.547300160534051505, 0.207734512035597202; ...
                0.2734375, -0.741782397916254554, 1.26941308635814942, -2.6593102175739185, 0, 2.6593102175739185, -1.26941308635814942, 0.741782397916254554, -0.2734375; ...
                -0.207734512035597202, 0.547300160534051505, -0.855726185092675284, 1.37696489376051212, -2.85191596846289519, 0, 2.83445891207942058, -1.28796075006390698, 0.444613449281090645; ...
                0.189655591978356464, -0.492350938315507636, 0.738349277190386122, -1.0798038112826307, 1.71783215719506299, -3.57668094012561522, 0, 3.48835875343445689, -0.985360090074507089; ...
                -0.215654018702498951, 0.555704981283716815, -0.816756381741385651, 1.14565373845513219, -1.66522164500538494, 2.69606544031405626, -5.78680581663731175, 0, 4.08701370203367009; ...
                0.5, -1.28483063269958864, 1.87444087344698285, -2.59074567655935484, 3.65714285714285714, -5.5449639069493788, 9.73870165721154457, -24.349745171593046, 18.0];
        case 10
            D = [ -22.5, 30.4381450292819089, -12.1779467074298182, 6.94378848513395308, -4.59935476110313246, 3.2946430337491841, -2.45288417544268667, 1.82956393190324662, -1.27595483609266358, 0.5; ...
                -5.07406470297807566, 0.000000000000000222044604925031308, 7.18550286970582075, -3.35166386274677253, 2.07820799403641665, -1.44494844875145501, 1.0591544636454413, -0.783239293137908521, 0.543753738235705608, -0.212702758009188808; ...
                1.20335199285220606, -4.25929735496521467, 0, 4.36867455701018681, -2.1043501794131565, 1.33491548387825132, -0.936603213139447388, 0.676797087196086111, -0.464274958908157453, 0.18078658548924989; ...
                -0.528369376820272629, 1.52990263818160321, -3.36412586829781857, 0, 3.3873181012024447, -1.64649408398705988, 1.04618936550249364, -0.721237312721604074, 0.483462326333948145, -0.186645789393735967; ...
                0.312047255608411234, -0.845813573406425045, 1.44485031560166122, -3.02021795819934713, 0, 3.02518848775197435, -1.46805550938999385, 0.91655518033643546, -0.588082143045169481, 0.223527944742453855; ...
                -0.223527944742453855, 0.588082143045169481, -0.91655518033643546, 1.46805550938999385, -3.02518848775197435, 0, 3.02021795819934713, -1.44485031560166122, 0.845813573406425045, -0.312047255608411234; ...
                0.186645789393735967, -0.483462326333948145, 0.721237312721604074, -1.04618936550249364, 1.64649408398705988, -3.3873181012024447, 0, 3.36412586829781857, -1.52990263818160321, 0.528369376820272629; ...
                -0.18078658548924989, 0.464274958908157453, -0.676797087196086111, 0.936603213139447388, -1.33491548387825132, 2.1043501794131565, -4.36867455701018681, 0.000000000000000222044604925031308, 4.25929735496521467, -1.20335199285220606; ...
                0.212702758009188808, -0.543753738235705608, 0.783239293137908521, -1.0591544636454413, 1.44494844875145501, -2.07820799403641665, 3.35166386274677253, -7.18550286970582075, -0.000000000000000222044604925031308, 5.07406470297807566; ...
                -0.5, 1.27595483609266358, -1.82956393190324662, 2.45288417544268667, -3.2946430337491841, 4.59935476110313246, -6.94378848513395308, 12.1779467074298182, -30.4381450292819089, 22.5];
    end

end

end