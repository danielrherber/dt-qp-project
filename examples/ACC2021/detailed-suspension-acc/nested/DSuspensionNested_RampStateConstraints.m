%--------------------------------------------------------------------------
% DSuspensionNested_RampStateConstraints.m
% State-dependent constraints for ramp load case
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function setup = DSuspensionNested_RampStateConstraints(xd,param,setup)

% extract plant variables
d = xd(1);  % spring wire diameter
D = xd(2);  % spring helix diameter
p = xd(3);  % spring pitch
Na = xd(4); % number of active spring coils

% extract problem parameters
ms = param.ms; g = param.g; G = param.G; Q = param.Q; LB = param.LB;

% rough state constraints
setup = DSuspensionNested_RoughStateConstraints(xd,param,setup);

% remove unnecessary constraints
setup.UB(end) = [];
setup.LB(end) = [];

%--------------------------------------------------------------------------
% rattlespace constraint (gi1)
% In DTCD: g(1) = - L0 + Ls + LB + delta_g + x(3);
% NOTE: make UB/LB
b1 = -(Na*p+2*d)+(d*(Na+Q-1))+0.02+(4*ms*g*D*Na*(2*D^2+d^2)/(G*d^4)); % same value as above
% Z(end+1).linear(1).right = 2;
% Z(end).linear(1).matrix = [0,0,1,0];
% Z(end).b = -b1;

setup.UB(end+1).right = 2;
setup.UB(end).matrix = [inf;inf;-b1;inf];

% validation
% x3 = 2*(rand(1000,1)-0.5);
% C = D/d;
% L0 = p*Na + 2*d;
% Ls = d*(Na + 1.75 - 1);
% ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
% delta_g = param.ms*9.81/ks;
% gDTCD = - L0 + Ls + LB + delta_g + x3;
% gDTQP = x3+b1;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)

%--------------------------------------------------------------------------
% control force bounds
% Z(end+1).linear(1).right = 1;
% Z(end).linear(1).matrix = 1;
% Z(end).b = ;

end