%--------------------------------------------------------------------------
% BrysonHo109_solution.m
% Solution function for BrysonHo109 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function sol = BrysonHo109_solution(T,Yguess,p)

% symbolic variables
syms y yf t a

% saturation function
sat = @(x) (sign(x-1) + sign(x+1) - (x*(sign(x-1)-sign(x+1))))/2;

% optimal control functions
u = sat(a^2*p.g*yf);
gu = p.g*u;
ufun = matlabFunction(u,'Vars',{'a','t','yf'});
gufun = matlabFunction(gu,'Vars',{'a','t','yf'});

% solution function
solfun = @(y) y - ( p.x0 - integral(@(t) gufun(p.a,t,y),p.t0,p.tf,...
    'RelTol',1e-16,'AbsTol',1e-16));

% solve the implicit equation
yf = fzero(@(x) solfun(x), Yguess(end), optimset('Display','none','TolX',0));

% control
U = -ufun(p.a,T,yf);

% state
Y(1,1) = p.x0;
for i = 2:length(T)
    Y(i,1) =  Y(i-1,1) - integral(@(t) gufun(p.a,t,yf), T(i-1),T(i),...
        'RelTol',1e-16,'AbsTol',1e-16);
end

% objective function
Iu = integral(@(t) (ufun(p.a,t,yf)).^2,p.t0,p.tf,'AbsTol',1e-16,'RelTol',1e-16);
F = p.a^2/2*yf^2 + 1/2*Iu;

% save to structure
sol.T = T;
sol.U = U;
sol.Y = Y;
sol.F = F;