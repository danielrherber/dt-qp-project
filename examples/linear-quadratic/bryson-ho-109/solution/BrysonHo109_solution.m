%--------------------------------------------------------------------------
% BrysonHo109_solution.m
% Solution function for BrysonHo109 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function sol = BrysonHo109_solution(T,Yguess,in)

% extract parameters
auxdata = in.auxdata;

% symbolic variables
syms y yf t a

% saturation function
sat = @(x) (sign(x-1) + sign(x+1) - (x*(sign(x-1)-sign(x+1))))/2;

% optimal control functions
u = sat(a^2*BrysonHo109_g(t)*yf);
gu = BrysonHo109_g(t)*u;
ufun = matlabFunction(u,'Vars',{'a','t','yf'});
gufun = matlabFunction(gu,'Vars',{'a','t','yf'});

% solution function
solfun = @(y) y - ( auxdata.x0 - integral(@(t) gufun(auxdata.a,t,y),in.t0,in.tf,...
    'RelTol',1e-16,'AbsTol',1e-16));

% solve the implicit equation
yf = fzero(@(x) solfun(x), Yguess(end), optimset('Display','none','TolX',0));

% control
U = -ufun(auxdata.a,T,yf);

% state
Y(1,1) = auxdata.x0;
for i = 2:length(T)
    Y(i,1) =  Y(i-1,1) - integral(@(t) gufun(auxdata.a,t,yf), T(i-1),T(i),...
        'RelTol',1e-16,'AbsTol',1e-16);
    disp(T(i-1))
end

% objective function
Iu = integral(@(t) (ufun(auxdata.a,t,yf)).^2,in.t0,in.tf,'AbsTol',1e-16,'RelTol',1e-16);
F = auxdata.a^2/2*yf^2 + 1/2*Iu;

% save to structure
sol.T = T;
sol.U = U;
sol.Y = Y;
sol.F = F;