function F = TurnerChunJuang1_F(S,a,em,eta,f,t0,um,y0)

% constants
Q = em^(-2);
R = um^(-2);

% suppress warning
warning('off','MATLAB:integral:MaxIntervalCountReached')

% calculate the integral term
I = integral(@(t) Q.*(TurnerChunJuang1_Y(S,a,em,eta,f,t,t0,um,y0)-eta).^2+...
    R.*TurnerChunJuang1_U(S,a,em,eta,f,t,t0,um,y0).^2,t0,f,...
    'RelTol',0,'AbsTol',0);

% unsuppress warning
warning('on','MATLAB:integral:MaxIntervalCountReached')

% final state
yf = TurnerChunJuang1_Y(S,a,em,eta,f,f,t0,um,y0);

% calculate final value
F = S.*(yf-eta).^2 + I;