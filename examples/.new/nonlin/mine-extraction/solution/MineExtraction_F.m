function I = MineExtraction_F(t,a,T,x0)
%MineExtraction_F
%    I = MineExtraction_F(t,A,T,X0)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    11-Apr-2024 11:46:18

t2 = T.*a;
I = -(a.*t2.*x0)./(t2+4.0);
end
