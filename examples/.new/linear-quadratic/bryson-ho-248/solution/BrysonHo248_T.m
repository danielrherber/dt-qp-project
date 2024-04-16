function [t1,t2] = BrysonHo248_T(a,b,c,tf)
syms t1 t2

% find transition times for initially positive control
eqs_pos = [ c*t2 - c*tf + (c*t2^2)/2 + (c*tf^2)/2 + a*exp(t1 - t2) + b*t1*exp(t1 - t2) - b*t2*exp(t2 - t1) + b*tf*exp(t2 - t1) + c*t1*exp(t1 - t2) + a*t2*sinh(t1 - t2) - a*tf*sinh(t1 - t2) - (c*t1^2*exp(t1 - t2))/2 - c*t2*tf - (c*t1^2*t2*sinh(t1 - t2))/2 + (c*t1^2*tf*sinh(t1 - t2))/2 + c*t1*t2*cosh(t1 - t2) - c*t1*tf*cosh(t1 - t2) + b*t1*t2*sinh(t1 - t2) - b*t1*tf*sinh(t1 - t2) == 0, c*tf - c*t2 + b*(cosh(t1 - t2) - sinh(t1 - t2)) - a*sinh(t1 - t2) - c*t1*cosh(t1 - t2) - b*t1*sinh(t1 - t2) + (c*t1^2*sinh(t1 - t2))/2 == 0];
sol_pos = vpasolve(eqs_pos,[t1,t2],[tf/3,tf*2/3]);

% convert to double
t1_pos = double(sol_pos.t1);
t2_pos = double(sol_pos.t2);

% find transition times for initially negative control
c = -c;
eqs_neg = [ c*t2 - c*tf + (c*t2^2)/2 + (c*tf^2)/2 + a*exp(t1 - t2) + b*t1*exp(t1 - t2) - b*t2*exp(t2 - t1) + b*tf*exp(t2 - t1) + c*t1*exp(t1 - t2) + a*t2*sinh(t1 - t2) - a*tf*sinh(t1 - t2) - (c*t1^2*exp(t1 - t2))/2 - c*t2*tf - (c*t1^2*t2*sinh(t1 - t2))/2 + (c*t1^2*tf*sinh(t1 - t2))/2 + c*t1*t2*cosh(t1 - t2) - c*t1*tf*cosh(t1 - t2) + b*t1*t2*sinh(t1 - t2) - b*t1*tf*sinh(t1 - t2) == 0, c*tf - c*t2 + b*(cosh(t1 - t2) - sinh(t1 - t2)) - a*sinh(t1 - t2) - c*t1*cosh(t1 - t2) - b*t1*sinh(t1 - t2) + (c*t1^2*sinh(t1 - t2))/2 == 0];
sol_neg = vpasolve(eqs_neg,[t1,t2],[tf/3,tf*2/3]);

% convert to double
t1_neg = double(sol_neg.t1);
t2_neg = double(sol_neg.t2);

% determine if the positive or negative solution is correct
try
    if (t1_pos >= 0) && (t2_pos <= tf)
        t1 = t1_pos;
        t2 = t2_pos;
    elseif (t1_neg >= 0) && (t2_neg <= tf)
        t1 = t1_neg;
        t2 = t2_neg;
    else
        error('could not determine the transition times or infeasible problem')
    end
catch
    error('could not determine the transition times or infeasible problem')
end

end