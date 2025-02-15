function out1 = GKSF0_2_M_12(T)
%GKSF0_2_M_12
%    OUT1 = GKSF0_2_M_12(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-04 21:29:51

t2 = ((T < 7.0./1.6e+1) & (3.0./8.0 <= T));
t3 = ((T < 1.0./2.0) & (7.0./1.6e+1 <= T));
if ~all(cellfun(@isscalar,{T,t2,t3}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 3.0./8.0)
    out1 = 0.0;
elseif (t2)
    out1 = -2.82842712474619;
elseif (t3)
    out1 = 2.82842712474619;
elseif (1.0./2.0 <= T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
