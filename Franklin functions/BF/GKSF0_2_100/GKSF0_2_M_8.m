function out1 = GKSF0_2_M_8(T)
%GKSF0_2_M_8
%    OUT1 = GKSF0_2_M_8(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-04 21:29:51

t2 = ((T < 7.0./8.0) & (3.0./4.0 <= T));
t3 = ((T <= 1.0) & (7.0./8.0 <= T));
if ~all(cellfun(@isscalar,{T,t2,t3}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 3.0./4.0)
    out1 = 0.0;
elseif (t2)
    out1 = -2.0;
elseif (t3)
    out1 = 2.0;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
