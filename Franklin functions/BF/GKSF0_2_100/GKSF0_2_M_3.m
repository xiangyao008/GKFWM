function out1 = GKSF0_2_M_3(T)
%GKSF0_2_M_3
%    OUT1 = GKSF0_2_M_3(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-04 21:29:50

t2 = ((T < 1.0./4.0) & (0.0 <= T));
t3 = ((T < 1.0./2.0) & (1.0./4.0 <= T));
t4 = ((T <= 1.0) & (1.0./2.0 <= T));
if ~all(cellfun(@isscalar,{T,t2,t3,t4}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = -1.414213562373095;
elseif (t3)
    out1 = 1.414213562373095;
elseif (t4)
    out1 = 0.0;
elseif (T == 1.0)
    out1 = 4.242640687119285;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
