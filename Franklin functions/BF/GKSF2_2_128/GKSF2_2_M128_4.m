function out1 = GKSF2_2_M128_4(T)
%GKSF2_2_M128_4
%    OUT1 = GKSF2_2_M128_4(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:45:34

t2 = ((0.0 <= T) & (T < 5.0e-1));
t3 = ((5.0e-1 <= T) & (T <= 1.0));
t4 = T.^2;
if ~all(cellfun(@isscalar,{T,t2,t3}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*2.23606797749979e+1-t4.*3.577708763999664e+1-2.23606797749979;
elseif (t3)
    out1 = T.*-4.919349550499537e+1+t4.*3.577708763999664e+1+1.565247584249853e+1;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
