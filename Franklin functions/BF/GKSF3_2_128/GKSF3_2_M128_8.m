function out1 = GKSF3_2_M128_8(T)
%GKSF3_2_M128_8
%    OUT1 = GKSF3_2_M128_8(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 10:03:29

t2 = ((0.0 <= T) & (T < 1.25e-1));
t3 = ((T < 2.5e-1) & (1.25e-1 <= T));
t4 = ((T < 5.0e-1) & (2.5e-1 <= T));
t5 = ((T < 7.5e-1) & (5.0e-1 <= T));
t6 = ((7.5e-1 <= T) & (T <= 1.0));
t7 = T.^2;
t8 = T.^3;
if ~all(cellfun(@isscalar,{T,t2,t3,t4,t5,t6}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*-3.499098503524895e+2+t7.*4.116718897739717e+3-t8.*1.328298180248964e+4+6.216122808901451;
elseif (t3)
    out1 = T.*4.064515806696038e+2-t7.*1.93417255043703e+3+t8.*2.852728725981686e+3-2.529893681701911e+1;
elseif (t4)
    out1 = T.*-2.000040965688617e+2+t7.*4.916501585168319e+2-t8.*3.81701552623464e+2+2.523903628618636e+1;
elseif (t5)
    out1 = T.*2.274427534853492e+2-t7.*3.6324354159159e+2+t8.*1.882275807821506e+2-4.600210538951547e+1;
elseif (t6)
    out1 = T.*-5.092918125119663e+2+t7.*6.190692130714973e+2-t8.*2.483558657347771e+2+1.381815361098134e+2;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
