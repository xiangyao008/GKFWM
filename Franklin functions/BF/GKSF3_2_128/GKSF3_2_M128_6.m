function out1 = GKSF3_2_M128_6(T)
%GKSF3_2_M128_6
%    OUT1 = GKSF3_2_M128_6(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 10:03:28

t2 = ((0.0 <= T) & (T < 2.5e-1));
t3 = ((T < 5.0e-1) & (2.5e-1 <= T));
t4 = ((5.0e-1 <= T) & (T <= 1.0));
t5 = T.^2;
t6 = T.^3;
if ~all(cellfun(@isscalar,{T,t2,t3,t4}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*-1.262337757315563e+2+t5.*7.747932419901005e+2-t6.*1.291605832485045e+3+4.27670627478365;
elseif (t3)
    out1 = T.*1.85094226072574e+2-t5.*4.705187652264206e+2+t6.*3.68810177136983e+2-2.166729387556054e+1;
elseif (t4)
    out1 = T.*-1.92390976114396e+2+t5.*2.844516391475194e+2-t6.*1.345034257789769e+2+4.124690648893445e+1;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
