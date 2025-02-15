function out1 = GKSF1_2_M128_22(T)
%GKSF1_2_M128_22
%    OUT1 = GKSF1_2_M128_22(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:41:49

t19 = ((7.5e-1 <= T) & (T < 8.125e-1));
t4 = ((T < 9.375e-2) & (6.25e-2 <= T));
t5 = ((9.375e-2 <= T) & (T < 1.25e-1));
t8 = ((1.875e-1 <= T) & (T < 2.1875e-1));
t13 = ((3.75e-1 <= T) & (T < 4.375e-1));
t17 = ((6.25e-1 <= T) & (T < 6.875e-1));
t3 = ((T < 6.25e-2) & (3.125e-2 <= T));
t20 = ((8.125e-1 <= T) & (T < 8.75e-1));
t7 = ((1.5625e-1 <= T) & (T < 1.875e-1));
t9 = ((2.1875e-1 <= T) & (T < 2.5e-1));
t12 = ((3.125e-1 <= T) & (T < 3.75e-1));
t14 = ((4.375e-1 <= T) & (T < 5.0e-1));
t10 = ((2.5e-1 <= T) & (T < 2.8125e-1));
t15 = ((5.0e-1 <= T) & (T < 5.625e-1));
t18 = ((6.875e-1 <= T) & (T < 7.5e-1));
t6 = ((T < 1.5625e-1) & (1.25e-1 <= T));
t22 = ((9.375e-1 <= T) & (T <= 1.0));
t21 = ((8.75e-1 <= T) & (T < 9.375e-1));
t2 = ((0.0 <= T) & (T < 3.125e-2));
t11 = ((T < 3.125e-1) & (2.8125e-1 <= T));
t16 = ((T < 6.25e-1) & (5.625e-1 <= T));
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t3,t4,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*-2.187223310654233e-2+2.278357201571527e-4;
elseif (t3)
    out1 = T.*6.561671520481313e-2-2.506193914572705e-3;
elseif (t4)
    out1 = T.*-2.405946435976576e-1+1.663201601058171e-2;
elseif (t5)
    out1 = T.*8.967618750680592e-1-8.999515761432924e-2;
elseif (t6)
    out1 = T.*-3.346452872548387+4.404066858377265e-1;
elseif (t7)
    out1 = T.*1.248904963099277e+1-2.033890580340579;
elseif (t8)
    out1 = T.*-4.660974566775205e+1+9.047133538174074;
elseif (t9)
    out1 = T.*1.739499330556665e+2-3.920029618257373e+1;
elseif (t10)
    out1 = T.*-3.463387349182568e+2+9.08718708109071e+1;
elseif (t11)
    out1 = T.*3.028512516776132e+2-9.171281291918136e+1;
elseif (t12)
    out1 = T.*-5.94050067335622e+1+2.149226783431097e+1;
elseif (t13)
    out1 = T.*1.591752358082408e+1-6.753681033583888;
elseif (t14)
    out1 = T.*-4.265087589734132+2.076211353535331;
elseif (t15)
    out1 = T.*1.142826778112448-6.27745830387959e-1;
elseif (t16)
    out1 = T.*-3.062195227156577e-1+1.873427138278502e-1;
elseif (t17)
    out1 = T.*8.205131275018327e-2-5.532655833830038e-2;
elseif (t18)
    out1 = T.*-2.198572828507536e-2+1.619890737343992e-2;
elseif (t19)
    out1 = T.*5.891600390118172e-3-4.709089132955228e-3;
elseif (t20)
    out1 = T.*-1.580673275397394e-3+1.362133220276169e-3;
elseif (t21)
    out1 = T.*4.310927114717641e-4-3.981620182343443e-4;
elseif (t22)
    out1 = T.*-1.436975704903484e-4+1.407038711051363e-4;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
