function out1 = GKSF1_2_M128_24(T)
%GKSF1_2_M128_24
%    OUT1 = GKSF1_2_M128_24(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:41:49

t14 = ((3.75e-1 <= T) & (T < 4.0625e-1));
t21 = ((7.5e-1 <= T) & (T < 8.125e-1));
t4 = ((T < 9.375e-2) & (6.25e-2 <= T));
t5 = ((9.375e-2 <= T) & (T < 1.25e-1));
t8 = ((1.875e-1 <= T) & (T < 2.1875e-1));
t12 = ((3.125e-1 <= T) & (T < 3.4375e-1));
t19 = ((6.25e-1 <= T) & (T < 6.875e-1));
t3 = ((T < 6.25e-2) & (3.125e-2 <= T));
t15 = ((4.0625e-1 <= T) & (T < 4.375e-1));
t22 = ((8.125e-1 <= T) & (T < 8.75e-1));
t7 = ((1.5625e-1 <= T) & (T < 1.875e-1));
t9 = ((2.1875e-1 <= T) & (T < 2.5e-1));
t16 = ((4.375e-1 <= T) & (T < 5.0e-1));
t10 = ((2.5e-1 <= T) & (T < 2.8125e-1));
t17 = ((5.0e-1 <= T) & (T < 5.625e-1));
t13 = ((3.4375e-1 <= T) & (T < 3.75e-1));
t20 = ((6.875e-1 <= T) & (T < 7.5e-1));
t6 = ((T < 1.5625e-1) & (1.25e-1 <= T));
t24 = ((9.375e-1 <= T) & (T <= 1.0));
t23 = ((8.75e-1 <= T) & (T < 9.375e-1));
t2 = ((0.0 <= T) & (T < 3.125e-2));
t11 = ((T < 3.125e-1) & (2.8125e-1 <= T));
t18 = ((T < 6.25e-1) & (5.625e-1 <= T));
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t3,t4,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*-1.127369275945655e-4+1.174293106342717e-6;
elseif (t3)
    out1 = T.*3.382299403230738e-4-1.291842151608351e-5;
elseif (t4)
    out1 = T.*-1.240201991235877e-3+8.573357420635089e-5;
elseif (t5)
    out1 = T.*4.622597182144634e-3-4.63903848298072e-4;
elseif (t6)
    out1 = T.*-1.725020589482338e-2+2.27019653632293e-3;
elseif (t7)
    out1 = T.*6.437824555459599e-2-1.048424900264885e-2;
elseif (t8)
    out1 = T.*-2.402627954833895e-1+4.663594619197344e-2;
elseif (t9)
    out1 = T.*8.966729555352947e-1-2.020687493433637e-1;
elseif (t10)
    out1 = T.*-3.346429045814208+8.587067509940119e-1;
elseif (t11)
    out1 = T.*1.248904324684056e+1-3.595019831315141;
elseif (t12)
    out1 = T.*-4.660974396054944e+1+1.487335117099423e+1;
elseif (t13)
    out1 = T.*1.739499326142677e+2-6.094403765159918e+1;
elseif (t14)
    out1 = T.*-3.463387348333829e+2+1.341642126412698e+2;
elseif (t15)
    out1 = T.*3.028512516873363e+2-1.295692193827724e+2;
elseif (t16)
    out1 = T.*-5.940500674914166e+1+2.891789368318675e+1;
elseif (t17)
    out1 = T.*1.591752362372978e+1-8.743371503248968;
elseif (t18)
    out1 = T.*-4.265087745777463+2.609347392098857;
elseif (t19)
    out1 = T.*1.142827359380071-7.70599548624602e-1;
elseif (t20)
    out1 = T.*-3.062216917428209e-1+2.256216740223862e-1;
elseif (t21)
    out1 = T.*8.205940759121259e-2-6.558915047813895e-2;
elseif (t22)
    out1 = T.*-2.201593862203036e-2+1.897206832012094e-2;
elseif (t23)
    out1 = T.*6.004346896913854e-3-5.545681508955248e-3;
elseif (t24)
    out1 = T.*-2.001448965634614e-3+1.95975211218394e-3;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
