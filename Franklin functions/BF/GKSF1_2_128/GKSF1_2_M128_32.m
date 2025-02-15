function out1 = GKSF1_2_M128_32(T)
%GKSF1_2_M128_32
%    OUT1 = GKSF1_2_M128_32(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:41:51

t14 = ((3.75e-1 <= T) & (T < 4.0625e-1));
t4 = ((T < 9.375e-2) & (6.25e-2 <= T));
t5 = ((9.375e-2 <= T) & (T < 1.25e-1));
t27 = ((7.8125e-1 <= T) & (T < 8.125e-1));
t22 = ((6.25e-1 <= T) & (T < 6.5625e-1));
t8 = ((1.875e-1 <= T) & (T < 2.1875e-1));
t12 = ((3.125e-1 <= T) & (T < 3.4375e-1));
t3 = ((T < 6.25e-2) & (3.125e-2 <= T));
t28 = ((8.125e-1 <= T) & (T < 8.4375e-1));
t23 = ((6.5625e-1 <= T) & (T < 6.875e-1));
t18 = ((5.0e-1 <= T) & (T < 5.3125e-1));
t15 = ((4.0625e-1 <= T) & (T < 4.375e-1));
t7 = ((1.5625e-1 <= T) & (T < 1.875e-1));
t9 = ((2.1875e-1 <= T) & (T < 2.5e-1));
t10 = ((2.5e-1 <= T) & (T < 2.8125e-1));
t29 = ((8.4375e-1 <= T) & (T < 8.75e-1));
t24 = ((6.875e-1 <= T) & (T < 7.1875e-1));
t19 = ((5.3125e-1 <= T) & (T < 5.625e-1));
t13 = ((3.4375e-1 <= T) & (T < 3.75e-1));
t17 = ((4.6875e-1 <= T) & (T < 5.0e-1));
t6 = ((T < 1.5625e-1) & (1.25e-1 <= T));
t32 = ((9.375e-1 <= T) & (T <= 1.0));
t30 = ((8.75e-1 <= T) & (T < 9.0625e-1));
t25 = ((7.1875e-1 <= T) & (T < 7.5e-1));
t20 = ((5.625e-1 <= T) & (T < 5.9375e-1));
t16 = ((4.375e-1 <= T) & (T < 4.6875e-1));
t2 = ((0.0 <= T) & (T < 3.125e-2));
t11 = ((T < 3.125e-1) & (2.8125e-1 <= T));
t31 = ((9.0625e-1 <= T) & (T < 9.375e-1));
t26 = ((7.5e-1 <= T) & (T < 7.8125e-1));
t21 = ((T < 6.25e-1) & (5.9375e-1 <= T));
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t4,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*1.652293714583049e-9-2.581750123431763e-11;
elseif (t3)
    out1 = T.*-1.65213552807933e-9+7.744591259888172e-11;
elseif (t4)
    out1 = T.*1.651502782022375e-9-1.290314817824749e-10;
elseif (t5)
    out1 = T.*-1.64912998427424e-9+1.804028400578328e-10;
elseif (t6)
    out1 = T.*1.640271539302615e-9-2.30772350389274e-10;
elseif (t7)
    out1 = T.*-1.607210557116185e-9+2.766467271761635e-10;
elseif (t8)
    out1 = T.*1.483825073406544e-9-3.029224535468482e-10;
elseif (t9)
    out1 = T.*-1.023344120901132e-9+2.455208077079558e-10;
elseif (t10)
    out1 = T.*-6.951942054049244e-10+1.63483328833904e-10;
elseif (t11)
    out1 = T.*7.108866558193682e-9-2.031408760928204e-9;
elseif (t12)
    out1 = T.*-3.104501764306836e-8+9.891680051966184e-9;
elseif (t13)
    out1 = T.*1.203759496298348e-7-4.215927744809427e-8;
elseif (t14)
    out1 = T.*-4.537635264964384e-7+1.731430260992582e-7;
elseif (t15)
    out1 = T.*1.697982901969171e-6-7.010039604648955e-7;
elseif (t16)
    out1 = T.*-6.341472826993617e-6+2.816257920956324e-6;
elseif (t17)
    out1 = T.*2.367121315154795e-5-1.125218863148504e-5;
elseif (t18)
    out1 = T.*-8.834668452451811e-5+4.475676020654799e-5;
elseif (t19)
    out1 = T.*3.297188296916723e-4-1.773405442208032e-4;
elseif (t20)
    out1 = T.*-1.230531938999518e-3+7.003005131679913e-4;
elseif (t21)
    out1 = T.*4.592412231045839e-3-2.757072587796439e-3;
elseif (t22)
    out1 = T.*-1.713912028992372e-2+1.082513523780953e-2;
elseif (t23)
    out1 = T.*6.396407223319737e-2-4.239883485548868e-2;
elseif (t24)
    out1 = T.*-2.387171719468118e-1+1.656945205182676e-1;
elseif (t25)
    out1 = T.*8.909046188575304e-1-6.462211416223533e-1;
elseif (t26)
    out1 = T.*-3.324901306794967+2.51563330261702;
elseif (t27)
    out1 = T.*1.24087006116219e+1-9.776243196146156;
elseif (t28)
    out1 = T.*-4.630990114316509e+1+3.793262072961827e+1;
elseif (t29)
    out1 = T.*1.728309039637766e+2-1.469674335793638e+2;
elseif (t30)
    out1 = T.*-3.455509778848626e+2+3.066167130381955e+2;
elseif (t31)
    out1 = T.*3.109847970920917e+2-2.883688330346693e+2;
elseif (t32)
    out1 = T.*-7.630594174001389e+1+7.471623462042968e+1;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
