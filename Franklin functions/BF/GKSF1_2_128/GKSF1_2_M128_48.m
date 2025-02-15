function out1 = GKSF1_2_M128_48(T)
%GKSF1_2_M128_48
%    OUT1 = GKSF1_2_M128_48(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:41:56

t14 = ((1.875e-1 <= T) & (T < 2.03125e-1));
t4 = ((T < 4.6875e-2) & (3.125e-2 <= T));
t5 = ((4.6875e-2 <= T) & (T < 6.25e-2));
t47 = ((9.375e-1 <= T) & (T < 9.6875e-1));
t27 = ((3.90625e-1 <= T) & (T < 4.0625e-1));
t42 = ((7.8125e-1 <= T) & (T < 8.125e-1));
t22 = ((3.125e-1 <= T) & (T < 3.28125e-1));
t37 = ((6.25e-1 <= T) & (T < 6.5625e-1));
t8 = ((9.375e-2 <= T) & (T < 1.09375e-1));
t12 = ((1.5625e-1 <= T) & (T < 1.71875e-1));
t3 = ((T < 3.125e-2) & (1.5625e-2 <= T));
t28 = ((4.0625e-1 <= T) & (T < 4.21875e-1));
t43 = ((8.125e-1 <= T) & (T < 8.4375e-1));
t23 = ((3.28125e-1 <= T) & (T < 3.4375e-1));
t38 = ((6.5625e-1 <= T) & (T < 6.875e-1));
t18 = ((2.5e-1 <= T) & (T < 2.65625e-1));
t33 = ((5.0e-1 <= T) & (T < 5.3125e-1));
t15 = ((2.03125e-1 <= T) & (T < 2.1875e-1));
t7 = ((7.8125e-2 <= T) & (T < 9.375e-2));
t9 = ((1.09375e-1 <= T) & (T < 1.25e-1));
t10 = ((1.25e-1 <= T) & (T < 1.40625e-1));
t29 = ((4.21875e-1 <= T) & (T < 4.375e-1));
t44 = ((8.4375e-1 <= T) & (T < 8.75e-1));
t24 = ((3.4375e-1 <= T) & (T < 3.59375e-1));
t39 = ((6.875e-1 <= T) & (T < 7.1875e-1));
t19 = ((2.65625e-1 <= T) & (T < 2.8125e-1));
t34 = ((5.3125e-1 <= T) & (T < 5.625e-1));
t13 = ((1.71875e-1 <= T) & (T < 1.875e-1));
t17 = ((2.34375e-1 <= T) & (T < 2.5e-1));
t32 = ((4.6875e-1 <= T) & (T < 5.0e-1));
t6 = ((T < 7.8125e-2) & (6.25e-2 <= T));
t30 = ((4.375e-1 <= T) & (T < 4.53125e-1));
t45 = ((8.75e-1 <= T) & (T < 9.0625e-1));
t25 = ((3.59375e-1 <= T) & (T < 3.75e-1));
t40 = ((7.1875e-1 <= T) & (T < 7.5e-1));
t20 = ((2.8125e-1 <= T) & (T < 2.96875e-1));
t48 = ((9.6875e-1 <= T) & (T <= 1.0));
t35 = ((5.625e-1 <= T) & (T < 5.9375e-1));
t16 = ((2.1875e-1 <= T) & (T < 2.34375e-1));
t2 = ((0.0 <= T) & (T < 1.5625e-2));
t11 = ((T < 1.5625e-1) & (1.40625e-1 <= T));
t31 = ((4.53125e-1 <= T) & (T < 4.6875e-1));
t46 = ((9.0625e-1 <= T) & (T < 9.375e-1));
t26 = ((3.75e-1 <= T) & (T < 3.90625e-1));
t41 = ((7.5e-1 <= T) & (T < 7.8125e-1));
t21 = ((T < 3.125e-1) & (2.96875e-1 <= T));
t36 = ((T < 6.25e-1) & (5.9375e-1 <= T));
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*2.984911400059972e-8-1.554641354077162e-10;
elseif (t3)
    out1 = T.*-4.97500276053029e-8+1.088272452184512e-9;
elseif (t4)
    out1 = T.*4.975905321090448e-8-2.021386323321969e-9;
elseif (t5)
    out1 = T.*-2.989424202873573e-8+1.712361891036167e-9;
elseif (t6)
    out1 = T.*3.002059803180842e-8-2.032315612747843e-9;
elseif (t7)
    out1 = T.*-5.039082826180807e-8+4.249827066440946e-9;
elseif (t8)
    out1 = T.*5.214404486431055e-8-5.362817289132675e-9;
elseif (t9)
    out1 = T.*-3.877324703180557e-8+4.581261512005026e-9;
elseif (t10)
    out1 = T.*4.497542492592719e-8-5.887322482711569e-9;
elseif (t11)
    out1 = T.*-4.68562823396755e-8+7.026511351513809e-9;
elseif (t12)
    out1 = T.*-3.14104445166445e-8+4.613099191665215e-9;
elseif (t13)
    out1 = T.*3.100462569041898e-7-5.407477136504068e-8;
elseif (t14)
    out1 = T.*-1.248759285958644e-6+2.382012679217406e-7;
elseif (t15)
    out1 = T.*4.725349812936164e-6-9.752896427912672e-7;
elseif (t16)
    out1 = T.*-1.78794735761749e-5+3.969515473576778e-6;
elseif (t17)
    out1 = T.*6.723276398813317e-5-1.597866520555793e-5;
elseif (t18)
    out1 = T.*-2.514120200761659e-4+6.368253081051685e-5;
elseif (t19)
    out1 = T.*9.385619937554001e-4-2.524043166134929e-4;
elseif (t20)
    out1 = T.*-3.502955346744807e-3+9.967724354021904e-4;
elseif (t21)
    out1 = T.*1.307337878493389e-2-3.924326759939924e-3;
elseif (t22)
    out1 = T.*-4.879059959162256e-2+1.540816648273397e-2;
elseif (t23)
    out1 = T.*1.820890593762807e-1-6.034922161610928e-2;
elseif (t24)
    out1 = T.*-6.795657573354289e-1+2.358446216285409e-1;
elseif (t25)
    out1 = T.*2.536174089283175-9.1981188575002e-1;
elseif (t26)
    out1 = T.*-9.465130657494809+3.580677394291724;
elseif (t27)
    out1 = T.*3.532434863692032e+1-1.391521295508919e+1;
elseif (t28)
    out1 = T.*-1.318322640657764e+2+5.399216095538136e+1;
elseif (t29)
    out1 = T.*4.920047077478581e+2-2.091890615284957e+2;
elseif (t30)
    out1 = T.*-9.795938719472011e+2+4.346353170880927e+2;
elseif (t31)
    out1 = T.*8.565926950117246e+2-3.973867210651705e+2;
elseif (t32)
    out1 = T.*-1.680227323985024e+2+8.290176053337339e+1;
elseif (t33)
    out1 = T.*4.502155545382426e+1-2.362038339278993e+1;
elseif (t34)
    out1 = T.*-1.206348942064578e+1+6.706046696772282;
elseif (t35)
    out1 = T.*3.232402237998901-1.897892361215353;
elseif (t36)
    out1 = T.*-8.661195587929327e-1+5.356049556297984e-1;
elseif (t37)
    out1 = T.*2.320760291185219e-1-1.507672868148607e-1;
elseif (t38)
    out1 = T.*-6.21845804874943e-2+4.234123823908739e-2;
elseif (t39)
    out1 = T.*1.666230651920194e-2-1.186599657801627e-2;
elseif (t40)
    out1 = T.*-4.464650151547214e-3+3.319003529084679e-3;
elseif (t41)
    out1 = T.*1.196296368263839e-3-9.267063607736106e-4;
elseif (t42)
    out1 = T.*-3.205421653108893e-4+2.583237435816457e-4;
elseif (t43)
    out1 = T.*8.588369930701711e-5-7.189727142040322e-5;
elseif (t44)
    out1 = T.*-2.300860077286867e-5+1.998060677200041e-5;
elseif (t45)
    out1 = T.*6.164391375057538e-6-5.545761357435028e-6;
elseif (t46)
    out1 = T.*-1.653527257409357e-6+1.539227403238096e-6;
elseif (t47)
    out1 = T.*4.502879700497299e-7-4.33099372504798e-7;
elseif (t48)
    out1 = T.*-1.493355693172644e-7+1.477859312569778e-7;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
