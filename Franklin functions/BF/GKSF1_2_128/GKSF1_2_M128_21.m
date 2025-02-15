function out1 = GKSF1_2_M128_21(T)
%GKSF1_2_M128_21
%    OUT1 = GKSF1_2_M128_21(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:41:49

t18 = ((7.5e-1 <= T) & (T < 8.125e-1));
t4 = ((T < 9.375e-2) & (6.25e-2 <= T));
t5 = ((9.375e-2 <= T) & (T < 1.25e-1));
t8 = ((1.875e-1 <= T) & (T < 2.1875e-1));
t12 = ((3.75e-1 <= T) & (T < 4.375e-1));
t16 = ((6.25e-1 <= T) & (T < 6.875e-1));
t3 = ((T < 6.25e-2) & (3.125e-2 <= T));
t19 = ((8.125e-1 <= T) & (T < 8.75e-1));
t7 = ((1.5625e-1 <= T) & (T < 1.875e-1));
t9 = ((2.1875e-1 <= T) & (T < 2.5e-1));
t11 = ((3.125e-1 <= T) & (T < 3.75e-1));
t13 = ((4.375e-1 <= T) & (T < 5.0e-1));
t14 = ((5.0e-1 <= T) & (T < 5.625e-1));
t17 = ((6.875e-1 <= T) & (T < 7.5e-1));
t6 = ((T < 1.5625e-1) & (1.25e-1 <= T));
t10 = ((T < 3.125e-1) & (2.5e-1 <= T));
t21 = ((9.375e-1 <= T) & (T <= 1.0));
t20 = ((8.75e-1 <= T) & (T < 9.375e-1));
t2 = ((0.0 <= T) & (T < 3.125e-2));
t15 = ((T < 6.25e-1) & (5.625e-1 <= T));
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t3,t4,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*-3.046410019002222e-1+3.173343719875313e-3;
elseif (t3)
    out1 = T.*9.13923024860569e-1-3.490678211639941e-2;
elseif (t4)
    out1 = T.*-3.351051116698631+2.316541017310506e-1;
elseif (t5)
    out1 = T.*1.249028146105285e+1-1.253470827433151;
elseif (t6)
    out1 = T.*-4.66100747465142e+1+6.13407369851273;
elseif (t7)
    out1 = T.*1.739500175439145e+2-2.832844072186675e+1;
elseif (t8)
    out1 = T.*-3.46338750641644e+2+6.922570331292546e+1;
elseif (t9)
    out1 = T.*3.028512506176508e+2-7.278460946254527e+1;
elseif (t10)
    out1 = T.*-5.94050059584437e+1+1.777945468147835e+1;
elseif (t11)
    out1 = T.*1.591752337294543e+1-5.758835734580751;
elseif (t12)
    out1 = T.*-4.265087533338022+1.809643355275544;
elseif (t13)
    out1 = T.*1.142826760406657-5.56319148237753e-1;
elseif (t14)
    out1 = T.*-3.06219508288607e-1+1.682039861098792e-1;
elseif (t15)
    out1 = T.*8.205127274777086e-2-5.019832822308338e-2;
elseif (t16)
    out1 = T.*-2.198558270247638e-2+1.482470643332114e-2;
elseif (t17)
    out1 = T.*5.891058062134646e-3-4.340484092348934e-3;
elseif (t18)
    out1 = T.*-1.578649546062209e-3+1.261796613798707e-3;
elseif (t19)
    out1 = T.*4.235401221142072e-4-3.649824915946308e-4;
elseif (t20)
    out1 = T.*-1.155109423947161e-4+1.066871898506771e-4;
elseif (t21)
    out1 = T.*3.850364746484118e-5-3.77014881426579e-5;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
