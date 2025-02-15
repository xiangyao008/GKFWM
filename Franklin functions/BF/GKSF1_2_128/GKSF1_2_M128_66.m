function out1 = GKSF1_2_M128_66(T)
%GKSF1_2_M128_66
%    OUT1 = GKSF1_2_M128_66(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:42:05

t52 = ((7.65625e-1 <= T) & (T < 7.8125e-1));
t15 = ((1.875e-1 <= T) & (T < 2.03125e-1));
t42 = ((6.09375e-1 <= T) & (T < 6.25e-1));
t63 = ((9.375e-1 <= T) & (T < 9.53125e-1));
t5 = ((T < 4.6875e-2) & (3.125e-2 <= T));
t6 = ((4.6875e-2 <= T) & (T < 6.25e-2));
t53 = ((7.8125e-1 <= T) & (T < 7.96875e-1));
t33 = ((4.6875e-1 <= T) & (T < 4.84375e-1));
t43 = ((6.25e-1 <= T) & (T < 6.40625e-1));
t28 = ((3.90625e-1 <= T) & (T < 4.0625e-1));
t64 = ((9.53125e-1 <= T) & (T < 9.6875e-1));
t23 = ((3.125e-1 <= T) & (T < 3.28125e-1));
t54 = ((7.96875e-1 <= T) & (T < 8.125e-1));
t9 = ((9.375e-2 <= T) & (T < 1.09375e-1));
t44 = ((6.40625e-1 <= T) & (T < 6.5625e-1));
t65 = ((9.6875e-1 <= T) & (T < 9.84375e-1));
t13 = ((1.5625e-1 <= T) & (T < 1.71875e-1));
t3 = ((T < 1.5625e-2) & (7.8125e-3 <= T));
t4 = ((T < 3.125e-2) & (1.5625e-2 <= T));
t55 = ((8.125e-1 <= T) & (T < 8.28125e-1));
t45 = ((6.5625e-1 <= T) & (T < 6.71875e-1));
t29 = ((4.0625e-1 <= T) & (T < 4.21875e-1));
t35 = ((5.0e-1 <= T) & (T < 5.15625e-1));
t24 = ((3.28125e-1 <= T) & (T < 3.4375e-1));
t56 = ((8.28125e-1 <= T) & (T < 8.4375e-1));
t19 = ((2.5e-1 <= T) & (T < 2.65625e-1));
t16 = ((2.03125e-1 <= T) & (T < 2.1875e-1));
t46 = ((6.71875e-1 <= T) & (T < 6.875e-1));
t8 = ((7.8125e-2 <= T) & (T < 9.375e-2));
t10 = ((1.09375e-1 <= T) & (T < 1.25e-1));
t36 = ((5.15625e-1 <= T) & (T < 5.3125e-1));
t57 = ((8.4375e-1 <= T) & (T < 8.59375e-1));
t11 = ((1.25e-1 <= T) & (T < 1.40625e-1));
t47 = ((6.875e-1 <= T) & (T < 7.03125e-1));
t30 = ((4.21875e-1 <= T) & (T < 4.375e-1));
t37 = ((5.3125e-1 <= T) & (T < 5.46875e-1));
t25 = ((3.4375e-1 <= T) & (T < 3.59375e-1));
t58 = ((8.59375e-1 <= T) & (T < 8.75e-1));
t20 = ((2.65625e-1 <= T) & (T < 2.8125e-1));
t48 = ((7.03125e-1 <= T) & (T < 7.1875e-1));
t14 = ((1.71875e-1 <= T) & (T < 1.875e-1));
t38 = ((5.46875e-1 <= T) & (T < 5.625e-1));
t18 = ((2.34375e-1 <= T) & (T < 2.5e-1));
t59 = ((8.75e-1 <= T) & (T < 8.90625e-1));
t7 = ((T < 7.8125e-2) & (6.25e-2 <= T));
t49 = ((7.1875e-1 <= T) & (T < 7.34375e-1));
t31 = ((4.375e-1 <= T) & (T < 4.53125e-1));
t39 = ((5.625e-1 <= T) & (T < 5.78125e-1));
t26 = ((3.59375e-1 <= T) & (T < 3.75e-1));
t34 = ((4.84375e-1 <= T) & (T < 5.0e-1));
t60 = ((8.90625e-1 <= T) & (T < 9.0625e-1));
t21 = ((2.8125e-1 <= T) & (T < 2.96875e-1));
t17 = ((2.1875e-1 <= T) & (T < 2.34375e-1));
t50 = ((7.34375e-1 <= T) & (T < 7.5e-1));
t66 = ((9.84375e-1 <= T) & (T <= 1.0));
t40 = ((5.78125e-1 <= T) & (T < 5.9375e-1));
t61 = ((9.0625e-1 <= T) & (T < 9.21875e-1));
t2 = ((0.0 <= T) & (T < 7.8125e-3));
t12 = ((T < 1.5625e-1) & (1.40625e-1 <= T));
t51 = ((7.5e-1 <= T) & (T < 7.65625e-1));
t32 = ((4.53125e-1 <= T) & (T < 4.6875e-1));
t41 = ((T < 6.09375e-1) & (5.9375e-1 <= T));
t27 = ((3.75e-1 <= T) & (T < 3.90625e-1));
t62 = ((9.21875e-1 <= T) & (T < 9.375e-1));
t22 = ((T < 3.125e-1) & (2.96875e-1 <= T));
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*-3.597970554028411e+3+1.483850895135036e+1;
elseif (t3)
    out1 = T.*2.393844126169429e+3-3.197254323769527e+1;
elseif (t4)
    out1 = T.*-4.407408623151084e+2+1.231784720737564e+1;
elseif (t5)
    out1 = T.*1.180961580082468e+2-5.145809677729214;
elseif (t6)
    out1 = T.*-3.164377008543042e+1+1.873249451661908;
elseif (t7)
    out1 = T.*8.478922608798549-6.344188417274028e-1;
elseif (t8)
    out1 = T.*-2.271920448610811+2.054907721327035e-1;
elseif (t9)
    out1 = T.*6.087592350326979e-1-6.457294820887549e-2;
elseif (t10)
    out1 = T.*-1.631165404628488e-1+1.985096473594992e-2;
elseif (t11)
    out1 = T.*4.370694750863024e-2-6.001971260484952e-3;
elseif (t12)
    out1 = T.*-1.171126264859121e-2+1.791214542874316e-3;
elseif (t13)
    out1 = T.*3.138116164458176e-3-5.290008966646517e-4;
elseif (t14)
    out1 = T.*-8.410319923864861e-4+1.549151927930246e-4;
elseif (t15)
    out1 = T.*2.254886605787994e-4-4.505742963796642e-5;
elseif (t16)
    out1 = T.*-6.043219597561715e-5+1.302024434964944e-5;
elseif (t17)
    out1 = T.*1.616414288905821e-5-3.735204776998294e-6;
elseif (t18)
    out1 = T.*-4.323405681898006e-6+1.066564419319568e-6;
elseif (t19)
    out1 = T.*1.163115883859574e-6-3.050659721198267e-7;
elseif (t20)
    out1 = T.*-3.457200101902901e-7+9.571856223716831e-8;
elseif (t21)
    out1 = T.*2.688105651124155e-7-7.711816206671762e-8;
elseif (t22)
    out1 = T.*-4.297986780320737e-7+1.302814569918026e-7;
elseif (t23)
    out1 = T.*4.20446631195797e-7-1.35420202141907e-7;
elseif (t24)
    out1 = T.*-2.057466889029508e-7+7.004948101549463e-8;
elseif (t25)
    out1 = T.*5.383801750090372e-8-1.918276181083036e-8;
elseif (t26)
    out1 = T.*-9.627895542318784e-9+3.625300689077725e-9;
elseif (t27)
    out1 = T.*-1.532638541405424e-8+5.762234390978522e-9;
elseif (t28)
    out1 = T.*3.823627830190562e-8-1.51606811230683e-8;
elseif (t29)
    out1 = T.*-3.952568814094396e-8+1.643011774433934e-8;
elseif (t30)
    out1 = T.*2.177020884797828e-8-9.429088797862227e-9;
elseif (t31)
    out1 = T.*-1.485316516406715e-8+6.593637332407647e-9;
elseif (t32)
    out1 = T.*2.129226422813501e-8-9.784760360933956e-9;
elseif (t33)
    out1 = T.*-2.12701747839517e-8+1.016638292598169e-8;
elseif (t34)
    out1 = T.*1.474272031643069e-8-7.277363138266033e-9;
elseif (t35)
    out1 = T.*-2.135201954516246e-8+1.077000679253054e-8;
elseif (t36)
    out1 = T.*7.066501047585642e-8-3.667627431205731e-8;
elseif (t37)
    out1 = T.*-9.631311305959455e-8+5.203085381615101e-8;
elseif (t38)
    out1 = T.*-6.011601202378281e-9+2.647214519235865e-9;
elseif (t39)
    out1 = T.*9.218892788555854e-8-5.25905830927286e-8;
elseif (t40)
    out1 = T.*-3.728153877153126e-8+2.225953044340144e-8;
elseif (t41)
    out1 = T.*-1.935609630859637e-8+1.161629898103385e-8;
elseif (t42)
    out1 = T.*1.661447229787827e-8-1.030326626353663e-8;
elseif (t43)
    out1 = T.*-1.222562292565546e-8+7.721793251171948e-9;
elseif (t44)
    out1 = T.*9.402431713268428e-9-6.133679251888666e-9;
elseif (t45)
    out1 = T.*3.019861768941564e-8-1.978117629873528e-8;
elseif (t46)
    out1 = T.*-9.422923929711054e-8+6.381879011408701e-8;
elseif (t47)
    out1 = T.*1.014897733481882e-7-7.07380310795559e-8;
elseif (t48)
    out1 = T.*-5.015271801080734e-8+3.588559565723786e-8;
elseif (t49)
    out1 = T.*1.208516654648007e-8-8.847883868312469e-9;
elseif (t50)
    out1 = T.*1.342715139257363e-9-9.588961161332905e-10;
elseif (t51)
    out1 = T.*-1.698669101012026e-8+1.278815849589992e-8;
elseif (t52)
    out1 = T.*5.827331799279531e-8-4.483278589695731e-8;
elseif (t53)
    out1 = T.*-1.090864433432896e-7+8.591702764685903e-8;
elseif (t54)
    out1 = T.*1.060577558294833e-7-8.55260060689444e-8;
elseif (t55)
    out1 = T.*-5.9467219262074e-8+4.896303619294593e-8;
elseif (t56)
    out1 = T.*3.293157232294778e-8-2.755471308840023e-8;
elseif (t57)
    out1 = T.*-2.317965985077065e-8+1.978913905817469e-8;
elseif (t58)
    out1 = T.*1.073012144241871e-8-9.352079240659913e-9;
elseif (t59)
    out1 = T.*-2.302753424803275e-9+2.051686268159326e-9;
elseif (t60)
    out1 = T.*-4.787611860854946e-9+4.264763312767845e-9;
elseif (t61)
    out1 = T.*2.472251142696188e-8-2.247878591681615e-8;
elseif (t62)
    out1 = T.*-5.159816489504594e-8+4.78793375675348e-8;
elseif (t63)
    out1 = T.*5.088399853721668e-8-4.81976906502114e-8;
elseif (t64)
    out1 = T.*-2.319363814861828e-8+2.240755681597505e-8;
elseif (t65)
    out1 = T.*4.444634255934535e-9-4.367019575935494e-9;
elseif (t66)
    out1 = T.*-7.742800606370333e-10+7.703492044396438e-10;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
