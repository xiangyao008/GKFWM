function out1 = FGKSF1_1D5_M_31(T)
%FGKSF1_1D5_M_31
%    OUT1 = FGKSF1_1D5_M_31(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-05 16:48:13

t30 = ((T < 7.18e+2./7.29e+2) & (9.654016156073769e-1 <= T));
t15 = ((T < 3.34e+2./7.29e+2) & (3.801249809480262e-1 <= T));
t29 = ((2.6e+1./2.7e+1 <= T) & (T < 9.654016156073769e-1));
t9 = ((T < 1.42e+2./7.29e+2) & (1.655235482395976e-1 <= T));
t24 = ((5.26e+2./7.29e+2 <= T) & (T < 7.312909617436366e-1));
t19 = ((T < 4.6e+1./8.1e+1) & (5.16689529035208e-1 <= T));
t16 = ((3.34e+2./7.29e+2 <= T) & (T < 4.971803078798964e-1));
t5 = ((T < 2.6e+1./2.43e+2) & (2.0./2.7e+1 <= T));
t17 = ((T < 1.22e+2./2.43e+2) & (4.971803078798964e-1 <= T));
t27 = ((T < 2.18e+2./2.43e+2) & (8.678555098308185e-1 <= T));
t4 = ((T < 2.0./2.7e+1) & (5.497129502616471e-2 <= T));
t26 = ((T < 8.678555098308185e-1) & (8.483462886755068e-1 <= T));
t6 = ((2.6e+1./2.43e+2 <= T) & (T < 1.330081796474115e-1));
t3 = ((2.895900015241579e-2 <= T) & (T < 5.497129502616471e-2));
t14 = ((3.411065386374028e-1 <= T) & (T < 3.801249809480262e-1));
t10 = ((1.42e+2./7.29e+2 <= T) & (T < 2.110450642686582e-1));
t11 = ((T < 2.630696540161561e-1) & (2.110450642686582e-1 <= T));
t21 = ((T < 2.0./3.0) & (6.142356348117665e-1 <= T));
t23 = ((T < 5.26e+2./7.29e+2) & (6.922725194330133e-1 <= T));
t2 = ((0.0 <= T) & (T < 2.895900015241579e-2));
t31 = ((T <= 1.0) & (7.18e+2./7.29e+2 <= T));
t13 = ((T < 3.411065386374028e-1) & (2.89081948889905e-1 <= T));
t12 = ((2.630696540161561e-1 <= T) & (T < 2.89081948889905e-1));
t20 = ((4.6e+1./8.1e+1 <= T) & (T < 6.142356348117665e-1));
t18 = ((1.22e+2./2.43e+2 <= T) & (T < 5.16689529035208e-1));
t28 = ((T < 2.6e+1./2.7e+1) & (2.18e+2./2.43e+2 <= T));
t7 = ((T < 1.460143270842859e-1) & (1.330081796474115e-1 <= T));
t22 = ((2.0./3.0 <= T) & (T < 6.922725194330133e-1));
t8 = ((T < 1.655235482395976e-1) & (1.460143270842859e-1 <= T));
t25 = ((7.312909617436366e-1 <= T) & (T < 8.483462886755068e-1));
et1 = T;
et2 = -1.534475339142464e-3;
et3 = et1.*et2;
et4 = 1.481231608449883e-5;
et5 = T;
et6 = 5.31841629805609e-3;
et7 = et5.*et6;
et8 = -1.83640573881623e-4;
et9 = T;
et10 = -3.046187823865001e-2;
et11 = et9.*et10;
et12 = 1.783248553218718e-3;
et13 = T;
et14 = 5.788280548142398e-2;
et15 = et13.*et14;
et16 = -4.760802092712689e-3;
et17 = T;
et18 = -2.815666813235241e-1;
et19 = et17.*et18;
et20 = 3.155889608394842e-2;
et21 = T;
et22 = 2.950726964357282;
et23 = et21.*et22;
et24 = -3.983625978139508e-1;
et25 = T;
et26 = -7.014324689479607;
et27 = et25.*et26;
et28 = 1.056677713781193;
et29 = T;
et30 = 1.471309873069872e+1;
et31 = et29.*et30;
et32 = -2.539722504830857;
et33 = T;
et34 = -1.208719499695619e+2;
et35 = et33.*et34;
et36 = 2.387053389494555e+1;
et37 = T;
et38 = 1.122364291864564e+2;
et39 = et37.*et38;
et40 = -2.532583896559907e+1;
et41 = T;
et42 = -3.813883675554435e+2;
et43 = et41.*et42;
et44 = 1.045318655270879e+2;
et45 = T;
et46 = 1.658931328002589e+2;
et47 = et45.*et46;
et48 = -5.367733718712979e+1;
et49 = T;
et50 = -8.870307403730135e+1;
et51 = et49.*et50;
et52 = 3.316709367744267e+1;
et53 = T;
et54 = 9.60693195886924;
et55 = et53.*et54;
et56 = -4.202995478852022;
et57 = T;
et58 = -7.366245537801439;
et59 = et57.*et58;
et60 = 3.573467187661019;
et61 = T;
et62 = 2.061396103429012e+1;
et63 = et61.*et62;
et64 = -1.033774053039456e+1;
et65 = T;
et66 = -8.95894426282883e-1;
et67 = et65.*et66;
et68 = 4.614461617449707e-1;
et69 = T;
et70 = 3.63103336922932e-2;
et71 = et69.*et70;
et72 = -2.021427665098229e-2;
et73 = T;
et74 = -1.103989152311552e-2;
et75 = et73.*et74;
et76 = 6.675974705916499e-3;
et77 = T;
et78 = 2.708391784259558e-3;
et79 = et77.*et78;
et80 = -1.768710818961047e-3;
et81 = T;
et82 = -1.814142104840504e-3;
et83 = et81.*et82;
et84 = 1.246311773772328e-3;
et85 = T;
et86 = 4.503555978310867e-4;
et87 = et85.*et86;
et88 = -3.213377561064051e-4;
et89 = T;
et90 = -3.879151315761479e-4;
et91 = et89.*et90;
et92 = 2.835050472793362e-4;
et93 = T;
et94 = 2.130645926813845e-6;
et95 = et93.*et94;
et96 = -1.731904474849196e-6;
et97 = T;
et98 = -4.697419467270469e-6;
et99 = et97.*et98;
et100 = 4.060659461055894e-6;
et101 = T;
et102 = 6.501085664660044e-7;
et103 = et101.*et102;
et104 = -5.802222069970679e-7;
et105 = T;
et106 = -6.879341912460081e-8;
et107 = et105.*et106;
et108 = 6.47186689648742e-8;
et109 = T;
et110 = 7.702507497097196e-7;
et111 = et109.*et110;
et112 = -7.432497899126195e-7;
et113 = T;
et114 = -2.38889854334079e-8;
et115 = et113.*et114;
et116 = 2.341399341257016e-8;
et117 = T;
et118 = 1.104824631113294e-8;
et119 = et117.*et118;
et120 = -1.099606473911753e-8;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t4,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = et3+et4;
elseif (t3)
    out1 = et7+et8;
elseif (t4)
    out1 = et11+et12;
elseif (t5)
    out1 = et15+et16;
elseif (t6)
    out1 = et19+et20;
elseif (t7)
    out1 = et23+et24;
elseif (t8)
    out1 = et27+et28;
elseif (t9)
    out1 = et31+et32;
elseif (t10)
    out1 = et35+et36;
elseif (t11)
    out1 = et39+et40;
elseif (t12)
    out1 = et43+et44;
elseif (t13)
    out1 = et47+et48;
elseif (t14)
    out1 = et51+et52;
elseif (t15)
    out1 = et55+et56;
elseif (t16)
    out1 = et59+et60;
elseif (t17)
    out1 = et63+et64;
elseif (t18)
    out1 = et67+et68;
elseif (t19)
    out1 = et71+et72;
elseif (t20)
    out1 = et75+et76;
elseif (t21)
    out1 = et79+et80;
elseif (t22)
    out1 = et83+et84;
elseif (t23)
    out1 = et87+et88;
elseif (t24)
    out1 = et91+et92;
elseif (t25)
    out1 = et95+et96;
elseif (t26)
    out1 = et99+et100;
elseif (t27)
    out1 = et103+et104;
elseif (t28)
    out1 = et107+et108;
elseif (t29)
    out1 = et111+et112;
elseif (t30)
    out1 = et115+et116;
elseif (t31)
    out1 = et119+et120;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
