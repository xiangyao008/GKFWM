function out1 = FGKSF1_1D5_M_36(T)
%FGKSF1_1D5_M_36
%    OUT1 = FGKSF1_1D5_M_36(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-05 16:48:14

t35 = ((T < 7.18e+2./7.29e+2) & (9.654016156073769e-1 <= T));
t23 = ((4.6e+1./8.1e+1 <= T) & (T < 6.01229487374892e-1));
t26 = ((2.0./3.0 <= T) & (T < 6.792663719961388e-1));
t34 = ((2.6e+1./2.7e+1 <= T) & (T < 9.654016156073769e-1));
t9 = ((T < 1.42e+2./7.29e+2) & (1.655235482395976e-1 <= T));
t29 = ((5.26e+2./7.29e+2 <= T) & (T < 7.312909617436366e-1));
t18 = ((3.34e+2./7.29e+2 <= T) & (T < 4.971803078798964e-1));
t5 = ((T < 2.6e+1./2.43e+2) & (2.0./2.7e+1 <= T));
t19 = ((T < 1.22e+2./2.43e+2) & (4.971803078798964e-1 <= T));
t32 = ((T < 2.18e+2./2.43e+2) & (8.678555098308185e-1 <= T));
t15 = ((T < 3.801249809480262e-1) & (3.671188335111518e-1 <= T));
t4 = ((T < 2.0./2.7e+1) & (5.497129502616471e-2 <= T));
t31 = ((T < 8.678555098308185e-1) & (8.483462886755068e-1 <= T));
t6 = ((2.6e+1./2.43e+2 <= T) & (T < 1.330081796474115e-1));
t3 = ((2.895900015241579e-2 <= T) & (T < 5.497129502616471e-2));
t10 = ((1.42e+2./7.29e+2 <= T) & (T < 2.110450642686582e-1));
t11 = ((T < 2.630696540161561e-1) & (2.110450642686582e-1 <= T));
t25 = ((T < 2.0./3.0) & (6.142356348117665e-1 <= T));
t28 = ((T < 5.26e+2./7.29e+2) & (6.922725194330133e-1 <= T));
t27 = ((T < 6.922725194330133e-1) & (6.792663719961388e-1 <= T));
t2 = ((0.0 <= T) & (T < 2.895900015241579e-2));
t36 = ((T <= 1.0) & (7.18e+2./7.29e+2 <= T));
t13 = ((T < 3.411065386374028e-1) & (2.89081948889905e-1 <= T));
t16 = ((3.801249809480262e-1 <= T) & (T < 4.451557181323985e-1));
t24 = ((T < 6.142356348117665e-1) & (6.01229487374892e-1 <= T));
t22 = ((T < 4.6e+1./8.1e+1) & (5.231926027536453e-1 <= T));
t21 = ((5.16689529035208e-1 <= T) & (T < 5.231926027536453e-1));
t12 = ((2.630696540161561e-1 <= T) & (T < 2.89081948889905e-1));
t20 = ((1.22e+2./2.43e+2 <= T) & (T < 5.16689529035208e-1));
t33 = ((T < 2.6e+1./2.7e+1) & (2.18e+2./2.43e+2 <= T));
t7 = ((T < 1.460143270842859e-1) & (1.330081796474115e-1 <= T));
t8 = ((T < 1.655235482395976e-1) & (1.460143270842859e-1 <= T));
t30 = ((7.312909617436366e-1 <= T) & (T < 8.483462886755068e-1));
t14 = ((3.411065386374028e-1 <= T) & (T < 3.671188335111518e-1));
t17 = ((T < 3.34e+2./7.29e+2) & (4.451557181323985e-1 <= T));
et1 = T;
et2 = -2.358508331150002e-8;
et3 = et1.*et2;
et4 = 3.110548417017601e-10;
et5 = T;
et6 = 3.268136223187462e-8;
et7 = et5.*et6;
et8 = -1.318365163364721e-9;
et9 = T;
et10 = -3.91485572973419e-8;
et11 = et9.*et10;
et12 = 2.630218534781511e-9;
et13 = T;
et14 = 1.224614096113259e-8;
et15 = et13.*et14;
et16 = -1.176796151031415e-9;
et17 = T;
et18 = -8.99640024878735e-9;
et19 = et17.*et18;
et20 = 1.096068340564958e-9;
et21 = T;
et22 = -3.703763712849857e-8;
et23 = et21.*et22;
et24 = 4.825782212997207e-9;
et25 = T;
et26 = 7.173711585778488e-8;
et27 = et25.*et26;
et28 = -1.105689014805439e-8;
et29 = T;
et30 = -5.02801389699368e-8;
et31 = et29.*et30;
et32 = 9.139838817485274e-9;
et33 = T;
et34 = -2.242026544908155e-8;
et35 = et33.*et34;
et36 = 3.713087047990835e-9;
et37 = T;
et38 = 8.641466187972917e-8;
et39 = et37.*et38;
et40 = -1.925598718479277e-8;
et41 = T;
et42 = -8.300840621238658e-7;
et43 = et41.*et42;
et44 = 2.218470150450815e-7;
et45 = T;
et46 = 1.353882637092319e-6;
et47 = et45.*et46;
et48 = -4.094983346749862e-7;
et49 = T;
et50 = -1.26898147875374e-5;
et51 = et49.*et50;
et52 = 4.380898683511467e-6;
et53 = T;
et54 = 1.410926876527052e-4;
et55 = et53.*et54;
et56 = -5.207555422677626e-5;
et57 = T;
et58 = -8.054964313271432e-5;
et59 = et57.*et58;
et60 = 3.217623254030747e-5;
et61 = T;
et62 = 3.080338183506484e-3;
et63 = et61.*et62;
et64 = -1.374911057863221e-3;
et65 = T;
et66 = -3.387492351326097e-3;
et67 = et65.*et66;
et68 = 1.588402246161583e-3;
et69 = T;
et70 = 3.134924051172926e-1;
et71 = et69.*et70;
et72 = -1.559580427382363e-1;
et73 = T;
et74 = -3.569704589244602e-1;
et75 = et73.*et74;
et76 = 1.806529425008331e-1;
et77 = T;
et78 = 3.875083155844292;
et79 = et77.*et78;
et80 = -2.006004846565783;
et81 = T;
et82 = -1.563614873555483;
et83 = et81.*et82;
et84 = 8.394817310270064e-1;
et85 = T;
et86 = 7.407669021911973;
et87 = et85.*et86;
et88 = -4.255321468868092;
et89 = T;
et90 = -1.143786283924759e+2;
et91 = et89.*et90;
et92 = 6.896619169487047e+1;
et93 = T;
et94 = 8.502832025695749e+1;
et95 = et93.*et94;
et96 = -5.351666199469157e+1;
et97 = T;
et98 = -1.064603447470723e+3;
et99 = et97.*et98;
et100 = 7.129045164904291e+2;
et101 = T;
et102 = 1.145433425982662e+3;
et103 = et101.*et102;
et104 = -7.882992125179419e+2;
et105 = T;
et106 = -2.189821943055242e+2;
et107 = et105.*et106;
et108 = 1.562482264927184e+2;
et109 = T;
et110 = 1.886153439657832e+2;
et111 = et109.*et110;
et112 = -1.378482140157969e+2;
et113 = T;
et114 = -1.035791824925372;
et115 = et113.*et114;
et116 = 8.419474723633726e-1;
et117 = T;
et118 = 2.283358099600979;
et119 = et117.*et118;
et120 = -1.973841047666146;
et121 = T;
et122 = -3.151196406036649e-1;
et123 = et121.*et122;
et124 = 2.81262176343189e-1;
et125 = T;
et126 = 3.237594221086678e-2;
et127 = et125.*et126;
et128 = -3.048283210770774e-2;
et129 = T;
et130 = -3.02981342133999e-1;
et131 = et129.*et130;
et132 = 2.924538120762371e-1;
et133 = T;
et134 = 3.027395799761191e-3;
et135 = et133.*et134;
et136 = -2.967517914989404e-3;
et137 = T;
et138 = -1.411307456410475e-3;
et139 = et137.*et138;
et140 = 1.404209023187901e-3;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t4,t5,t6,t7,t8,t9}))
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
elseif (t32)
    out1 = et123+et124;
elseif (t33)
    out1 = et127+et128;
elseif (t34)
    out1 = et131+et132;
elseif (t35)
    out1 = et135+et136;
elseif (t36)
    out1 = et139+et140;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
