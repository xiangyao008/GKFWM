function out1 = FGKSF1_1D5_M_28(T)
%FGKSF1_1D5_M_28
%    OUT1 = FGKSF1_1D5_M_28(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-05 16:48:12

t27 = ((T < 7.18e+2./7.29e+2) & (9.654016156073769e-1 <= T));
t12 = ((T < 3.34e+2./7.29e+2) & (3.801249809480262e-1 <= T));
t26 = ((2.6e+1./2.7e+1 <= T) & (T < 9.654016156073769e-1));
t8 = ((T < 1.42e+2./7.29e+2) & (1.655235482395976e-1 <= T));
t21 = ((5.26e+2./7.29e+2 <= T) & (T < 7.312909617436366e-1));
t16 = ((T < 4.6e+1./8.1e+1) & (5.16689529035208e-1 <= T));
t10 = ((T < 3.411065386374028e-1) & (2.630696540161561e-1 <= T));
t13 = ((3.34e+2./7.29e+2 <= T) & (T < 4.971803078798964e-1));
t5 = ((T < 2.6e+1./2.43e+2) & (2.0./2.7e+1 <= T));
t14 = ((T < 1.22e+2./2.43e+2) & (4.971803078798964e-1 <= T));
t24 = ((T < 2.18e+2./2.43e+2) & (8.678555098308185e-1 <= T));
t4 = ((T < 2.0./2.7e+1) & (5.497129502616471e-2 <= T));
t23 = ((T < 8.678555098308185e-1) & (8.483462886755068e-1 <= T));
t9 = ((1.42e+2./7.29e+2 <= T) & (T < 2.630696540161561e-1));
t3 = ((2.895900015241579e-2 <= T) & (T < 5.497129502616471e-2));
t11 = ((3.411065386374028e-1 <= T) & (T < 3.801249809480262e-1));
t18 = ((T < 2.0./3.0) & (6.142356348117665e-1 <= T));
t20 = ((T < 5.26e+2./7.29e+2) & (6.922725194330133e-1 <= T));
t6 = ((2.6e+1./2.43e+2 <= T) & (T < 1.460143270842859e-1));
t2 = ((0.0 <= T) & (T < 2.895900015241579e-2));
t28 = ((T <= 1.0) & (7.18e+2./7.29e+2 <= T));
t17 = ((4.6e+1./8.1e+1 <= T) & (T < 6.142356348117665e-1));
t15 = ((1.22e+2./2.43e+2 <= T) & (T < 5.16689529035208e-1));
t25 = ((T < 2.6e+1./2.7e+1) & (2.18e+2./2.43e+2 <= T));
t19 = ((2.0./3.0 <= T) & (T < 6.922725194330133e-1));
t7 = ((T < 1.655235482395976e-1) & (1.460143270842859e-1 <= T));
t22 = ((7.312909617436366e-1 <= T) & (T < 8.483462886755068e-1));
et1 = T;
et2 = 2.223243563910894e+2;
et3 = et1.*et2;
et4 = -2.146097023694175;
et5 = T;
et6 = -4.604737822368706e+2;
et7 = et5.*et6;
et8 = 1.762705437690214e+1;
et9 = T;
et10 = 6.38844004439532e+2;
et11 = et9.*et10;
et12 = -4.280386800199679e+1;
et13 = T;
et14 = -1.718776194061006e+2;
et15 = et13.*et14;
et16 = 1.724958561619822e+1;
et17 = T;
et18 = 3.93309836460504e+1;
et19 = et17.*et18;
et20 = -5.348865739175967;
et21 = T;
et22 = -2.444612767916653e+1;
et23 = et21.*et22;
et24 = 3.963506254355174;
et25 = T;
et26 = 3.299901385136275;
et27 = et25.*et26;
et28 = -6.29114925927229e-1;
et29 = T;
et30 = -2.515258425145665e-1;
et31 = et29.*et30;
et32 = 6.265827891010919e-2;
et33 = T;
et34 = 6.04696260286535e-2;
et35 = et33.*et34;
et36 = -1.941826105416421e-2;
et37 = T;
et38 = -3.683364546725958e-2;
et39 = et37.*et38;
et40 = 1.377252108390217e-2;
et41 = T;
et42 = 3.989253403729201e-3;
et43 = et41.*et42;
et44 = -1.745282571675647e-3;
et45 = T;
et46 = -3.05887614476472e-3;
et47 = et45.*et46;
et48 = 1.483901611036245e-3;
et49 = T;
et50 = 8.561144443559184e-3;
et51 = et49.*et50;
et52 = -4.293343802637368e-3;
et53 = T;
et54 = -3.723270565707845e-4;
et55 = et53.*et54;
et56 = 1.917735760287069e-4;
et57 = T;
et58 = 1.508952903703186e-5;
et59 = et57.*et58;
et60 = -8.400517129424112e-6;
et61 = T;
et62 = -4.584891691527973e-6;
et63 = et61.*et62;
et64 = 2.772610691733326e-6;
et65 = T;
et66 = 1.122381918645662e-6;
et67 = et65.*et66;
et68 = -7.330001372561201e-7;
et69 = T;
et70 = -7.533362739156212e-7;
et71 = et69.*et70;
et72 = 5.174786577847357e-7;
et73 = T;
et74 = 1.950723332612815e-7;
et75 = et73.*et74;
et76 = -1.390785581575737e-7;
et77 = T;
et78 = -1.851287439596377e-7;
et79 = et77.*et78;
et80 = 1.352503398097836e-7;
et81 = T;
et82 = 1.963147895710507e-9;
et83 = et81.*et82;
et84 = -1.568269719550423e-9;
et85 = T;
et86 = -1.083497730617847e-8;
et87 = et85.*et86;
et88 = 9.288972297476562e-9;
et89 = T;
et90 = 8.535789409705118e-9;
et91 = et89.*et90;
et92 = -7.522054326550439e-9;
et93 = T;
et94 = -8.048053749346686e-9;
et95 = et93.*et94;
et96 = 7.355632128895212e-9;
et97 = T;
et98 = 6.793652293733749e-7;
et99 = et97.*et98;
et100 = -6.545978997670589e-7;
et101 = T;
et102 = -8.850085879311541e-8;
et103 = et101.*et102;
et104 = 8.670126231898737e-8;
et105 = T;
et106 = 4.729866147959172e-8;
et107 = et105.*et106;
et108 = -4.704915682477631e-8;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t3,t4,t5,t6,t7,t8,t9}))
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
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
