function out1 = FGKSF1_1D5_M_26(T)
%FGKSF1_1D5_M_26
%    OUT1 = FGKSF1_1D5_M_26(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-05 16:48:11

t11 = ((T < 3.34e+2./7.29e+2) & (3.801249809480262e-1 <= T));
t7 = ((T < 1.42e+2./7.29e+2) & (1.655235482395976e-1 <= T));
t20 = ((5.26e+2./7.29e+2 <= T) & (T < 7.312909617436366e-1));
t15 = ((T < 4.6e+1./8.1e+1) & (5.16689529035208e-1 <= T));
t9 = ((T < 3.411065386374028e-1) & (2.630696540161561e-1 <= T));
t12 = ((3.34e+2./7.29e+2 <= T) & (T < 4.971803078798964e-1));
t4 = ((T < 2.6e+1./2.43e+2) & (2.0./2.7e+1 <= T));
t13 = ((T < 1.22e+2./2.43e+2) & (4.971803078798964e-1 <= T));
t23 = ((T < 2.18e+2./2.43e+2) & (8.678555098308185e-1 <= T));
t22 = ((T < 8.678555098308185e-1) & (8.483462886755068e-1 <= T));
t8 = ((1.42e+2./7.29e+2 <= T) & (T < 2.630696540161561e-1));
t10 = ((3.411065386374028e-1 <= T) & (T < 3.801249809480262e-1));
t17 = ((T < 2.0./3.0) & (6.142356348117665e-1 <= T));
t19 = ((T < 5.26e+2./7.29e+2) & (6.922725194330133e-1 <= T));
t5 = ((2.6e+1./2.43e+2 <= T) & (T < 1.460143270842859e-1));
t2 = ((0.0 <= T) & (T < 2.895900015241579e-2));
t26 = ((T <= 1.0) & (7.18e+2./7.29e+2 <= T));
t3 = ((T < 2.0./2.7e+1) & (2.895900015241579e-2 <= T));
t25 = ((2.6e+1./2.7e+1 <= T) & (T < 7.18e+2./7.29e+2));
t16 = ((4.6e+1./8.1e+1 <= T) & (T < 6.142356348117665e-1));
t14 = ((1.22e+2./2.43e+2 <= T) & (T < 5.16689529035208e-1));
t24 = ((T < 2.6e+1./2.7e+1) & (2.18e+2./2.43e+2 <= T));
t18 = ((2.0./3.0 <= T) & (T < 6.922725194330133e-1));
t6 = ((T < 1.655235482395976e-1) & (1.460143270842859e-1 <= T));
t21 = ((7.312909617436366e-1 <= T) & (T < 8.483462886755068e-1));
et1 = T;
et2 = -4.547310121180632e-9;
et3 = et1.*et2;
et4 = 8.520364089602875e-11;
et5 = T;
et6 = 5.554948428417733e-10;
et7 = et5.*et6;
et8 = -6.256848883484415e-11;
et9 = T;
et10 = 7.140703198032577e-9;
et11 = et9.*et10;
et12 = -5.503617003304592e-10;
et13 = T;
et14 = -2.419296120471515e-8;
et15 = et13.*et14;
et16 = 2.802211445642549e-9;
et17 = T;
et18 = 2.392286383119939e-7;
et19 = et17.*et18;
et20 = -3.566111614725599e-8;
et21 = T;
et22 = -5.736679274006073e-7;
et23 = et21.*et22;
et24 = 9.8892407761277e-8;
et25 = T;
et26 = 7.025927093638292e-7;
et27 = et25.*et26;
et28 = -1.497070578361853e-7;
et29 = T;
et30 = -1.993959617763793e-6;
et31 = et29.*et30;
et32 = 5.596740298977389e-7;
et33 = T;
et34 = 1.981178506606184e-5;
et35 = et33.*et34;
et36 = -6.878408061612972e-6;
et37 = T;
et38 = -3.267534068987043e-5;
et39 = et37.*et38;
et40 = 1.307325961637744e-5;
et41 = T;
et42 = 3.069349257652883e-4;
et43 = et41.*et42;
et44 = -1.42523213629196e-4;
et45 = T;
et46 = -3.615079864277961e-2;
et47 = et45.*et46;
et48 = 1.798354398658319e-2;
et49 = T;
et50 = 4.142903238973691e-2;
et51 = et49.*et50;
et52 = -2.096600081163498e-2;
et53 = T;
et54 = -2.975374680844489e-2;
et55 = et53.*et54;
et56 = 1.581339584769078e-2;
et57 = T;
et58 = 1.113833099965459e-1;
et59 = et57.*et58;
et60 = -6.433851295514349e-2;
et61 = T;
et62 = -3.524533942151378e-1;
et63 = et61.*et62;
et64 = 2.205665195053176e-1;
et65 = T;
et66 = 3.664795366056007;
et67 = et65.*et66;
et68 = -2.457599320675445;
et69 = T;
et70 = -1.246342062891947e+1;
et71 = et69.*et70;
et72 = 8.70752140012605;
et73 = T;
et74 = 2.38788967589139e+2;
et75 = et73.*et74;
et76 = -1.725802100164703e+2;
et77 = T;
et78 = -4.856052366532054e+1;
et79 = et77.*et78;
et80 = 3.75558757995481e+1;
et81 = T;
et82 = 5.241219705227748e+2;
et83 = et81.*et82;
et84 = -4.482771927343551e+2;
et85 = T;
et86 = -2.650594247587998e+2;
et87 = et85.*et86;
et88 = 2.366182294167325e+2;
et89 = T;
et90 = 2.475556274138454e+1;
et91 = et89.*et90;
et92 = -2.338040134474971e+1;
et93 = T;
et94 = -2.776928342482187e+1;
et95 = et93.*et94;
et96 = 2.719908014863424e+1;
et97 = T;
et98 = 1.502946418414517e+1;
et99 = et97.*et98;
et100 = -1.49538701713086e+1;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t3,t4,t5,t6,t7,t8,t9}))
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
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
