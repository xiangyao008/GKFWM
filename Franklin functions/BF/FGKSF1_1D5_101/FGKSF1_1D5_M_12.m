function out1 = FGKSF1_1D5_M_12(T)
%FGKSF1_1D5_M_12
%    OUT1 = FGKSF1_1D5_M_12(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-05 16:48:09

t2 = ((T < 2.0./2.7e+1) & (0.0 <= T));
t3 = ((T < 2.6e+1./2.43e+2) & (2.0./2.7e+1 <= T));
t4 = ((2.6e+1./2.43e+2 <= T) & (T < 1.42e+2./7.29e+2));
t5 = ((T < 3.34e+2./7.29e+2) & (1.42e+2./7.29e+2 <= T));
t6 = ((T < 1.22e+2./2.43e+2) & (3.34e+2./7.29e+2 <= T));
t7 = ((T < 4.6e+1./8.1e+1) & (1.22e+2./2.43e+2 <= T));
t8 = ((T < 2.0./3.0) & (4.6e+1./8.1e+1 <= T));
t9 = ((2.0./3.0 <= T) & (T < 5.26e+2./7.29e+2));
t10 = ((T < 2.18e+2./2.43e+2) & (5.26e+2./7.29e+2 <= T));
t11 = ((T < 2.6e+1./2.7e+1) & (2.18e+2./2.43e+2 <= T));
t12 = ((T <= 1.0) & (2.6e+1./2.7e+1 <= T));
et1 = T;
et2 = 1.546075771780396e-2;
et3 = et1.*et2;
et4 = -3.817471042192711e-4;
et5 = T;
et6 = -1.478434962655749e-1;
et7 = et5.*et6;
et8 = 1.171486430195694e-2;
et9 = T;
et10 = 1.720311273395778e-1;
et11 = et9.*et10;
et12 = -2.251040406731866e-2;
et13 = T;
et14 = -1.479340869006295e-1;
et15 = et13.*et14;
et16 = 3.981478169689181e-2;
et17 = T;
et18 = 8.052035136847497;
et19 = et17.*et18;
et20 = -3.717098415466173;
et21 = T;
et22 = -2.113799638442729e+1;
et23 = et21.*et22;
et24 = 1.093797913842487e+1;
et25 = T;
et26 = 4.458753710977687e+1;
et27 = et25.*et26;
et28 = -2.638763247556761e+1;
et29 = T;
et30 = -1.187184070959757e+2;
et31 = et29.*et30;
et32 = 8.24829969949341e+1;
et33 = T;
et34 = 2.788999482858988e+1;
et35 = et33.*et34;
et36 = -2.330029437999251e+1;
et37 = T;
et38 = -3.53158122322953e+1;
et39 = et37.*et38;
et40 = 3.34028576334765e+1;
et41 = T;
et42 = 2.450094290381658e+1;
et43 = et41.*et42;
et44 = -2.419846212722382e+1;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t2,t3,t4,t5,t6,t7,t8,t9}))
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
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
