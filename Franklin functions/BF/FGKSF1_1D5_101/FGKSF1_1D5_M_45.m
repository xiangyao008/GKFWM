function out1 = FGKSF1_1D5_M_45(T)
%FGKSF1_1D5_M_45
%    OUT1 = FGKSF1_1D5_M_45(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-06-05 16:48:18

t13 = ((1.42e+2./7.29e+2 <= T) & (T < 2.023742993107419e-1));
t43 = ((T < 7.18e+2./7.29e+2) & (9.654016156073769e-1 <= T));
t28 = ((4.6e+1./8.1e+1 <= T) & (T < 6.01229487374892e-1));
t11 = ((T < 1.655235482395976e-1) & (1.503497095632441e-1 <= T));
t31 = ((2.0./3.0 <= T) & (T < 6.792663719961388e-1));
t14 = ((T < 2.110450642686582e-1) & (2.023742993107419e-1 <= T));
t42 = ((2.6e+1./2.7e+1 <= T) & (T < 9.654016156073769e-1));
t12 = ((T < 1.42e+2./7.29e+2) & (1.655235482395976e-1 <= T));
t37 = ((T < 8.483462886755068e-1) & (8.353401412386323e-1 <= T));
t10 = ((1.460143270842859e-1 <= T) & (T < 1.503497095632441e-1));
t34 = ((5.26e+2./7.29e+2 <= T) & (T < 7.312909617436366e-1));
t6 = ((2.0./2.7e+1 <= T) & (T < 9.832511981574624e-2));
t23 = ((3.34e+2./7.29e+2 <= T) & (T < 4.971803078798964e-1));
t41 = ((T < 2.6e+1./2.7e+1) & (9.133770258598791e-1 <= T));
t4 = ((T < 5.497129502616471e-2) & (4.63005300682484e-2 <= T));
t24 = ((T < 1.22e+2./2.43e+2) & (4.971803078798964e-1 <= T));
t39 = ((T < 2.18e+2./2.43e+2) & (8.678555098308185e-1 <= T));
t20 = ((T < 3.801249809480262e-1) & (3.671188335111518e-1 <= T));
t5 = ((T < 2.0./2.7e+1) & (5.497129502616471e-2 <= T));
t38 = ((T < 8.678555098308185e-1) & (8.483462886755068e-1 <= T));
t8 = ((2.6e+1./2.43e+2 <= T) & (T < 1.330081796474115e-1));
t40 = ((2.18e+2./2.43e+2 <= T) & (T < 9.133770258598791e-1));
t30 = ((T < 2.0./3.0) & (6.142356348117665e-1 <= T));
t33 = ((T < 5.26e+2./7.29e+2) & (6.922725194330133e-1 <= T));
t45 = ((T <= 1.0) & (9.914139104811258e-1 <= T));
t32 = ((T < 6.922725194330133e-1) & (6.792663719961388e-1 <= T));
t2 = ((0.0 <= T) & (T < 2.895900015241579e-2));
t44 = ((7.18e+2./7.29e+2 <= T) & (T < 9.914139104811258e-1));
t3 = ((2.895900015241579e-2 <= T) & (T < 4.63005300682484e-2));
t35 = ((7.312909617436366e-1 <= T) & (T < 7.573032566173856e-1));
t18 = ((T < 3.411065386374028e-1) & (2.89081948889905e-1 <= T));
t21 = ((3.801249809480262e-1 <= T) & (T < 4.451557181323985e-1));
t29 = ((T < 6.142356348117665e-1) & (6.01229487374892e-1 <= T));
t27 = ((T < 4.6e+1./8.1e+1) & (5.231926027536453e-1 <= T));
t16 = ((T < 2.630696540161561e-1) & (2.543988890582398e-1 <= T));
t36 = ((7.573032566173856e-1 <= T) & (T < 8.353401412386323e-1));
t26 = ((5.16689529035208e-1 <= T) & (T < 5.231926027536453e-1));
t17 = ((2.630696540161561e-1 <= T) & (T < 2.89081948889905e-1));
t15 = ((2.110450642686582e-1 <= T) & (T < 2.543988890582398e-1));
t25 = ((1.22e+2./2.43e+2 <= T) & (T < 5.16689529035208e-1));
t9 = ((T < 1.460143270842859e-1) & (1.330081796474115e-1 <= T));
t7 = ((T < 2.6e+1./2.43e+2) & (9.832511981574624e-2 <= T));
t19 = ((3.411065386374028e-1 <= T) & (T < 3.671188335111518e-1));
t22 = ((T < 3.34e+2./7.29e+2) & (4.451557181323985e-1 <= T));
et1 = T;
et2 = 5.192124021236523e-6;
et3 = et1.*et2;
et4 = -4.989139272674972e-8;
et5 = T;
et6 = -3.204733242214273e-5;
et7 = et5.*et6;
et8 = 1.028526032092951e-6;
et9 = T;
et10 = 3.45306481313865e-4;
et11 = et9.*et10;
et12 = -1.644315556715928e-5;
et13 = T;
et14 = -5.087179840975218e-4;
et15 = et13.*et14;
et16 = 3.050367528053266e-5;
et17 = T;
et18 = 1.272033261281719e-3;
et19 = et17.*et18;
et20 = -1.014038243771889e-4;
et21 = T;
et22 = -2.114311058464576e-2;
et23 = et21.*et22;
et24 = 2.102567879960818e-3;
et25 = T;
et26 = 2.220191875910874e-2;
et27 = et25.*et26;
et28 = -2.535171885214561e-3;
et29 = T;
et30 = -2.00346365617113e-1;
et31 = et29.*et30;
et32 = 2.706557030332115e-2;
et33 = T;
et34 = 4.252729936377674;
et35 = et33.*et34;
et36 = -6.231473693874281e-1;
et37 = T;
et38 = -3.78337053406969;
et39 = et37.*et38;
et40 = 5.850780023653822e-1;
et41 = T;
et42 = 5.390082619753577;
et43 = et41.*et42;
et44 = -9.333445132651722e-1;
et45 = T;
et46 = -1.437036543900883e+2;
et47 = et45.*et46;
et48 = 2.810823388920059e+1;
et49 = T;
et50 = 5.216420552955357e+2;
et51 = et49.*et50;
et52 = -1.065406379078359e+2;
et53 = T;
et54 = -2.325455642533856e+2;
et55 = et53.*et54;
et56 = 5.262693673049257e+1;
et57 = T;
et58 = 1.549549379466676e+3;
et59 = et57.*et58;
et60 = -4.007360371481974e+2;
et61 = T;
et62 = -3.155250116293835e+2;
et63 = et61.*et62;
et64 = 8.990843763183593e+1;
et65 = T;
et66 = 3.392433516316469e+1;
et67 = et65.*et66;
et68 = -1.111106057725818e+1;
et69 = T;
et70 = -2.371833997664497e+1;
et71 = et69.*et70;
et72 = 8.551232817482567;
et73 = T;
et74 = 1.322539848618911e+1;
et75 = et73.*et74;
et76 = -5.01150935253415;
et77 = T;
et78 = -3.45370312252886e-1;
et79 = et77.*et78;
et80 = 1.47078878423681e-1;
et81 = T;
et82 = 5.894142000476282e-1;
et83 = et81.*et82;
et84 = -2.690457924484984e-1;
et85 = T;
et86 = -3.716019946646453e-2;
et87 = et85.*et86;
et88 = 1.802670335082535e-2;
et89 = T;
et90 = 1.052157624830436e-1;
et91 = et89.*et90;
et92 = -5.275982124592749e-2;
et93 = T;
et94 = -5.956434530290257e-3;
et95 = et93.*et94;
et96 = 3.055026637310066e-3;
et97 = T;
et98 = 3.730671951184018e-3;
et99 = et97.*et98;
et100 = -1.950199848316789e-3;
et101 = T;
et102 = -4.86817031485791e-5;
et103 = et101.*et102;
et104 = 2.713002681998304e-5;
et105 = T;
et106 = 2.120157923498954e-5;
et107 = et105.*et106;
et108 = -1.255677552130285e-5;
et109 = T;
et110 = -1.643636186887979e-5;
et111 = et109.*et110;
et112 = 1.007226451442289e-5;
et113 = T;
et114 = 6.369921468220249e-7;
et115 = et113.*et114;
et116 = -4.147979277777414e-7;
et117 = T;
et118 = -1.001951119095317e-6;
et119 = et117.*et118;
et120 = 6.778309161671533e-7;
et121 = T;
et122 = 2.330671117514297e-7;
et123 = et121.*et122;
et124 = -1.610754368492063e-7;
et125 = T;
et126 = 8.77062652735712e-9;
et127 = et125.*et126;
et128 = -5.801143923167954e-9;
et129 = T;
et130 = -1.278421880541785e-7;
et131 = et129.*et130;
et132 = 9.276996783250797e-8;
et133 = T;
et134 = 3.003320043115146e-8;
et135 = et133.*et134;
et136 = -2.268287684857925e-8;
et137 = T;
et138 = 1.345722418820241e-10;
et139 = et137.*et138;
et140 = -4.054835245314251e-11;
et141 = T;
et142 = -4.001682038132702e-9;
et143 = et141.*et142;
et144 = 3.414630882013256e-9;
et145 = T;
et146 = -1.943863236951488e-8;
et147 = et145.*et146;
et148 = 1.651051040410947e-8;
et149 = T;
et150 = 2.724434078844952e-8;
et151 = et149.*et150;
et152 = -2.400356506631415e-8;
et153 = T;
et154 = -4.438078724602839e-8;
et155 = et153.*et154;
et156 = 4.025272263539855e-8;
et157 = T;
et158 = 1.327655930613306e-8;
et159 = et157.*et158;
et160 = -1.241017307738704e-8;
et161 = T;
et162 = -4.360434971007703e-7;
et163 = et161.*et162;
et164 = 4.202683997588903e-7;
et165 = T;
et166 = 5.047215123396006e-8;
et167 = et165.*et166;
et168 = -4.941459316172885e-8;
et169 = T;
et170 = -4.909116433758626e-8;
et171 = et169.*et170;
et172 = 4.864639528871046e-8;
et173 = T;
et174 = -5.276521736674765e-9;
et175 = et173.*et174;
et176 = 5.207949131407868e-9;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t5,t6,t7,t8,t9}))
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
elseif (t37)
    out1 = et143+et144;
elseif (t38)
    out1 = et147+et148;
elseif (t39)
    out1 = et151+et152;
elseif (t40)
    out1 = et155+et156;
elseif (t41)
    out1 = et159+et160;
elseif (t42)
    out1 = et163+et164;
elseif (t43)
    out1 = et167+et168;
elseif (t44)
    out1 = et171+et172;
elseif (t45)
    out1 = et175+et176;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
