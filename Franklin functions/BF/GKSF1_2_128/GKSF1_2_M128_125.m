function out1 = GKSF1_2_M128_125(T)
%GKSF1_2_M128_125
%    OUT1 = GKSF1_2_M128_125(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:43:03

t51 = ((3.828125e-1 <= T) & (T < 3.90625e-1));
t90 = ((6.875e-1 <= T) & (T < 6.953125e-1));
t111 = ((8.515625e-1 <= T) & (T < 8.59375e-1));
t14 = ((9.375e-2 <= T) & (T < 1.015625e-1));
t41 = ((3.046875e-1 <= T) & (T < 3.125e-1));
t70 = ((5.3125e-1 <= T) & (T < 5.390625e-1));
t91 = ((6.953125e-1 <= T) & (T < 7.03125e-1));
t62 = ((4.6875e-1 <= T) & (T < 4.765625e-1));
t112 = ((8.59375e-1 <= T) & (T < 8.671875e-1));
t4 = ((T < 2.34375e-2) & (1.5625e-2 <= T));
t5 = ((2.34375e-2 <= T) & (T < 3.125e-2));
t122 = ((9.375e-1 <= T) & (T < 9.53125e-1));
t71 = ((5.390625e-1 <= T) & (T < 5.46875e-1));
t52 = ((3.90625e-1 <= T) & (T < 3.984375e-1));
t92 = ((7.03125e-1 <= T) & (T < 7.109375e-1));
t32 = ((2.34375e-1 <= T) & (T < 2.421875e-1));
t113 = ((8.671875e-1 <= T) & (T < 8.75e-1));
t42 = ((3.125e-1 <= T) & (T < 3.203125e-1));
t72 = ((5.46875e-1 <= T) & (T < 5.546875e-1));
t27 = ((1.953125e-1 <= T) & (T < 2.03125e-1));
t93 = ((7.109375e-1 <= T) & (T < 7.1875e-1));
t63 = ((4.765625e-1 <= T) & (T < 4.84375e-1));
t114 = ((8.75e-1 <= T) & (T < 8.828125e-1));
t123 = ((9.53125e-1 <= T) & (T < 9.6875e-1));
t22 = ((1.5625e-1 <= T) & (T < 1.640625e-1));
t73 = ((5.546875e-1 <= T) & (T < 5.625e-1));
t53 = ((3.984375e-1 <= T) & (T < 4.0625e-1));
t94 = ((7.1875e-1 <= T) & (T < 7.265625e-1));
t115 = ((8.828125e-1 <= T) & (T < 8.90625e-1));
t8 = ((4.6875e-2 <= T) & (T < 5.46875e-2));
t43 = ((3.203125e-1 <= T) & (T < 3.28125e-1));
t74 = ((5.625e-1 <= T) & (T < 5.703125e-1));
t95 = ((7.265625e-1 <= T) & (T < 7.34375e-1));
t64 = ((4.84375e-1 <= T) & (T < 4.921875e-1));
t116 = ((8.90625e-1 <= T) & (T < 8.984375e-1));
t124 = ((9.6875e-1 <= T) & (T < 9.84375e-1));
t12 = ((7.8125e-2 <= T) & (T < 8.59375e-2));
t75 = ((5.703125e-1 <= T) & (T < 5.78125e-1));
t3 = ((T < 1.5625e-2) & (7.8125e-3 <= T));
t54 = ((4.0625e-1 <= T) & (T < 4.140625e-1));
t96 = ((7.34375e-1 <= T) & (T < 7.421875e-1));
t117 = ((8.984375e-1 <= T) & (T < 9.0625e-1));
t44 = ((3.28125e-1 <= T) & (T < 3.359375e-1));
t76 = ((5.78125e-1 <= T) & (T < 5.859375e-1));
t28 = ((2.03125e-1 <= T) & (T < 2.109375e-1));
t97 = ((7.421875e-1 <= T) & (T < 7.5e-1));
t118 = ((9.0625e-1 <= T) & (T < 9.140625e-1));
t34 = ((2.5e-1 <= T) & (T < 2.578125e-1));
t23 = ((1.640625e-1 <= T) & (T < 1.71875e-1));
t77 = ((5.859375e-1 <= T) & (T < 5.9375e-1));
t55 = ((4.140625e-1 <= T) & (T < 4.21875e-1));
t98 = ((7.5e-1 <= T) & (T < 7.578125e-1));
t119 = ((9.140625e-1 <= T) & (T < 9.21875e-1));
t18 = ((1.25e-1 <= T) & (T < 1.328125e-1));
t15 = ((1.015625e-1 <= T) & (T < 1.09375e-1));
t45 = ((3.359375e-1 <= T) & (T < 3.4375e-1));
t78 = ((T < 6.015625e-1) & (5.9375e-1 <= T));
t99 = ((7.578125e-1 <= T) & (T < 7.65625e-1));
t120 = ((9.21875e-1 <= T) & (T < 9.296875e-1));
t7 = ((3.90625e-2 <= T) & (T < 4.6875e-2));
t9 = ((5.46875e-2 <= T) & (T < 6.25e-2));
t35 = ((2.578125e-1 <= T) & (T < 2.65625e-1));
t79 = ((6.015625e-1 <= T) & (T < 6.09375e-1));
t56 = ((4.21875e-1 <= T) & (T < 4.296875e-1));
t100 = ((7.65625e-1 <= T) & (T < 7.734375e-1));
t121 = ((9.296875e-1 <= T) & (T < 9.375e-1));
t10 = ((6.25e-2 <= T) & (T < 7.03125e-2));
t46 = ((3.4375e-1 <= T) & (T < 3.515625e-1));
t80 = ((6.09375e-1 <= T) & (T < 6.171875e-1));
t29 = ((2.109375e-1 <= T) & (T < 2.1875e-1));
t101 = ((7.734375e-1 <= T) & (T < 7.8125e-1));
t36 = ((2.65625e-1 <= T) & (T < 2.734375e-1));
t24 = ((1.71875e-1 <= T) & (T < 1.796875e-1));
t81 = ((6.171875e-1 <= T) & (T < 6.25e-1));
t57 = ((4.296875e-1 <= T) & (T < 4.375e-1));
t102 = ((7.8125e-1 <= T) & (T < 7.890625e-1));
t19 = ((1.328125e-1 <= T) & (T < 1.40625e-1));
t47 = ((3.515625e-1 <= T) & (T < 3.59375e-1));
t82 = ((6.25e-1 <= T) & (T < 6.328125e-1));
t103 = ((7.890625e-1 <= T) & (T < 7.96875e-1));
t13 = ((8.59375e-2 <= T) & (T < 9.375e-2));
t37 = ((2.734375e-1 <= T) & (T < 2.8125e-1));
t17 = ((1.171875e-1 <= T) & (T < 1.25e-1));
t83 = ((6.328125e-1 <= T) & (T < 6.40625e-1));
t58 = ((4.375e-1 <= T) & (T < 4.453125e-1));
t104 = ((7.96875e-1 <= T) & (T < 8.046875e-1));
t6 = ((T < 3.90625e-2) & (3.125e-2 <= T));
t48 = ((3.59375e-1 <= T) & (T < 3.671875e-1));
t84 = ((6.40625e-1 <= T) & (T < 6.484375e-1));
t30 = ((2.1875e-1 <= T) & (T < 2.265625e-1));
t105 = ((8.046875e-1 <= T) & (T < 8.125e-1));
t38 = ((2.8125e-1 <= T) & (T < 2.890625e-1));
t25 = ((1.796875e-1 <= T) & (T < 1.875e-1));
t33 = ((2.421875e-1 <= T) & (T < 2.5e-1));
t85 = ((6.484375e-1 <= T) & (T < 6.5625e-1));
t59 = ((4.453125e-1 <= T) & (T < 4.53125e-1));
t106 = ((8.125e-1 <= T) & (T < 8.203125e-1));
t20 = ((1.40625e-1 <= T) & (T < 1.484375e-1));
t16 = ((1.09375e-1 <= T) & (T < 1.171875e-1));
t49 = ((3.671875e-1 <= T) & (T < 3.75e-1));
t86 = ((6.5625e-1 <= T) & (T < 6.640625e-1));
t65 = ((4.921875e-1 <= T) & (T < 5.0e-1));
t107 = ((8.203125e-1 <= T) & (T < 8.28125e-1));
t39 = ((2.890625e-1 <= T) & (T < 2.96875e-1));
t66 = ((5.0e-1 <= T) & (T < 5.078125e-1));
t125 = ((9.84375e-1 <= T) & (T <= 1.0));
t87 = ((6.640625e-1 <= T) & (T < 6.71875e-1));
t60 = ((4.53125e-1 <= T) & (T < 4.609375e-1));
t108 = ((8.28125e-1 <= T) & (T < 8.359375e-1));
t2 = ((0.0 <= T) & (T < 7.8125e-3));
t11 = ((T < 7.8125e-2) & (7.03125e-2 <= T));
t67 = ((5.078125e-1 <= T) & (T < 5.15625e-1));
t50 = ((3.75e-1 <= T) & (T < 3.828125e-1));
t88 = ((6.71875e-1 <= T) & (T < 6.796875e-1));
t31 = ((2.265625e-1 <= T) & (T < 2.34375e-1));
t109 = ((8.359375e-1 <= T) & (T < 8.4375e-1));
t40 = ((T < 3.046875e-1) & (2.96875e-1 <= T));
t68 = ((5.15625e-1 <= T) & (T < 5.234375e-1));
t26 = ((1.875e-1 <= T) & (T < 1.953125e-1));
t89 = ((6.796875e-1 <= T) & (T < 6.875e-1));
t61 = ((4.609375e-1 <= T) & (T < 4.6875e-1));
t110 = ((8.4375e-1 <= T) & (T < 8.515625e-1));
t21 = ((T < 1.5625e-1) & (1.484375e-1 <= T));
t69 = ((5.234375e-1 <= T) & (T < 5.3125e-1));
if ~all(cellfun(@isscalar,{T,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12,t120,t121,t122,t123,t124,t125,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*3.376467202709134e-7-1.318592942028946e-9;
elseif (t3)
    out1 = T.*-3.381682821779969e-7+3.961211764603166e-9;
elseif (t4)
    out1 = T.*3.402545286745024e-7-6.639144654967137e-9;
elseif (t5)
    out1 = T.*-3.433516548516111e-7+9.382875271426148e-9;
elseif (t6)
    out1 = T.*3.442013172327823e-7-1.210315510621114e-8;
elseif (t7)
    out1 = T.*-3.425269976636164e-7+1.472216969442943e-8;
elseif (t8)
    out1 = T.*3.404809681460343e-7-1.729382870289795e-8;
elseif (t9)
    out1 = T.*-3.386974680223571e-7+1.984874202506096e-8;
elseif (t10)
    out1 = T.*3.375611826555269e-7-2.241742364230679e-8;
elseif (t11)
    out1 = T.*-3.37450152874709e-7+2.504431088716292e-8;
elseif (t12)
    out1 = T.*3.394918554915575e-7-2.784178351645165e-8;
elseif (t13)
    out1 = T.*-3.430433966507639e-7+3.08135897145291e-8;
elseif (t14)
    out1 = T.*3.444057253455118e-7-3.363476547262176e-8;
elseif (t15)
    out1 = T.*-3.456286261173925e-7+3.644684834782947e-8;
elseif (t16)
    out1 = T.*3.58610287975421e-7-4.057928288107201e-8;
elseif (t17)
    out1 = T.*-4.140403334461191e-7+4.996571181676472e-8;
elseif (t18)
    out1 = T.*6.227790631214732e-7-7.963671275418431e-8;
elseif (t19)
    out1 = T.*-9.177481206946375e-7+1.249645538463929e-7;
elseif (t20)
    out1 = T.*9.197737802211644e-7-1.334369634698918e-7;
elseif (t21)
    out1 = T.*-6.281810612039944e-7+9.633758330415525e-8;
elseif (t22)
    out1 = T.*4.194437686886625e-7-6.735379636657239e-8;
elseif (t23)
    out1 = T.*-3.606432402766742e-7+6.062922854180317e-8;
elseif (t24)
    out1 = T.*3.436310196300749e-7-6.041790987966934e-8;
elseif (t25)
    out1 = T.*-3.391089660058765e-7+6.226193128929067e-8;
elseif (t26)
    out1 = T.*3.380329741164798e-7-6.470218248365113e-8;
elseif (t27)
    out1 = T.*-3.382510591887433e-7+6.738454277127526e-8;
elseif (t28)
    out1 = T.*3.401993876215741e-7-7.042570423707045e-8;
elseif (t29)
    out1 = T.*-3.430483170946219e-7+7.369685847650213e-8;
elseif (t30)
    out1 = T.*3.430431101188018e-7-7.63856412264343e-8;
elseif (t31)
    out1 = T.*-3.401733519115369e-7+7.84055884523143e-8;
elseif (t32)
    out1 = T.*3.38152122642192e-7-8.05769446462159e-8;
elseif (t33)
    out1 = T.*-3.37663261498327e-7+8.309709370031604e-8;
elseif (t34)
    out1 = T.*3.377290448289791e-7-8.575098288151049e-8;
elseif (t35)
    out1 = T.*-3.384810388940414e-7+8.858442932833075e-8;
elseif (t36)
    out1 = T.*3.414231749925399e-7-9.201512748529242e-8;
elseif (t37)
    out1 = T.*-3.477133138198841e-7+9.642063117435478e-8;
elseif (t38)
    out1 = T.*3.604411256051232e-7-1.027478049139285e-7;
elseif (t39)
    out1 = T.*-4.030104808855104e-7+1.179374250872702e-7;
elseif (t40)
    out1 = T.*4.633551970853335e-7-1.392648855603241e-7;
elseif (t41)
    out1 = T.*-4.615553165791989e-7+1.425437865718382e-7;
elseif (t42)
    out1 = T.*3.980007906241923e-7-1.260674969292216e-7;
elseif (t43)
    out1 = T.*-3.536448962768632e-7+1.146940121562728e-7;
elseif (t44)
    out1 = T.*3.438312253624789e-7-1.141653402566364e-7;
elseif (t45)
    out1 = T.*-3.442061291635299e-7+1.169722085294447e-7;
elseif (t46)
    out1 = T.*3.447172801366036e-7-1.198452134174762e-7;
elseif (t47)
    out1 = T.*-3.457121080149498e-7+1.228838683545543e-7;
elseif (t48)
    out1 = T.*3.586326580044223e-7-1.302400319336575e-7;
elseif (t49)
    out1 = T.*-4.140463287329458e-7+1.534780335089698e-7;
elseif (t50)
    out1 = T.*6.227806710715932e-7-2.353320914177323e-7;
elseif (t51)
    out1 = T.*-9.177485531231193e-7+3.544017522193061e-7;
elseif (t52)
    out1 = T.*9.197738964202641e-7-3.633804546335781e-7;
elseif (t53)
    out1 = T.*-6.281810915340945e-7+2.533828608794867e-7;
elseif (t54)
    out1 = T.*4.19443776805495e-7-1.722147418834715e-7;
elseif (t55)
    out1 = T.*-3.606432438146524e-7+1.507900400920582e-7;
elseif (t56)
    out1 = T.*3.436310229381826e-7-1.46325666194294e-7;
elseif (t57)
    out1 = T.*-3.391089696949472e-7+1.470391743902539e-7;
elseif (t58)
    out1 = T.*3.380329756884634e-7-1.492104267149882e-7;
elseif (t59)
    out1 = T.*-3.38251053328863e-7+1.5194730495679e-7;
elseif (t60)
    out1 = T.*3.401993588575594e-7-1.554755380651827e-7;
elseif (t61)
    out1 = T.*-3.430482038437576e-7+1.594588853674556e-7;
elseif (t62)
    out1 = T.*3.430426801150576e-7-1.62146216488239e-7;
elseif (t63)
    out1 = T.*-3.401717398914997e-7+1.63448155546136e-7;
elseif (t64)
    out1 = T.*3.381461000468271e-7-1.651120481739911e-7;
elseif (t65)
    out1 = T.*-3.376407798649623e-7+1.675018067825928e-7;
elseif (t66)
    out1 = T.*3.376451389856375e-7-1.701411526427071e-7;
elseif (t67)
    out1 = T.*-3.381678955659523e-7+1.73045153965522e-7;
elseif (t68)
    out1 = T.*3.402545628944545e-7-1.767664261781252e-7;
elseif (t69)
    out1 = T.*-3.433521761071594e-7+1.810589762680321e-7;
elseif (t70)
    out1 = T.*3.44203362185731e-7-1.842049034500659e-7;
elseif (t71)
    out1 = T.*-3.425346511023786e-7+1.859898068380556e-7;
elseif (t72)
    out1 = T.*3.405095351360644e-7-1.875499825110929e-7;
elseif (t73)
    out1 = T.*-3.388040878517405e-7+1.892567927399551e-7;
elseif (t74)
    out1 = T.*3.379591130139479e-7-1.914225077469946e-7;
elseif (t75)
    out1 = T.*-3.382605039756594e-7+1.942339925673908e-7;
elseif (t76)
    out1 = T.*3.4031102917355e-7-1.980651750344959e-7;
elseif (t77)
    out1 = T.*-3.434854353812263e-7+2.025968159155684e-7;
elseif (t78)
    out1 = T.*3.446799389980613e-7-2.060013751221337e-7;
elseif (t79)
    out1 = T.*-3.462834379782787e-7+2.096562813401959e-7;
elseif (t80)
    out1 = T.*3.60949856843754e-7-2.213140076919803e-7;
elseif (t81)
    out1 = T.*-4.227274184371938e-7+2.623618106454797e-7;
elseif (t82)
    out1 = T.*6.404235592097156e-7-4.021075503838387e-7;
elseif (t83)
    out1 = T.*-9.353898644833334e-7+5.950868817969189e-7;
elseif (t84)
    out1 = T.*9.284525601500356e-7-5.989371714838331e-7;
elseif (t85)
    out1 = T.*-6.305065343739094e-7+4.119503663715375e-7;
elseif (t86)
    out1 = T.*4.200668788494622e-7-2.774884360563001e-7;
elseif (t87)
    out1 = T.*-3.608102036733227e-7+2.410627515564867e-7;
elseif (t88)
    out1 = T.*3.436757556554045e-7-2.322637523675019e-7;
elseif (t89)
    out1 = T.*-3.391209417575976e-7+2.31824627905398e-7;
elseif (t90)
    out1 = T.*3.380361423379532e-7-2.337208674102932e-7;
elseif (t91)
    out1 = T.*-3.38251758474935e-7+2.365105636236681e-7;
elseif (t92)
    out1 = T.*3.401990146046176e-7-2.405251361978923e-7;
elseif (t93)
    out1 = T.*-3.430461182220728e-7+2.452194504210829e-7;
elseif (t94)
    out1 = T.*3.430346749902736e-7-2.47901119700291e-7;
elseif (t95)
    out1 = T.*-3.401417986458204e-7+2.484692869259335e-7;
elseif (t96)
    out1 = T.*3.380343378785345e-7-2.495663133341397e-7;
elseif (t97)
    out1 = T.*-3.372236713086865e-7+2.516017403595009e-7;
elseif (t98)
    out1 = T.*3.36088465358911e-7-2.533823621411972e-7;
elseif (t99)
    out1 = T.*-3.323583071653259e-7+2.531749576623261e-7;
elseif (t100)
    out1 = T.*3.185728786459048e-7-2.451942314743975e-7;
elseif (t101)
    out1 = T.*-2.671583992340141e-7+2.078323037608523e-7;
elseif (t102)
    out1 = T.*7.528006434736409e-8-5.969774591209936e-8;
elseif (t103)
    out1 = T.*6.427946384114435e-7-5.07502214509537e-7;
elseif (t104)
    out1 = T.*-3.327160952971367e-6+2.65605614737364e-6;
elseif (t105)
    out1 = T.*1.334654858246577e-5-1.076106949442343e-5;
elseif (t106)
    out1 = T.*-5.073578111015568e-5+4.13058233808315e-5;
elseif (t107)
    out1 = T.*1.90271347739586e-4-1.563953370037222e-4;
elseif (t108)
    out1 = T.*-7.110243816140213e-4+5.899901888672338e-4;
elseif (t109)
    out1 = T.*2.654500953950863e-3-2.223378646331537e-3;
elseif (t110)
    out1 = T.*-9.907654214839125e-3+8.375939777335015e-3;
elseif (t111)
    out1 = T.*3.697679068006205e-2-3.154909532847926e-2;
elseif (t112)
    out1 = T.*-1.380001832784836e-1+1.188217416671459e-1;
elseif (t113)
    out1 = T.*5.150244484945662e-1-4.474730561985457e-1;
elseif (t114)
    out1 = T.*-1.922097779340827+1.685008893157424;
elseif (t115)
    out1 = T.*7.17336666925631-6.344580815369737;
elseif (t116)
    out1 = T.*-2.67713689000778e+1+2.388744930106846e+1;
elseif (t117)
    out1 = T.*9.991210892295774e+1-8.992973780556504e+1;
elseif (t118)
    out1 = T.*-3.728770667917532e+2+3.385354526858917e+2;
elseif (t119)
    out1 = T.*1.391596158356834e+3-1.274303354676488e+3;
elseif (t120)
    out1 = T.*-2.770707512436749e+3+2.562820341836346e+3;
elseif (t121)
    out1 = T.*2.422833728787135e+3-2.265550030863984e+3;
elseif (t122)
    out1 = T.*-4.752894380190612e+2+4.514404380168248e+2;
elseif (t123)
    out1 = T.*1.275166784935943e+2-1.231091417842999e+2;
elseif (t124)
    out1 = T.*-3.477727595272368e+1+3.411312658557064e+1;
elseif (t125)
    out1 = T.*1.159242531761673e+1-1.15320481024207e+1;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
