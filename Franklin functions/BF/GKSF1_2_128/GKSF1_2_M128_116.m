function out1 = GKSF1_2_M128_116(T)
%GKSF1_2_M128_116
%    OUT1 = GKSF1_2_M128_116(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:42:51

t51 = ((3.828125e-1 <= T) & (T < 3.90625e-1));
t90 = ((6.875e-1 <= T) & (T < 6.953125e-1));
t14 = ((9.375e-2 <= T) & (T < 1.015625e-1));
t41 = ((3.046875e-1 <= T) & (T < 3.125e-1));
t70 = ((5.3125e-1 <= T) & (T < 5.390625e-1));
t91 = ((6.953125e-1 <= T) & (T < 7.03125e-1));
t62 = ((4.6875e-1 <= T) & (T < 4.765625e-1));
t4 = ((T < 2.34375e-2) & (1.5625e-2 <= T));
t5 = ((2.34375e-2 <= T) & (T < 3.125e-2));
t113 = ((9.375e-1 <= T) & (T < 9.53125e-1));
t71 = ((5.390625e-1 <= T) & (T < 5.46875e-1));
t52 = ((3.90625e-1 <= T) & (T < 3.984375e-1));
t92 = ((7.03125e-1 <= T) & (T < 7.109375e-1));
t32 = ((2.34375e-1 <= T) & (T < 2.421875e-1));
t42 = ((3.125e-1 <= T) & (T < 3.203125e-1));
t72 = ((5.46875e-1 <= T) & (T < 5.546875e-1));
t27 = ((1.953125e-1 <= T) & (T < 2.03125e-1));
t93 = ((7.109375e-1 <= T) & (T < 7.1875e-1));
t63 = ((4.765625e-1 <= T) & (T < 4.84375e-1));
t114 = ((9.53125e-1 <= T) & (T < 9.6875e-1));
t22 = ((1.5625e-1 <= T) & (T < 1.640625e-1));
t73 = ((5.546875e-1 <= T) & (T < 5.625e-1));
t53 = ((3.984375e-1 <= T) & (T < 4.0625e-1));
t94 = ((7.1875e-1 <= T) & (T < 7.265625e-1));
t104 = ((7.96875e-1 <= T) & (T < 8.125e-1));
t8 = ((4.6875e-2 <= T) & (T < 5.46875e-2));
t43 = ((3.203125e-1 <= T) & (T < 3.28125e-1));
t74 = ((5.625e-1 <= T) & (T < 5.703125e-1));
t95 = ((7.265625e-1 <= T) & (T < 7.34375e-1));
t64 = ((4.84375e-1 <= T) & (T < 4.921875e-1));
t115 = ((9.6875e-1 <= T) & (T < 9.84375e-1));
t12 = ((7.8125e-2 <= T) & (T < 8.59375e-2));
t75 = ((5.703125e-1 <= T) & (T < 5.78125e-1));
t3 = ((T < 1.5625e-2) & (7.8125e-3 <= T));
t54 = ((4.0625e-1 <= T) & (T < 4.140625e-1));
t96 = ((7.34375e-1 <= T) & (T < 7.421875e-1));
t105 = ((8.125e-1 <= T) & (T < 8.28125e-1));
t44 = ((3.28125e-1 <= T) & (T < 3.359375e-1));
t76 = ((5.78125e-1 <= T) & (T < 5.859375e-1));
t28 = ((2.03125e-1 <= T) & (T < 2.109375e-1));
t97 = ((7.421875e-1 <= T) & (T < 7.5e-1));
t34 = ((2.5e-1 <= T) & (T < 2.578125e-1));
t23 = ((1.640625e-1 <= T) & (T < 1.71875e-1));
t77 = ((5.859375e-1 <= T) & (T < 5.9375e-1));
t55 = ((4.140625e-1 <= T) & (T < 4.21875e-1));
t98 = ((7.5e-1 <= T) & (T < 7.578125e-1));
t106 = ((8.28125e-1 <= T) & (T < 8.4375e-1));
t18 = ((1.25e-1 <= T) & (T < 1.328125e-1));
t15 = ((1.015625e-1 <= T) & (T < 1.09375e-1));
t45 = ((3.359375e-1 <= T) & (T < 3.4375e-1));
t78 = ((T < 6.015625e-1) & (5.9375e-1 <= T));
t99 = ((7.578125e-1 <= T) & (T < 7.65625e-1));
t7 = ((3.90625e-2 <= T) & (T < 4.6875e-2));
t9 = ((5.46875e-2 <= T) & (T < 6.25e-2));
t35 = ((2.578125e-1 <= T) & (T < 2.65625e-1));
t79 = ((6.015625e-1 <= T) & (T < 6.09375e-1));
t56 = ((4.21875e-1 <= T) & (T < 4.296875e-1));
t100 = ((7.65625e-1 <= T) & (T < 7.734375e-1));
t107 = ((8.4375e-1 <= T) & (T < 8.59375e-1));
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
t108 = ((8.59375e-1 <= T) & (T < 8.75e-1));
t19 = ((1.328125e-1 <= T) & (T < 1.40625e-1));
t47 = ((3.515625e-1 <= T) & (T < 3.59375e-1));
t82 = ((6.25e-1 <= T) & (T < 6.328125e-1));
t103 = ((7.890625e-1 <= T) & (T < 7.96875e-1));
t13 = ((8.59375e-2 <= T) & (T < 9.375e-2));
t37 = ((2.734375e-1 <= T) & (T < 2.8125e-1));
t17 = ((1.171875e-1 <= T) & (T < 1.25e-1));
t83 = ((6.328125e-1 <= T) & (T < 6.40625e-1));
t58 = ((4.375e-1 <= T) & (T < 4.453125e-1));
t109 = ((8.75e-1 <= T) & (T < 8.90625e-1));
t6 = ((T < 3.90625e-2) & (3.125e-2 <= T));
t48 = ((3.59375e-1 <= T) & (T < 3.671875e-1));
t84 = ((6.40625e-1 <= T) & (T < 6.484375e-1));
t30 = ((2.1875e-1 <= T) & (T < 2.265625e-1));
t38 = ((2.8125e-1 <= T) & (T < 2.890625e-1));
t25 = ((1.796875e-1 <= T) & (T < 1.875e-1));
t33 = ((2.421875e-1 <= T) & (T < 2.5e-1));
t85 = ((6.484375e-1 <= T) & (T < 6.5625e-1));
t59 = ((4.453125e-1 <= T) & (T < 4.53125e-1));
t110 = ((8.90625e-1 <= T) & (T < 9.0625e-1));
t20 = ((1.40625e-1 <= T) & (T < 1.484375e-1));
t16 = ((1.09375e-1 <= T) & (T < 1.171875e-1));
t49 = ((3.671875e-1 <= T) & (T < 3.75e-1));
t86 = ((6.5625e-1 <= T) & (T < 6.640625e-1));
t65 = ((4.921875e-1 <= T) & (T < 5.0e-1));
t39 = ((2.890625e-1 <= T) & (T < 2.96875e-1));
t66 = ((5.0e-1 <= T) & (T < 5.078125e-1));
t116 = ((9.84375e-1 <= T) & (T <= 1.0));
t87 = ((6.640625e-1 <= T) & (T < 6.71875e-1));
t60 = ((4.53125e-1 <= T) & (T < 4.609375e-1));
t111 = ((9.0625e-1 <= T) & (T < 9.21875e-1));
t2 = ((0.0 <= T) & (T < 7.8125e-3));
t11 = ((T < 7.8125e-2) & (7.03125e-2 <= T));
t67 = ((5.078125e-1 <= T) & (T < 5.15625e-1));
t50 = ((3.75e-1 <= T) & (T < 3.828125e-1));
t88 = ((6.71875e-1 <= T) & (T < 6.796875e-1));
t31 = ((2.265625e-1 <= T) & (T < 2.34375e-1));
t40 = ((T < 3.046875e-1) & (2.96875e-1 <= T));
t68 = ((5.15625e-1 <= T) & (T < 5.234375e-1));
t26 = ((1.875e-1 <= T) & (T < 1.953125e-1));
t89 = ((6.796875e-1 <= T) & (T < 6.875e-1));
t61 = ((4.609375e-1 <= T) & (T < 4.6875e-1));
t112 = ((9.21875e-1 <= T) & (T < 9.375e-1));
t21 = ((T < 1.5625e-1) & (1.484375e-1 <= T));
t69 = ((5.234375e-1 <= T) & (T < 5.3125e-1));
if ~all(cellfun(@isscalar,{T,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t11,t110,t111,t112,t113,t114,t115,t116,t12,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*5.288747876100288e-7-1.377275936984333e-9;
elseif (t3)
    out1 = T.*-9.118728112020791e-7+9.87856467873526e-9;
elseif (t4)
    out1 = T.*1.094355180408897e-6-2.146874769018624e-8;
elseif (t5)
    out1 = T.*-9.668031144936394e-7+2.683964984659196e-8;
elseif (t6)
    out1 = T.*6.746551665009901e-7-2.445592143449022e-8;
elseif (t7)
    out1 = T.*-3.083218558800193e-7+1.394161850226796e-8;
elseif (t8)
    out1 = T.*8.400670978937722e-8-4.448783013485e-9;
elseif (t9)
    out1 = T.*-2.756295619633141e-8+1.652683095108441e-9;
elseif (t10)
    out1 = T.*2.619997921628107e-8-1.707500368179839e-9;
elseif (t11)
    out1 = T.*-7.695654753796758e-8+5.545692919228269e-9;
elseif (t12)
    out1 = T.*2.807783506223251e-7-2.240234599954459e-8;
elseif (t13)
    out1 = T.*-5.708254743765457e-7+5.078235771129586e-8;
elseif (t14)
    out1 = T.*5.787925735500679e-7-5.699433428182416e-8;
elseif (t15)
    out1 = T.*-3.021010478378168e-7+3.247142414038288e-8;
elseif (t16)
    out1 = T.*9.874184255729838e-8-1.137076699658285e-8;
elseif (t17)
    out1 = T.*-3.648006793319774e-8+4.475550639022169e-9;
elseif (t18)
    out1 = T.*2.838521767130815e-8-3.632610061541067e-9;
elseif (t19)
    out1 = T.*-7.705715590992387e-8+1.037145517971631e-8;
elseif (t20)
    out1 = T.*2.798258463975475e-7-3.981521701977185e-8;
elseif (t21)
    out1 = T.*-5.677451489450556e-7+8.59961026013958e-8;
elseif (t22)
    out1 = T.*5.677005395261424e-7-9.141728622222889e-8;
elseif (t23)
    out1 = T.*-2.796097579029456e-7+4.759455944973085e-8;
elseif (t24)
    out1 = T.*7.62582861630624e-8-1.357026062411427e-8;
elseif (t25)
    out1 = T.*-2.54267034610669e-8+4.701260948971468e-9;
elseif (t26)
    out1 = T.*2.545184466526928e-8-4.838466824716568e-9;
elseif (t27)
    out1 = T.*-7.638399052429262e-8+1.505134473574474e-8;
elseif (t28)
    out1 = T.*2.800874291778671e-7-5.735691239125645e-8;
elseif (t29)
    out1 = T.*-5.69485517053819e-7+1.218498809544898e-7;
elseif (t30)
    out1 = T.*5.744073351821074e-7-1.283766804721191e-7;
elseif (t31)
    out1 = T.*-3.046934819097869e-7+7.079459840026325e-8;
elseif (t32)
    out1 = T.*1.697253946985297e-7-4.039732580481096e-8;
elseif (t33)
    out1 = T.*-3.737464416837132e-7+9.122475956901348e-8;
elseif (t34)
    out1 = T.*9.139802112429319e-7-2.307069036626478e-7;
elseif (t35)
    out1 = T.*-1.374797572227136e-6+3.59368618638229e-7;
elseif (t36)
    out1 = T.*1.328475566422861e-6-3.586883088156764e-7;
elseif (t37)
    out1 = T.*-1.029535483632707e-6+2.86080337683893e-7;
elseif (t38)
    out1 = T.*6.91464254438894e-7-1.979508386487448e-7;
elseif (t39)
    out1 = T.*-3.12825837036071e-7+9.23517659182373e-8;
elseif (t40)
    out1 = T.*8.521354434742156e-8-2.581617542998705e-8;
elseif (t41)
    out1 = T.*-2.788631908934779e-8+8.643939210903613e-9;
elseif (t42)
    out1 = T.*2.62866181605083e-8-8.285103679676416e-9;
elseif (t43)
    out1 = T.*-7.697976145091634e-8+2.479240853960804e-8;
elseif (t44)
    out1 = T.*2.807845724876935e-7-9.259901353399833e-8;
elseif (t45)
    out1 = T.*-5.708271438593063e-7+1.934892974263219e-7;
elseif (t46)
    out1 = T.*5.787930225058185e-7-2.016926347616897e-7;
elseif (t47)
    out1 = T.*-3.021011686576694e-7+1.07996729319224e-7;
elseif (t48)
    out1 = T.*9.874187533208771e-8-3.605623961459497e-8;
elseif (t49)
    out1 = T.*-3.648007709712225e-8+1.359557104300557e-8;
elseif (t50)
    out1 = T.*2.838521994861936e-8-1.072891534914754e-8;
elseif (t51)
    out1 = T.*-7.705715440501271e-8+2.963574358310224e-8;
elseif (t52)
    out1 = T.*2.798258364955371e-7-1.09771674737425e-7;
elseif (t53)
    out1 = T.*-5.677451128349202e-7+2.27932375386429e-7;
elseif (t54)
    out1 = T.*5.677004131937906e-7-2.333423695627347e-7;
elseif (t55)
    out1 = T.*-2.796092949311529e-7+1.174968064577497e-7;
elseif (t56)
    out1 = T.*7.625655863282371e-8-3.263410051455293e-8;
elseif (t57)
    out1 = T.*-2.542024974274373e-8+1.105515308432371e-8;
elseif (t58)
    out1 = T.*2.542774933955573e-8-1.119084651418231e-8;
elseif (t59)
    out1 = T.*-7.629405671898175e-8+3.410714524626017e-8;
elseif (t60)
    out1 = T.*2.797517871623505e-7-1.272261277624685e-7;
elseif (t61)
    out1 = T.*-5.682328811913572e-7+2.636418053068187e-7;
elseif (t62)
    out1 = T.*5.697324865176492e-7-2.697794608067781e-7;
elseif (t63)
    out1 = T.*-2.872468885425746e-7+1.386247726203598e-7;
elseif (t64)
    out1 = T.*1.047660811692399e-7-5.125650958380038e-8;
elseif (t65)
    out1 = T.*-1.318119699890564e-7+6.51842499706736e-8;
elseif (t66)
    out1 = T.*4.224821835427509e-7-2.1196282679523e-7;
elseif (t67)
    out1 = T.*-8.833652134843536e-7+4.511628045075964e-7;
elseif (t68)
    out1 = T.*1.086717396420583e-6-5.646610412246364e-7;
elseif (t69)
    out1 = T.*-9.647595749069522e-7+5.091589359546202e-7;
elseif (t70)
    out1 = T.*6.741187838698836e-7-3.614951921455738e-7;
elseif (t71)
    out1 = T.*-3.082198598042016e-7+1.680467329599877e-7;
elseif (t72)
    out1 = T.*8.413511134175412e-8-4.652239188545678e-8;
elseif (t73)
    out1 = T.*-2.817856588662521e-8+1.577660095215988e-8;
elseif (t74)
    out1 = T.*2.852857350548946e-8-1.612116495590462e-8;
elseif (t75)
    out1 = T.*-8.563896286745499e-8+4.899000813179027e-8;
elseif (t76)
    out1 = T.*2.984151313675661e-7-1.730412650978313e-7;
elseif (t77)
    out1 = T.*-5.884409420175272e-7+3.466009654012468e-7;
elseif (t78)
    out1 = T.*5.873738774616201e-7-3.515390836644969e-7;
elseif (t79)
    out1 = T.*-3.040628712720914e-7+1.847158354956264e-7;
elseif (t80)
    out1 = T.*9.80078355312325e-8-6.029600146264909e-8;
elseif (t81)
    out1 = T.*-3.158221669034612e-8+1.968535889285647e-8;
elseif (t82)
    out1 = T.*9.527819667962902e-9-6.008413831086669e-9;
elseif (t83)
    out1 = T.*-6.525414740455074e-9+4.150273567990331e-9;
elseif (t84)
    out1 = T.*1.655627879067862e-8-1.063643635039219e-8;
elseif (t85)
    out1 = T.*4.148013804173914e-7-2.688734944364638e-7;
elseif (t86)
    out1 = T.*-3.099216008085018e-6+2.037200416768242e-6;
elseif (t87)
    out1 = T.*1.340550990109415e-5-8.922969132296048e-6;
elseif (t88)
    out1 = T.*-5.099730381615392e-5+3.434767133398e-5;
elseif (t89)
    out1 = T.*1.90583702059295e-4-1.298519185969892e-4;
elseif (t90)
    out1 = T.*-7.113375010105561e-4+4.902189085135335e-4;
elseif (t91)
    out1 = T.*2.654766299013283e-3-1.850275139940542e-3;
elseif (t92)
    out1 = T.*-9.907727691423366e-3+6.982728447085227e-3;
elseif (t93)
    out1 = T.*3.697661894048012e-2-2.634911173653366e-2;
elseif (t94)
    out1 = T.*-1.380001715397028e-1+9.941545642109779e-2;
elseif (t95)
    out1 = T.*5.150254906651683e-1-3.750485012746289e-1;
elseif (t96)
    out1 = T.*-1.922102265289218+1.414717194504374;
elseif (t97)
    out1 = T.*7.173383569980496-5.335838698859867;
elseif (t98)
    out1 = T.*-2.677143201752188e+1+2.012277299176691e+1;
elseif (t99)
    out1 = T.*9.991234516946004e+1-7.587977690774282e+1;
elseif (t100)
    out1 = T.*-3.728779507512505e+2+2.861002934065512e+2;
elseif (t101)
    out1 = T.*1.391599460334307e+3-1.078612704229934e+3;
elseif (t102)
    out1 = T.*-2.770709878379645e+3+2.17319146664034e+3;
elseif (t103)
    out1 = T.*2.42281001324035e+3-1.924820322841062e+3;
elseif (t104)
    out1 = T.*-4.752400538566615e+2+3.845633243768691e+2;
elseif (t105)
    out1 = T.*1.273401887058032e+2-1.050331227051336e+2;
elseif (t106)
    out1 = T.*-3.412070085080371e+1+2.867667645893157e+1;
elseif (t107)
    out1 = T.*9.142614340336268-7.826745733592788;
elseif (t108)
    out1 = T.*-2.449756150850269+2.135447657270642;
elseif (t109)
    out1 = T.*6.564101640042456e-1-5.824478682270581e-1;
elseif (t110)
    out1 = T.*-1.758845684331377e-1+1.588146278499865e-1;
elseif (t111)
    out1 = T.*4.712817298947583e-2-4.329066906425705e-2;
elseif (t112)
    out1 = T.*-1.262809297389444e-2+1.179713862072492e-2;
elseif (t113)
    out1 = T.*3.384043995697459e-3-3.214239788267486e-3;
elseif (t114)
    out1 = T.*-9.079280972391901e-4+8.76546112812757e-4;
elseif (t115)
    out1 = T.*2.476167543706368e-4-2.428879621842628e-4;
elseif (t116)
    out1 = T.*-8.253891833602334e-5+8.210902813635577e-5;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
