function out1 = GKSF2_2_M128_20(T)
%GKSF2_2_M128_20
%    OUT1 = GKSF2_2_M128_20(T)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    2024-08-02 17:45:36

t16 = ((7.5e-1 <= T) & (T < 8.125e-1));
t4 = ((T < 9.375e-2) & (6.25e-2 <= T));
t5 = ((9.375e-2 <= T) & (T < 1.25e-1));
t6 = ((T < 1.875e-1) & (1.25e-1 <= T));
t7 = ((1.875e-1 <= T) & (T < 2.5e-1));
t10 = ((3.75e-1 <= T) & (T < 4.375e-1));
t14 = ((6.25e-1 <= T) & (T < 6.875e-1));
t3 = ((T < 6.25e-2) & (3.125e-2 <= T));
t17 = ((8.125e-1 <= T) & (T < 8.75e-1));
t9 = ((3.125e-1 <= T) & (T < 3.75e-1));
t11 = ((4.375e-1 <= T) & (T < 5.0e-1));
t12 = ((5.0e-1 <= T) & (T < 5.625e-1));
t15 = ((6.875e-1 <= T) & (T < 7.5e-1));
t8 = ((T < 3.125e-1) & (2.5e-1 <= T));
t19 = ((9.375e-1 <= T) & (T <= 1.0));
t18 = ((8.75e-1 <= T) & (T < 9.375e-1));
t2 = ((0.0 <= T) & (T < 3.125e-2));
t13 = ((T < 6.25e-1) & (5.625e-1 <= T));
t20 = T.^2;
if ~all(cellfun(@isscalar,{T,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2,t3,t4,t5,t6,t7,t8,t9}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (T < 0.0)
    out1 = 0.0;
elseif (t2)
    out1 = T.*3.848814643419858e+2-t20.*1.124285541579004e+4-1.908951306645231;
elseif (t3)
    out1 = T.*-1.084142887545522e+3+t20.*1.226153421441009e+4+2.104455419159709e+1;
elseif (t4)
    out1 = T.*2.25629145085557e+3-t20.*1.446194049279865e+4-8.334401888343705e+1;
elseif (t5)
    out1 = T.*-2.421157545558986e+3+t20.*1.048445415474565e+4+1.359114028234953e+2;
elseif (t6)
    out1 = T.*7.516943737401432e+2-t20.*2.206953522450865e+3-6.239184213270032e+1;
elseif (t7)
    out1 = T.*-4.003953070967491e+2+t20.*8.652856264475145e+2+4.561656544575834e+1;
elseif (t8)
    out1 = T.*2.167018083429983e+2-t20.*3.689086044319804e+2-3.152057398421009e+1;
elseif (t9)
    out1 = T.*-1.130442833249574e+2+t20.*1.586851422367488e+2+2.000225283890799e+1;
elseif (t10)
    out1 = T.*5.720891719434208e+1-t20.*6.831912512231727e+1-1.192022225846067e+1;
elseif (t11)
    out1 = T.*-2.830955335374018e+1+t20.*2.941626978977675e+1+6.786943173932325;
elseif (t12)
    out1 = T.*1.377269071491511e+1-t20.*1.266597427887854e+1-3.733617843231499;
elseif (t13)
    out1 = T.*-6.612055354601159+t20.*5.453800005135918+1.999591988819952;
elseif (t14)
    out1 = T.*3.140942490889359-t20.*2.348598271256496-1.048219837895835;
elseif (t15)
    out1 = T.*-1.479872788448445+t20.*1.011994659170997+5.401854143765352e-1;
elseif (t16)
    out1 = T.*6.943011504995027e-1-t20.*4.374546334609679e-1-2.751298127289453e-1;
elseif (t17)
    out1 = T.*-3.291116821046276e-1+t20.*1.923378789108046e-1+1.406316505164826e-1;
elseif (t18)
    out1 = T.*1.692560770676378e-1-t20.*9.244369775906136e-2-7.740424412138354e-2;
elseif (t19)
    out1 = T.*-1.377322714575625e-1+t20.*7.128342145437881e-2+6.649654424980412e-2;
elseif (1.0 < T)
    out1 = 0.0;
else
    out1 = NaN;
end
end
