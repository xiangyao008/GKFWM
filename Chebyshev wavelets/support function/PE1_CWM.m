function [output] = PE1_CWM(RF,input_coefficient,n)
% 用来描述函数误差
% RF表示理论函数, 输入类型为函数句柄
% input_coefficient表示近似解的系数，输入为1*n大小的矩阵
% 输出output为一个实数，表示绝对误差
nn=CP_CWM(n);
nn=[nn,0.5];
% nn=fliplr(nn);
chebyshev_basis = @(n, t) cos(n * acos(2*t-1)); % Chebyshev多项式T_n(t)
% chebyshev_basis = @(n, t) 2*cos(n * acos(2*t-1))./sqrt(sqrt(1-(2*t-1).^2)); % Chebyshev多项式T_n(t)\
P_ns=zeros(length(input_coefficient),length(nn));
for i=1:length(input_coefficient)
    P_ns(i,:)=chebyshev_basis(i-1,nn);% P_ns的每一行表示n个配置点的i阶s次Chebyshev的数值
end
P_ns=P_ns';
output=abs(RF(nn)'-P_ns*input_coefficient);
end

