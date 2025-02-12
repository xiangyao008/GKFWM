function [output] = CP_CWM(n)
% CP 生成 [0,1] 区间上的 Chebyshev 配点
% n: 配点的数量
% output: 1×n 的数组，包含 Chebyshev 配点

% cheb_nodes = cos((2*(1:n) - 1) * pi / (2*n)); 
% output = (cheb_nodes + 1) / 2;

% output = 1/2+1/2.*(cos(pi*(1:n)./n)); 

output = 1/2-1/2.*(cos(pi*(1:n)./n)); 
end
