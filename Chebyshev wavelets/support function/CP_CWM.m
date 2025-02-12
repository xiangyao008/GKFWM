function [output] = CP_CWM(n)
% CP ���� [0,1] �����ϵ� Chebyshev ���
% n: ��������
% output: 1��n �����飬���� Chebyshev ���

% cheb_nodes = cos((2*(1:n) - 1) * pi / (2*n)); 
% output = (cheb_nodes + 1) / 2;

% output = 1/2+1/2.*(cos(pi*(1:n)./n)); 

output = 1/2-1/2.*(cos(pi*(1:n)./n)); 
end
