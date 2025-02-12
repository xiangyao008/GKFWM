function [Array_ax,tk_xj] = LFIE_HAAR_S1(fx, kf, n)
% fx 为给定函数，接受x参数 定义为函数句柄
% kf 为已知积分核，接受x,t两个参数 定义为函数句柄
% n 为阶数
% OUTPUT_c1 为 常数，对应论文中的C1
% Array_fx 为 1*n，表示ai系数
%% 预分配参数
haar_wavelet_integral= @(x,alpha,beta,gamma,s)0.*((0<=x) < alpha)+((x - alpha).^s./ factorial(s)).*(alpha<=x&x< beta)+...
(((x - alpha).^s - 2*(x - beta).^s) ./ factorial(s)).* (beta<=x& x< gamma)+ ...
(((x - alpha).^s - 2*(x - beta).^s + (x - gamma).^s) ./ factorial(s)).*(gamma<=x&x <= 1);
P_ns=zeros(n,n);
j=ceil(log(n)/log(2));
s=1;
[alpha,beta,gamma]=generate_haar(j); 
tk_xj=CP(n);
s2=sum(kf(0,tk_xj)./n);%论文中的s2
for i=1:n
    P_ns(i,:)=haar_wavelet_integral(tk_xj,alpha(i),beta(i),gamma(i),s); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
end
P_ns=P_ns';
Constant_kf=zeros(n,n);
for rr=1:n
    for cc=1:n
        Constant_kf(rr,cc)=kf(tk_xj(rr),tk_xj(cc)); % 表示k(x,t)的矩阵，
    end
end
% C2 表示公式13中最后的常数项
s1=(((P_ns')*(kf(0,tk_xj)'))./n)';
s1_new=repmat(s1,n,1); % s1_new为n*n,其中行向量表示s1(1),s1(2)...s1(n)
C1=1/(1-s2)*fx(0);
% C1_new=repmat(C1,n,1);
C2=((ones(n,1)-(Constant_kf*ones(n,1))./n).*s1_new).*(1/(1-s2));
% Array_fx=zeros(1,n);
%% 主函数
Array_ax=pinv(P_ns-((Constant_kf*P_ns)./n)+C2)*((fx(tk_xj)')-C1.*(ones(n,1)-(Constant_kf*ones(n,1))./n));
OUTPUT_c1=1/(1-s2)*(s1*Array_ax+fx(0));
Array_ax=[OUTPUT_c1;Array_ax];
end

