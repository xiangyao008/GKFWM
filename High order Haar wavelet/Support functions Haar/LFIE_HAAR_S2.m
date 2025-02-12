function [Array_fx,tk_xj] = LFIE_HAAR_S2(fx, kf, n)
% fx 为给定函数，接受x参数 定义为函数句柄
% kf 为已知积分核，接受x,t两个参数 定义为函数句柄
% n 为阶数
% OUTPUT_c1 为 常数
% Array_fx 为 1*n，表示ai系数
%% 预分配参数
haar_wavelet_integral= @(x,alpha,beta,gamma,s)0.*((0<=x) < alpha)+((x - alpha).^s./ factorial(s)).*(alpha<=x&x< beta)+...
(((x - alpha).^s - 2*(x - beta).^s) ./ factorial(s)).* (beta<=x& x< gamma)+ ...
(((x - alpha).^s - 2*(x - beta).^s + (x - gamma).^s) ./ factorial(s)).*(gamma<=x&x <= 1);
P_ns=zeros(n,n);
P_ns_x1=zeros(n,1);
j=ceil(log(n)/log(2));
s=2;
[alpha,beta,gamma]=generate_haar(j); 
tk_xj=CP(n);
%% 计算参数
s2=sum(kf(0,tk_xj)./n);%论文中的s2
s3=sum(kf(1,tk_xj)./n);
s4=sum(kf(0,tk_xj).*tk_xj./n);
s5=sum(kf(1,tk_xj).*tk_xj./n);
for i=1:n
    P_ns(i,:)=haar_wavelet_integral(tk_xj,alpha(i),beta(i),gamma(i),s); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
    P_ns_x1(i)=haar_wavelet_integral(1,alpha(i),beta(i),gamma(i),s);
end
P_ns=P_ns';
Constant_kf=zeros(n,n);
for rr=1:n
    for cc=1:n
        Constant_kf(rr,cc)=kf(tk_xj(rr),tk_xj(cc)); % 表示k(x,t)的矩阵，
    end
end
% s1=(((P_ns')*(kf(0,tk_xj)'))./n)';
s6=(((P_ns')*(kf(0,tk_xj)'))./n)';
s6_arr=repmat(s6,n,1);
s7=(((P_ns')*(kf(1,tk_xj)'))./n)';
s7_arr=repmat(s7,n,1);
P_ns_x1_arr=repmat(P_ns_x1',n,1);
s8=1-s2+s4*(1-s3)-s5*(1-s2);
%% 化简
% C2 表示公式13中最后的常数项
Constant_xj=(fx(tk_xj)'-1/s8*(fx(0)*(1-s5)+fx(1)*s4)*(ones(n,1)-(Constant_kf*ones(n,1))./n))...
            -1/s8*(-fx(0)*(1-s3)+fx(1)*(1-s2))*(tk_xj'-(Constant_kf*(tk_xj'))./n);
Array_ax=(P_ns-(Constant_kf*P_ns)./n)+1/s8*((s6_arr.*(1-s5)-(s4.*(P_ns_x1_arr-s7_arr))).*(ones(n,1)-(Constant_kf*ones(n,1))./n))...
        +1/s8*(((-(1-s3).*s6_arr)-(1-s2).*(P_ns_x1_arr-s7_arr)).*(tk_xj'-(Constant_kf*(tk_xj'))./n));
%% 主函数
Array_fx=pinv(Array_ax)*Constant_xj;
OUTPUT_c2=1/s8*(-fx(0)*(1-s3)+fx(1)*(1-s2)-(1-s3)*(s6*Array_fx)-(1-s2)*(P_ns_x1'-s7)*Array_fx);
OUTPUT_c1=1/s8*(fx(0)*(1-s5)+fx(1)*s4+(1-s5)*(s6*Array_fx)-s4*(P_ns_x1'-s7)*Array_fx);
Array_fx=[OUTPUT_c1;Array_fx;OUTPUT_c2];
end

