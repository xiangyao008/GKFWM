function [ua_iter,tk_xj] = LFIE_CWM(fx, kf, n)
% linear Fredholm integral equation
% fx 为给定函数, 接受x参数 定义为函数句柄
% kf 为函数句柄，表示核函数，接受x,t两个参数
% n 为阶数，iter,tol表示迭代此时和容差
% ux_iter 为 每次迭代后的值
% cg_it表示收敛时的阶次
%% 预分配参数
chebyshev_basis = @(n, t) cos(n * acos(2*t-1)); % Chebyshev多项式T_n(t)
tk_xj=CP_CWM(n);
tk_xj_new=[0,tk_xj]; %0开头
length_tk_xj=tk_xj-tk_xj_new(:,1:length(tk_xj));
% tk_xj=fliplr(tk_xj);
CW_ns=zeros(n,n);
kf_xt=zeros(n,n);
for i=1:n
    CW_ns(i,:)=chebyshev_basis(i-1,tk_xj); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
end
CW_ns=CW_ns';
for i=1:n
    for j=1:n
        kf_xt(i,j)=kf(tk_xj(i),tk_xj(j));
    end
end
fx_x=fx(tk_xj)';
%% 主函数
ua_iter=pinv(CW_ns-kf_xt*CW_ns./length(tk_xj))*fx_x;
end
