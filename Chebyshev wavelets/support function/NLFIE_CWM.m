function [ua_iter,tk_xj,cg_it] = NLFIE_CWM(fx,dkut, n, iter,tol)
% nonlinear Fredholm integral equation (chebyshev wavelet method)
% fx 为给定函数, 接受x参数 定义为函数句柄
% get_kf 和 get_dkf 表示卷积核和卷积核对被解函数u(t)的偏导数
% n 为阶数，iter,tol表示迭代此时和容差
% ux_iter 为 每次迭代后的值
% tk_xj 表示[0,1]区间的配置点
% cg_it表示收敛时的阶次
%% 预分配参数
chebyshev_basis = @(n, t) cos(n * acos(2*t-1)); % Chebyshev多项式T_n(t)\
% chebyshev_basis = @(n, t) 2*cos(n * acos(2*t-1))./sqrt(sqrt(1-(2*t-1).^2)); % Chebyshev多项式T_n(t)\
tk_xj=CP_CWM(n);
% tk_xj=fliplr(tk_xj);
% tk_xj_new=[tk_xj,1];
% length_tk_xj=tk_xj_new(:,2:length(tk_xj_new))-tk_xj;
tk_xj_new=[0,tk_xj]; %0开头
length_tk_xj=tk_xj-tk_xj_new(:,1:length(tk_xj));
%% 11
BdJx=cell(1,n);
Init_ux = zeros(n, 1);
Init_ux (1,1)=0.5;
% Init_ux (3,1)=-0.5;
ua_iter=zeros(n,iter+1);
ua_iter(:,1)=Init_ux;
CW_ns=zeros(n,n);
for i=1:n
    CW_ns(i,:)=chebyshev_basis(i-1,tk_xj); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
end
CW_ns=CW_ns';
u1=(CW_ns*ua_iter(:,1));
% kf1=get_dkf(tk_xj,tk_xj,u1);
% temp=repmat(kf1, 1, n);
% BdJx{1}=pinv(CW_ns-temp.*CW_ns);
dkut_xt=zeros(n,n);
for i=1:length(u1)
    for j=1:(length(u1))
        dkut_xt(i,j)=dkut(tk_xj(i),tk_xj(j),u1(j));
    end
end
% dkut_xt(1:length(tk_xj),1:length(tk_xj))=dkut_xt(1:length(tk_xj),1:length(tk_xj)).*tril(ones(length(tk_xj)),-1);
BdJx{1}=pinv(CW_ns-dkut_xt*CW_ns.*length_tk_xj);
% epsilon = 1e-6;
% J_approx = (get_kf(tk_xj, tk_xj, u1 + epsilon) - get_kf(tk_xj, tk_xj, u1)) / epsilon;
% BdJx{1} = pinv(J_approx);
%% 主函数
for i=2:iter+1  
    u=(CW_ns*ua_iter(:,i-1));% (n+1*1)
    S2_x=u;
    kf=get_kf(tk_xj,tk_xj,u);%(n*1)
    g_x=S2_x-fx(tk_xj)'-kf;%(n*1)
    delta_a=BdJx{i-1}*(-g_x);%(n*1)
    ua_iter(:,i)=ua_iter(:,i-1)+delta_a;%(n*1)
    u_new=(CW_ns*ua_iter(:,i));%(n*1)
    kf_new=get_kf(tk_xj,tk_xj,u_new);
    gx_new=u_new-fx(tk_xj)'-kf_new;
    delta_gx=gx_new-g_x;%(n+1*1)
    BdJx{i}=BdJx{i-1}+((delta_a-BdJx{i-1}*delta_gx)*(delta_a')*BdJx{i-1})/((delta_a')*BdJx{i-1}*delta_gx);
    if norm(ua_iter(:,i)-ua_iter(:,i-1))<tol
        disp(['Converged in ', num2str(i), ' iterations.']);
        cg_it=i;
        break;
    end
end
end

function kf=get_kf(x,t,u)
% x_new=[x,1];
% length_int=x_new(:,2:length(x_new))-x;
x_new=[0,x];
length_int=x-x_new(:,1:length(x));
temp=ones(1,length(x));
kt=sin(pi.*t').*(u.^3).*length_int';
kf=temp*kt.*1/5.*cos(pi.*x');
end
function kf=get_dkf(x,t,u)
kf=3/5.*cos(pi.*x').*sum(sin(pi.*t').* (u.^2))./(length(x)-1);
% t_new=t.*u;
%     [X, T] = ndgrid(x, t_new);
%     Kernel = (1/5 .* cos(pi*X) .* sin(pi*T) .* (T.^3))./length(x);
%     kf=sum(Kernel, 2);
end
