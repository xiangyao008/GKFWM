function [ua_iter,tk_xj,cg_it] = NLFVIE_CWM(fx,dvkut,dfkut, n, iter,tol)
% nonlinear Fredholm integral equation
% fx 为给定函数, 接受x参数 定义为函数句柄
% kuf 为已知积分核, 接受x,t,u三个参数 定义为函数句柄,包含被解函数的函数
% dkut表示kuf对被解函数u(t)的偏导
% n 为阶数，iter,tol表示迭代此时和容差
% ux_iter 为 每次迭代后的值
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
%% 预计算参数
BdJx=cell(1,n);
Init_ux = zeros(n, 1);
Init_ux (1,1)=-1;
% Init_ux (3,1)=-0.5;
ua_iter=zeros(n,iter+1);
ua_iter(:,1)=Init_ux;
CW_ns=zeros(n,n+1);
for i=1:n
    CW_ns(i,:)=chebyshev_basis(i-1,tk_xj_new); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
end
CW_ns=CW_ns';
u1=(CW_ns*ua_iter(:,1));
dvkut_xt=zeros(n,n);
dfkut_xt=zeros(n,n);
for i=1:length(tk_xj)
    for j=1:(length(tk_xj))
        dvkut_xt(i,j)=dvkut(tk_xj(i),tk_xj(j),u1(j+1));
        dfkut_xt(i,j)=dfkut(tk_xj(i),tk_xj(j),u1(j+1));
    end
end
dvkut_xt(1:length(tk_xj),1:length(tk_xj))=dvkut_xt(1:length(tk_xj),1:length(tk_xj)).*tril(ones(length(tk_xj)));
BdJx{1}=pinv(CW_ns(2:n+1,:)-dvkut_xt*(CW_ns(2:n+1,:).*length_tk_xj')-dfkut_xt*(CW_ns(2:n+1,:).*length_tk_xj'));
%% 主函数
for i=2:iter+1  
    u=(CW_ns*ua_iter(:,i-1));% (n+1*1)
    S2_x=u;
    kf=get_kf(tk_xj_new,tk_xj_new,u);%(n*1)
    g_x=S2_x(2:length(tk_xj_new))-fx(tk_xj)'-kf;%(n*1)
    delta_a=BdJx{i-1}*(-g_x);%(n*1)
    ua_iter(:,i)=ua_iter(:,i-1)+delta_a;%(n*1)
    u_new=(CW_ns*ua_iter(:,i));%(n*1)
    kf_new=get_kf(tk_xj_new,tk_xj_new,u_new);
    gx_new=u_new(2:length(tk_xj_new))-fx(tk_xj)'-kf_new;
    delta_gx=gx_new-g_x;%(n+1*1)
    BdJx{i}=BdJx{i-1}+((delta_a-BdJx{i-1}*delta_gx)*(delta_a')*BdJx{i-1})/((delta_a')*BdJx{i-1}*delta_gx);
    if norm(ua_iter(:,i)-ua_iter(:,i-1))<tol
        disp(['Converged in ', num2str(i), ' iterations.']);
        cg_it=i;
        break;
    end
end
end
%% 问题1 Fredholm and Volterra NLVIE 
function kf=get_kf(x,t,u)
temp_Volt=tril(ones(length(x)-1));
length_int=t(2:length(x))-t(1:length(x)-1);
temp_Fred=ones(length(x)-1);
ktv=-t(1:length(x)-1)'.*(u(1:length(x)-1).^2).*length_int';
ktv1=-t(2:length(x))'.*(u(2:length(x)).^2).*length_int';
    kf_V1=(temp_Volt*ktv+temp_Volt*ktv1)./2;
ktv2=(u(1:length(x)-1).^2).*length_int';
ktv3=(u(2:length(x)).^2).*length_int';
    kf_V2=(temp_Volt*ktv2+temp_Volt*ktv3).*(x(2:length(x)))'./2;
    kf_V=kf_V1+kf_V2;
ktf=t(1:length(x)-1)'.*u(1:length(x)-1).*length_int';
ktf1=t(2:length(x))'.*u(2:length(x)).*length_int';
    kf_F1=(temp_Fred*ktf+temp_Fred*ktf1)./2;
ktf2=u(1:length(x)-1).*length_int';
ktf3=u(2:length(x)).*length_int';
    kf_F2=(temp_Fred*ktf2+temp_Fred*ktf3).*(x(2:length(x)))'./2;
    kf_F=kf_F1+kf_F2;
kf=kf_V+kf_F;
end