function [ua_iter,tk_xj,cg_it] = NLVIE_CWM_Simpson(fx, dkut, n, iter,tol)
% nonlinear Fredholm integral equation
% fx 为给定函数, 接受x参数 定义为函数句柄
% kuf 为已知积分核, 接受x,t,u三个参数 定义为函数句柄,包含被解函数的函数
% dkut表示kuf对被解函数u(t)的偏导
% n 为阶数，iter,tol表示迭代此时和容差
% ux_iter 为 每次迭代后的值
% cg_it表示收敛时的阶次
% load('true_coff.mat','Array_fx1');
%% 预分配参数
chebyshev_basis = @(n, t) cos(n * acos(2*t-1)); % Chebyshev多项式T_n(t)\
% chebyshev_basis = @(n, t) 2*cos(n * acos(2*t-1))./sqrt(sqrt(1-(2*t-1).^2)); % Chebyshev多项式T_n(t)\
tk_xj=CP_CWM(n);
tk_xj_new=[0,tk_xj]; %0开头
tk_xj_new2=(tk_xj_new(1:n)+tk_xj_new(2:n+1))./2;
length_tk_xj=tk_xj-tk_xj_new(:,1:length(tk_xj));
%% 预计算参数
Init_ux = zeros(n, 1);
Init_ux (1,1)=2;
% Init_ux (3,1)=-0.5;
ua_iter=zeros(n,iter+1);
ua_iter(:,1)=Init_ux;
for i=1:n
    CW_ns(i,:)=chebyshev_basis(i-1,tk_xj_new); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
end
CW_ns=CW_ns';% P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
for i=1:n
    CW_ns2(i,:)=chebyshev_basis(i-1,tk_xj_new2); % P_ns的每一行表示n个配置点的i阶s次高阶HAAR小波的数值
end
CW_ns2=CW_ns2';
u1=(CW_ns*ua_iter(:,1));
% kf1=get_dkf(tk_xj,tk_xj,u1);
% temp=repmat(kf1, 1, n);
% BdJx{1}=pinv(CW_ns-temp.*CW_ns);
dkut_xt=zeros(n,n);
for i=1:length(tk_xj)
    for j=1:(length(tk_xj))
        dkut_xt(i,j)=dkut(tk_xj(i),tk_xj(j),u1(j+1));
    end
end
dkut_xt(1:length(tk_xj),1:length(tk_xj))=dkut_xt(1:length(tk_xj),1:length(tk_xj)).*tril(ones(length(tk_xj)));
BdJx{1}=pinv(CW_ns(2:n+1,:)-dkut_xt*CW_ns(2:n+1,:).*length_tk_xj);
%% 主函数
for i=2:iter+1  
    u=(CW_ns*ua_iter(:,i-1));% (n+1*1)
    u2=(CW_ns2*ua_iter(:,i-1));
    S2_x=u;
    kf=get_kf(tk_xj_new,tk_xj_new,tk_xj_new2,u,u2);%(n*1)
    g_x=S2_x(2:length(tk_xj_new))-fx(tk_xj)'-kf;%(n*1)
    delta_a=-BdJx{i-1}*(g_x);%(n*1)
    ua_iter(:,i)=ua_iter(:,i-1)+delta_a;%(n*1)
    if norm(ua_iter(:,i)-ua_iter(:,i-1))<tol
        disp(['Converged in ', num2str(i), ' iterations.']);
        cg_it=i;
        break;
    end
    u_new=(CW_ns*ua_iter(:,i));%(n*1)
    u2_new=(CW_ns2*ua_iter(:,i));
    kf_new=get_kf(tk_xj_new,tk_xj_new,tk_xj_new2,u_new,u2_new);    
    gx_new=u_new(2:length(tk_xj_new))-fx(tk_xj)'-kf_new;
    delta_gx=gx_new-g_x;%(n+1*1)
    % muk=1+(delta_gx'*BdJx{i-1}*delta_gx)./(delta_a'*delta_gx);
    % BdJx{i}=BdJx{i-1}+(muk*(delta_a*delta_a')-(delta_a*delta_gx')*BdJx{i-1}-BdJx{i-1}*delta_gx*delta_a')...
    % /((delta_a')*delta_gx);% 秩2方法
    BdJx{i}=BdJx{i-1}+((delta_a-BdJx{i-1}*delta_gx)*(delta_a')*BdJx{i-1})/((delta_a')*BdJx{i-1}*delta_gx);% Sherman–Morrison formula
    % BdJx{i}=BdJx{i-1}+(delta_a-BdJx{i-1}*delta_gx)*(delta_gx')/(norm(delta_gx)^2);
    % BdJx{i}=BdJx{i-1}+((delta_a-BdJx{i-1}*delta_gx)*(delta_a-BdJx{i-1}*delta_gx)')/((delta_a-BdJx{i-1}*delta_gx)'*delta_gx);
    % BdJx{i}=BdJx{i-1}+((delta_gx-BdJx{i-1}*delta_a)*(delta_a'))./((delta_a')*delta_a);%不求逆
end
end
%% 问题1
% function kf=get_kf(x,t,u)
% temp=tril(ones(length(x)-1),-1);
% temp2=[temp;ones(1,length(x)-1)];
% kt=(1./(1+u(1:length(x)-1).^2)).*(t(2:length(x))-t(1:length(x)-1))';
% kt1=(1./(1+u(2:length(x)).^2)).*(t(2:length(x))-t(1:length(x)-1))';
% kf=(temp2*kt+temp2*kt1)./2;
% end
%% 问题2
% function kf=get_kf(x,t,t2,u,u2)
% temp=tril(ones(length(x)-1),-1);
% temp2=[temp;ones(1,length(x)-1)];
% kt=(t(1:length(x)-1)'.^2).*(u(1:length(x)-1).^2).*(t(2:length(x))-t(1:length(x)-1))';
% kt1=(t2'.^2).*(u2.^2).*(t(2:length(x))-t(1:length(x)-1))';
% kt2=(t(2:length(x))'.^2).*(u(2:length(x)).^2).*(t(2:length(x))-t(1:length(x)-1))';
% kf=(temp2*kt+4*temp2*kt1+temp2*kt2).*x'./6;
% end
%% 问题3
% function kf=get_kf(x,t,t2,u,u2)
% temp=tril(ones(length(x)-1),-1);
% temp2=[temp;ones(1,length(x)-1)];
% kt=(exp(-2.*t(1:length(x)-1)')).*(u(1:length(x)-1).^2).*(t(2:length(x))-t(1:length(x)-1))';
% kt1=exp(-2.*t2').*(u2.^2).*(t(2:length(x))-t(1:length(x)-1))';
% kt2=(exp(-2.*t(2:length(x))')).*(u(2:length(x)).^2).*(t(2:length(x))-t(1:length(x)-1))';
% kf=(temp2*kt+4*temp2*kt1+temp2*kt2).*exp(-x)'./6;
% end
%% 问题4 Third kind of NLVIE 
function kf=get_kf(x,t,t2,u,u2)
length_int=t(2:length(x))-t(1:length(x)-1);
temp=tril(ones(length(x)-1));
kt=(u(1:length(x)-1).^2+u(1:length(x)-1)).*length_int';
kt1=(u2.^2+u2).*length_int';
kt2=(u(2:length(x)).^2+u(2:length(x))).*length_int';
kf=(temp*kt+4*temp*kt1+temp*kt2).*(-1)./6;
end
%% 问题5 
% function kf=get_kf(x,t,t2,u,u2)
% temp=tril(ones(length(x)-1),-1);
% temp2=[temp;ones(1,length(x)-1)];
% kt=t(1:length(x)-1)'.*(u(1:length(x)-1).^2).*(t(2:length(x))-t(1:length(x)-1))';
% kt1=t2'.*(u2.^2).*(t(2:length(x))-t(1:length(x)-1))';
% kt2=t(2:length(x))'.*(u(2:length(x)).^2).*(t(2:length(x))-t(1:length(x)-1))';
% kf=(temp2*kt+4*temp2*kt1+temp2*kt2).*(x)'./6;
% end