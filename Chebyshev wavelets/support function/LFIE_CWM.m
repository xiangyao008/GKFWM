function [ua_iter,tk_xj] = LFIE_CWM(fx, kf, n)
% linear Fredholm integral equation
% fx Ϊ��������, ����x���� ����Ϊ�������
% kf Ϊ�����������ʾ�˺���������x,t��������
% n Ϊ������iter,tol��ʾ������ʱ���ݲ�
% ux_iter Ϊ ÿ�ε������ֵ
% cg_it��ʾ����ʱ�Ľ״�
%% Ԥ�������
chebyshev_basis = @(n, t) cos(n * acos(2*t-1)); % Chebyshev����ʽT_n(t)
tk_xj=CP_CWM(n);
tk_xj_new=[0,tk_xj]; %0��ͷ
length_tk_xj=tk_xj-tk_xj_new(:,1:length(tk_xj));
% tk_xj=fliplr(tk_xj);
CW_ns=zeros(n,n);
kf_xt=zeros(n,n);
for i=1:n
    CW_ns(i,:)=chebyshev_basis(i-1,tk_xj); % P_ns��ÿһ�б�ʾn�����õ��i��s�θ߽�HAARС������ֵ
end
CW_ns=CW_ns';
for i=1:n
    for j=1:n
        kf_xt(i,j)=kf(tk_xj(i),tk_xj(j));
    end
end
fx_x=fx(tk_xj)';
%% ������
ua_iter=pinv(CW_ns-kf_xt*CW_ns./length(tk_xj))*fx_x;
end
