%% Higher order Haar wavelet method for numerical solution of integral equations
clc;close all;clear all;
n=64;iter=100;tol = 1e-6; 
fx=@(x)sin(pi*x);
kf=@(x,t)1/5.*cos(pi*x).*sin(pi*t);
kut=@(x,t,u)1/5.*cos(pi*x).*sin(pi*t).*(u.^3);
dkut=@(x,t,u)3/5.*cos(pi*x).*sin(pi*t).*(u.^2);
RF=@(x)sin(pi*x)+1/3.*(20-sqrt(391)).*cos(pi*x);
%% 求解n配置点的解和误差
[ua_iter1,tk_xj1,cg_it1]=NLFIE_CWM(fx,dkut, n, iter,tol);
pointwise_error1 = PE1_CWM(RF,ua_iter1(:,cg_it1),256);
 %% 求解n/2配置点的解和误差6
[ua_iter2,tk_xj2,cg_it2]=NLFIE_CWM(fx,dkut, n/2, iter,tol);
pointwise_error2 = PE1_CWM(RF,ua_iter2(:,cg_it1),256);

%% 速度
SP=log(pointwise_error2./pointwise_error1)/log(2);