%% Higher order Haar wavelet method for numerical solution of integral equations
% 参数命名定义形式参考这篇文章
%%
clc;close all;clear all;
n=64;
fx=@(x)exp(x) + exp(-x);
kf=@(x,t)-exp(-x-t);
RF=@(x)exp(x);
%% 求解n配置点的解和误差
[Array_fx1,tk_xj1]=LFIE_HAAR_S2(fx, kf, n);
pointwise_error = PE2_HOHWM(RF,Array_fx1,256);
%% 求解n/2配置点的解和误差
[Array_fx2,tk_xj2]=LFIE_HAAR_S2(fx, kf, n/2);
pointwise_error2 = PE2_HOHWM(RF,Array_fx2,256);
% [Array_se2,Array_fx2,tk_xj2]=LFIE(fx, kf, n/2);
% % 伪逆求解线性方程
% Array_se_inv2 = pinv(Array_se2); 
% result2=Array_se_inv2*Array_fx2;
% % A\B求解
% % result2=Array_se2\Array_fx2;
% pointwise_error2 = PE(RF,result2,2*n);
% %% 速度
SP=log(pointwise_error2./pointwise_error)/log(2);