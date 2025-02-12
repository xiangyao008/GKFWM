clc;close all;clear all;
n=32;
fx=@(x)1/2*(x.^2).*exp(-x);
kf=@(x,t)1/2*((x-t).^2).*(exp(-x+t));
RF=@(x)1/3-1/3*exp(-(3/2)*x).*(cos(sqrt(3)/2*x)+sqrt(3)*sin(sqrt(3)/2*x));
%% 求解n配置点的解和误差
[ua_iter1,tk_xj]=LVIE_K1(fx, kf, n);
pointwise_error = PE1(RF,ua_iter1,256);
%% 求解n/2配置点的解和误差
[ua_iter2,tk_xj2]=LVIE_K1(fx, kf, n/2);
pointwise_error2 = PE1(RF,ua_iter2,256);
%% 速度
SP=log(pointwise_error2./pointwise_error)/log(2);