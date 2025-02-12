clc;close all;clear all;
n=128;
fx=@(x)exp(x) + exp(-x);
kf=@(x,t)-exp(-x-t);
RF=@(x)exp(x);
%% 求解n配置点的解和误差
[ua_iter1,tk_xj]=LFIE_CWM(fx, kf, n);
pointwise_error = PE1_CWM(RF,ua_iter1,256);
%% 求解n/2配置点的解和误差
[ua_iter2,tk_xj2]=LFIE_CWM(fx, kf, n/2);
pointwise_error2 = PE1_CWM(RF,ua_iter2,256);
%% 速度
SP=log(pointwise_error2./pointwise_error)/log(2);