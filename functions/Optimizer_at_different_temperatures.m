function [x,resnorm,residual,exitflag,output,lambda,jacobian] = Optimizer_at_different_temperatures(x0,lb,ub)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('lsqnonlin');
%% Modify options setting
options = optimoptions(options,'Display', 'off');

[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@Parameters_of_temp,x0,lb,ub,options);

