addpath('..')
addpath('../born/')
addpath('../../pointbem/')

% Initial parameters and bounds
x0 = [0.5 -60 -0.5];
lb = [0 -Inf -Inf];
ub = [Inf 0 0];

% Range of temperature values in Celcius
t = 0.1:.1:99;

options = optimoptions('lsqnonlin');
options = optimoptions(options,'Display', 'off');

x(length(t),3) = zeros;


for i = 1:length(t)
y = @(x)ObjectiveFunction_MSA(x,t(i));
[x(i,:),resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(y,x0,lb,ub,options);
display(['Optimized at ' , num2str(t(i)), ' degrees Celcius'])
end

save('Parameters','x','t')