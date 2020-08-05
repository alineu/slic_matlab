function Parameters = Parameters_of_temp(x)

addpath('..')
addpath('../born/')
addpath('../../pointbem/')

x0 = [0.5 -60 -0.5];
lb = [0 -500 -100];
ub = [200 0 0];

options = optimoptions('lsqnonlin');
options = optimoptions(options,'Display', 'off');

t = 5:6;


for i = 1:length(t)
   
Parameters = ObjectiveFunction_MSA(x,t(i));

[x(i,:),resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@Parameters_of_temp,x0,lb,ub,options);

end
end
