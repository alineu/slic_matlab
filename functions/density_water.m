function density = density_water(temp)

% Returns water density for temperature range 25-125 C
% data from Jorgensen 07

temp_list = [273.18	283.18 288.78 294.18 299.88 305.38 310.98 ... 
             322.08	333.18 344.28 355.38 366.48	373.18];
density_list  = [0.99987 0.99975 0.99907 0.99802 0.99669 0.9951	0.99318 ...
             0.9887	0.98338	0.97729	0.97056	0.96333	0.95865];

polycoeffs = polyfit(temp_list,density_list,2);
if temp < min(temp_list)
  fprintf('Warning in density_water: requested temp %f < T_min %f!\n',temp,min(temp_list));
end
if temp > max(temp_list)
  fprintf('Warning in density_water: requested temp %f > T_max %f!\n',temp,max(temp_list));
end
density = polyval(polycoeffs,temp);