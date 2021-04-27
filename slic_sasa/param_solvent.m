q% Path information 
% Note: the path information is set up assuming that the slic_matlab libraries is 
%       located at 'HOME' directory. If your library is located elswhere, you need
%       to change Home variable to point to the parent folder of the library 

clear all
Home = getenv('HOME');
repo_path = sprintf('%s/slic_matlab',Home);
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/slic_cdc',repo_path));
addpath(sprintf('%s/slic_sasa',repo_path));
addpath(sprintf('%s/functions',repo_path));
addpath(sprintf('%s/ref_data',repo_path));

% loadConstants includes a  bunch of useful variables and constants. also defining 
% the global variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = 0;
writeLogfile = 0;
logfileName = 'param_logfile';

% epsIn : dielectric constant of the solute
epsIn  =  1;
% epsOut : dielectric constant of the solvent
epsOut = 78.34;% from mnsol Database
KelvinOffset = 273.15;
% conv_factor is a pre-factor that results in solvation energies in kcal/mol units
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!

UsefulConstants = struct('epsIn', epsIn, 'epsOut', epsOut,...
                         'kappa', kappa, 'conv_factor', conv_factor,...
			                   'staticpotential',staticpotential);

training_set = {'methane','propane','n_butane','2_methylpropane','toluene',...
                'methanol','ethanol','p_cresol','methanethiol',...
                'methyl_ethyl_sulfide','ethanamide','3_methyl_1h_indole'};

Data = readtable('../ref_data/thermo.csv');
% Experimental data for 502 neutral molecules and 9 monovalent ions
% from multiple references - See '../ref_data/thermo_expt.xls' for more details
% Var1 = solutes | 2 = dG_expt | 3 = Surface_Area | 4 = Volume 
%    5 = dG_np   | 6 = dG_es   | 7 = dG_Mobley

mol_list = Data.Var1;
dG_list = Data.Var2;
surfArea_list = Data.Var3;
% [m, index] = ismember(training_set,mol_list);
% surfArea_list = all_surfAreas(index);             

curdir=pwd;
for i=1:length(training_set)
  dir=sprintf('%s/ref_data/nlbc_test/%s',repo_path,training_set{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;
  foo = strcmp(mol_list,training_set{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\n');
    keyboard
  end
  referenceData{i} = dG_list(index);
  surfArea{i} = surfArea_list(index);
  chdir(curdir);
  addProblemSA(training_set{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)  0 -0.03 1.6];
lb = [-2 -200 -100 -20  -0.1  -0.1  -4];
ub = [+2 +200 +100 +20  +0.1  +0.1  +4];

options = optimoptions('lsqnonlin','MaxIter',8);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);
rmse = rms(calc-ref);
save(OptFileName,'x','rmse','ref','calc','es','np','x0','calc0','es0','np0','epsOut');

