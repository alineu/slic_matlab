% Note: the path information is set up assuming that the slic_matlab libraries is 
%       located at 'HOME' directory. If your library is located elswhere, you need
%       to change Home variable to point to the parent folder of the library 
close all
clear all

% Path information
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
saveMemory = 1;
writeLogfile = 1;
logfileName = 'run_logfile';
% epsIn : dielectric constant of the solute
epsIn  =  1;
% epsOut : dielectric constant of the solvent
epsOut = 78.34;% from mnsol Database
% conv_factor is a pre-factor that results in solvation energies in kcal/mol units
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn', epsIn, 'epsOut', epsOut,...
                         'kappa', kappa, 'conv_factor', conv_factor,...
                         'staticpotential',staticpotential);

ParamWatInfo = load('OptWater_wo_ion.mat');
x = ParamWatInfo.x;

% load the test-set data; 
allData = readtable('../ref_data/all_data.csv'); 
all_solutes = allData.solute;
all_surfAreas = allData.SASA;
mol_list = all_solutes;
dG_list = allData.dG_expt;
[m, index] = ismember(mol_list,all_solutes);
surfArea_list = all_surfAreas(index);

curdir = pwd;
for i=1:length(mol_list)
  dir=sprintf('%s/ref_data/nlbc_test/%s',repo_path,mol_list{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;
  referenceData{i} = dG_list(i);
  surfArea{i} = surfArea_list(i);
  chdir(curdir);
  addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end

[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x);
save('RunSlicSasaWater','errfinal','calcE','refE','es','np');
