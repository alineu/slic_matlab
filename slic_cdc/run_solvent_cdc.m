% Path information 
% Note: the path information is set up assuming that the slic_matlab libraries is 
%       located at 'HOME' directory. If your library is located elswhere, you need
%       to change Home variable to point to the parent folder of the library 
% Running this script generates SLIC/CDC predictions of hydration free energies 
% Running this script requires the optimized parameter file 
% An example is included in the current folder (OptSlicCdcWater.mat)

clear all
Home = getenv('HOME');
% The line below assumes that the repository is located in your HOME
% directory (/Users/user in MacOS). Otherwise you need to replace Home with
% 'path/to/parent/fodler'
repo_path = sprintf('%s/slic_matlab', Home);
addpath(sprintf('%s/panelbem', repo_path));
addpath(sprintf('%s/slic_cdc', repo_path));
addpath(sprintf('%s/slic_sasa', repo_path));
addpath(sprintf('%s/functions', repo_path));
% Ask Ali for ref_data
% https://www.dropbox.com/sh/5okqykiw8dr6gmb/AAAgcYlkuqp0lQIWcuYzGYgZa?dl=0
addpath(sprintf('%s/ref_data', repo_path));

% loadConstants includes a  bunch of useful variables and constants. also defining 
% the global variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = 1;
writeLogfile = 0;
logfileName = 'run_logfile';

% allData includes atom parameters for 495 neutral small molecules that are REQUIRED
% for parameterization and prediction runs. This includes dispersion-atom-types, 
% Hbond-atom-types, surface-area fractions etc.
allData = readtable('all_data_2.csv');

% COSMO-SAC Dispersion atom types
% all_atom_types = {'br', 'c-sp', 'c-sp2', 'c-sp3', 'cl', ...
%                   'f', 'h', 'i', 'n-sp', 'n-sp3', ...
%                   'n-sp2', 'o-sp2', 'o-sp3-h', ...
%                   'o-sp2-n', 'o-sp3', 'p', 's'};
              
% COSMO-SAC H-bond atom types
% allHbondTypes = {'n_amine', 'n_amide', 'n_nitro', ...
%                  'n_other', 'o_carbonyl', 'o_ester', ...
%                  'o_nitro', 'o_hydroxyl', 'fluorine', ...
%                  'h_oh', 'h_nh', 'h_other'};
             
% epsIn : dielectric constant of the solute
epsIn  =  1;

% epsOut : dielectric constant of the solvent
epsOut = 78.34;% from mnsol Database

% conv_factor is a pre-factor that results in solvation energies in kcal/mol units
conv_factor = 332.112;
KelvinOffset = 273.15;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!

UsefulConstants = struct('epsIn', epsIn, 'epsOut', epsOut, ...
                         'kappa', kappa, 'conv_factor', conv_factor, ...
                         'staticpotential', staticpotential);

% load test-set data; mol_list is the test test solutes
allData = readtable('all_data.csv'); 
mol_list = allData.solute;
dG_list = allData.dG_expt;
solventAreas = allData{495, 9:79};
solventATypes = allData{495, 80:144};
solventVdWA = allData{495, 11};
solventVdWV = allData{495, 12};
temperature = 24.85 + KelvinOffset;
curdir = pwd;
for i=1:length(mol_list)
  dir=sprintf('%s/ref_data/nlbc_test/%s', repo_path, mol_list{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf', dir);
  chargeDist{i} = pqrData.q;
  soluteAtomAreas{i} = allData{i, 9:79};
  soluteAtomTypes{i} = {allData{i, 80:144}};
  hbondData{i} = allData{i, 145:154};
  solute_VdWA{i} = allData{i, 11};
  solute_VdWV{i} = allData{i, 12};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  spherocity{i} = allData{i, 155};
  atom_vols{i} = allData{i, 14};
  temp{i} = temperature;
  referenceData{i} = dG_list(i);
  chdir(curdir);
  addProblemCosmo_2(mol_list{i}, pqrAll{i}, srfFile{i}, chargeDist{i}, referenceData{i}, ...
                    soluteAtomAreas{i}, soluteAtomTypes{i}, hbondData{i}, ...
                    solute_VdWV{i}, solute_VdWA{i}, ...
                    solventAtomAreas{i}, solventAtomTypes{i}, ...
                    solvent_VdWV{i}, solvent_VdWA{i}, ...
                    atom_vols{i}, temp{i});
end

disp_mob = allData.disp_mobley; 
cav_mob = allData.cav_mobley; 
np_mob = allData.np_mobley; 
es_mob = allData.es_mobley; 
np_SLIC = allData.np_SLIC; 
es_SLIC= allData.es_SLIC;

% load the optimal parameters from the optimization run (param_solvent)
chdir(curdir);
file_name = 'OptSlicCdcWater.mat';
ParamInfo = load(file_name);
training_set = ParamInfo.training_set;

% load the optimal parameters that are obtained from the optimization process
x = ParamInfo.x;

% objective function
[err, calc, ref, es, np, hb, disp, disp_slsl, disp_svsl, disp_svsv, cav, comb] = ObjectiveFromBEMCosmo_2(x);
rmse = rms(ref-calc);
rmse_np = rms(np_mob-np);
rmse_disp = rms(disp_mob-disp);
rmse_cav = rms(cav_mob-cav);
rmse_es = rms(es_mob-es);
rmse_eshb = rms(es_mob-es-hb);

% save the results
save('RunSlicCdcWater.mat', 'mol_list', 'training_set', 'x', ...
    'err', 'calc', 'ref', 'es', 'np', 'hb', ...
    'disp', 'disp_slsl', 'disp_svsl', 'disp_svsv', 'cav', 'comb', ...
    'rmse', 'rmse_np', 'rmse_disp', 'rmse_cav', 'rmse_es', 'rmse_eshb', ...
    'disp_mob', 'cav_mob', 'np_mob', 'es_mob', 'np_SLIC', 'es_SLIC');
