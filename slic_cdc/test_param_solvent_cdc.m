% Path information 
% Note: the path information is set up assuming that the slic_matlab libraries is 
%       located at 'HOME' directory. If your library is located elswhere, you need
%       to change Home variable to point to the parent folder of the library 

clear all
Home = getenv('HOME');
repo_path = sprintf('%s/slic_matlab', Home);
addpath(sprintf('%s/panelbem', repo_path));
addpath(sprintf('%s/slic_cdc', repo_path));
addpath(sprintf('%s/slic_sasa', repo_path));
addpath(sprintf('%s/functions', repo_path));
addpath(sprintf('%s/ref_data', repo_path));

% loadConstants includes a bunch of useful variables and constants. also defining 
% the global variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = 0;
writeLogfile = 1;
logfileName = 'param_logfile_1';

% allData includes atom parameters for 495 neutral small molecules that are REQUIRED
% for parameterization and prediction runs. This includes dispersion-atom-types, 
% Hbond-atom-types, surface-area fractions etc.
allData = readtable('all_data.csv');

% COSMO-SAC Dispersion atom types
% all_atom_types = {'br', 'c-sp', 'c-sp2', 'c-sp3', 'cl',...
%                   'f', 'h', 'i', 'n-sp', 'n-sp3',...
%                   'n-sp2', 'o-sp2', 'o-sp3-h',...
%                   'o-sp2-n', 'o-sp3', 'p', 's'};
              
% COSMO-SAC H-bond atom types
% allHbondTypes = {'n_amine', 'n_amide', 'n_nitro',...
%                  'n_other', 'o_carbonyl', 'o_ester',...
%                  'o_nitro', 'o_hydroxyl', 'fluorine',...
%                  'h_oh', 'h_nh', 'h_other'};
             
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
			                   'staticpotential', staticpotential);

% slic-cdc calculations require a training-set with many compounds to ensure that
% all different atom-type-specific fitting parameters have enough representative
% compounds.

training_set  = ...
    {'4_bromophenol', 'ethanamide', 'teflurane', '4_chloroaniline', ...
     '2_methylpropane', '222_trifluoroethanol', '2_fluorophenol', ...
     '2_iodopropane', 'iodobenzene', '1_nitropentane', '3_cyanophenol', ...
     'pyridine', '4_nitroaniline', '14_dioxane', 'acetic_acid', ...
     'butan_1_ol', 'methyl_acetate', 'propanone', 'triethyl_phosphate', ...
     'trimethyl_phosphate', 'methanethiol', 'dimethyl_sulfate', ...
     'piperidine', 'ethylamine', 'N_methylacetamide', 'nitromethane', ...
     'nonanal', 'benzaldehyde', 'methanol', '3_methyl_1h_indole', ...
     'anthracene', '124_trimethylbenzene', '2_naphthylamine', ...
     '4_formylpyridine', 'cyclohexylamine', 'dimethyl_sulfide', ...
     'hex_1_ene', 'n_butanethiol', 'naphthalene', ...
     '33_dimethylbutan_2_one', '333_trimethoxypropionitrile', ...
     'chloroethane', 'diethyl_sulfide', 'ethene', 'imidazole', ...
     'methyl_octanoate', 'n_octane', 'n_propylbenzene', 'p_cresol', ...
     'propanoic_acid', 'tetrahydropyran', 'trichloroethene', ...
     '2_methoxyaniline', '2_methylhexane', '2_nitropropane', ...
     '26_dimethylpyridine', 'benzene', 'but_1_ene', 'but_1_yne', ...
     'm_xylene', 'methane', 'n_pentylamine', 'p_dibromobenzene'};
        
dG_list = allData.dG_expt; 
dG_disp_mob = allData.disp_mobley; 
dG_cav_mob = allData.cav_mobley; 
dG_np_mob = allData.np_mobley; 
dG_es_mob = allData.es_mobley; 
dG_np_SLIC = allData.np_SLIC; 
dG_es_SLIC = allData.es_SLIC; 
all_solutes = allData.solute;
mol_list = all_solutes;
solventAreas = allData{495, 9:79};
solventATypes = allData{495, 80:144};
solventHbData = allData{495, 145:end};
solventVdWA = allData{495, 11};
solventVdWV = allData{495, 12};
soluteAreas = allData{:,9:79};
soluteATypes = allData{:,80:144};
soluteHbData = allData{:,145:end};
soluteVdWA = allData{:,11};
soluteVdWV = allData{:,12};
temperature = 24.85 + KelvinOffset;
curdir=pwd;

for i=1:length(training_set)
  dir=sprintf('%s/ref_data/nlbc_test/%s', repo_path, training_set{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf', dir);
  chargeDist{i} = pqrData.q;
  foo = strcmp(mol_list, training_set{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\n');
    keyboard
  end
  soluteAtomAreas{i} = allData{index, 9:79};
  soluteAtomTypes{i} = {allData{index, 80:144}};
  soluteHbondData{i} = allData{index, 145:end};
  referenceData{i} = dG_list(index);
  solute_VdWA{i} = allData{index, 11};
  solute_VdWV{i} = allData{index, 12};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solventHbondData{i} = solventHbData;
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  atom_vols{i} = allData{index, 14};
  temp{i} = temperature;
  chdir(curdir);
  addProblemCosmo(training_set{i}, pqrAll{i}, srfFile{i}, chargeDist{i}, ...
                  referenceData{i}, soluteAtomAreas{i}, soluteAtomTypes{i}, ...
                  soluteHbondData{i}, solute_VdWV{i}, solute_VdWA{i}, ...
                  solventAtomAreas{i}, solventAtomTypes{i}, ...
                  solventHbondData{i}, solvent_VdWV{i}, solvent_VdWA{i}, ...
                  atom_vols{i}, temp{i});
end

% optimization

x0_1 = [0.227 -31.964 0.315 0.501 0.047 ... %slic es
        0.480 0.899 0.262 0.002 0.360 0.841 ... %disp
        0.466 1.000 0.788 0.990 0.601 0.178 ... %disp cont.
        0.562 0.008 0.853 0.840 0.429 0.408 ... %disp cont.
        0.906 0.369 0.529 0.575 0.568 0.229 ... %hbond
        0.397 0.407 0.780 0.446 0.671 0.614 ... %hbond 
        2 ... % combinatorial z
        1]; % cavity   
x0_2 = [0.461 -17.318 0.837 0.538 0.539 ... %slic es
        0.669 0.422 0.453 0.056 0.978 0.469 ... %disp
        0.517 0.539 0.310 0.721 0.486 0.430 ... %disp cont.
        0.788 0.380 0.552 0.653 0.889 0.895 ... %disp cont.
        0.626 0.454 0.000 0.481 0.870 0.277 ... %hbond
        0.436 0.608 0.933 0.160 0.559 0.627 ... %hbond 
        3 ... % combinatorial z
        1.2]; % cavity   
x0_3 = [0.582 -48.195 0.489 0.514 0.523 ... %slic es
        0.054 0.010 0.478 0.539 0.461 0.732 ... %disp
        0.415 0.185 0.695 0.824 0.562 0.824 ... %disp cont.
        0.081 0.566 0.399 0.989 0.197 0.674 ... %disp cont.
        0.173 0.027 0.224 0.438 0.926 0.031 ... %hbond
        0.820 0.715 0.254 0.161 0.585 0.204 ... %hbond 
        1.5 ... % combinatorial z
        0.25]; % cavity   
x0_4 = [0.154 -45.677 0.117 0.135 0.425 ... %slic es
        0.772 0.119 0.626 0.108 0.425 0.008 ... %disp
        0.495 0.826 0.122 0.835 0.884 0.271 ... %disp cont.
        0.624 0.432 0.230 0.357 0.520 0.049 ... %disp cont.
        0.765 0.636 0.714 0.459 0.777 0.828 ... %hbond
        0.930 0.928 0.885 0.637 0.748 0.879 ... %hbond 
        2 ... % combinatorial z
        1]; % cavity   
x0_5 = [0.444 -40.067 0.094 0.563 0.649 ... %slic es
        0.557 0.799 0.707 0.866 0.981 0.053 ... %disp
        0.854 0.858 0.718 0.585 0.809 0.413 ... %disp cont.
        0.454 0.763 0.948 0.789 0.753 0.575 ... %disp cont.
        0.875 0.951 0.498 0.921 0.888 0.522 ... %hbond
        0.451 0.186 0.350 0.853 0.611 0.339 ... %hbond 
        3 ... % combinatorial z
        0.5]; % cavity   
x0_6 = [0.361 -95.930 0.933 0.378 0.567 ... %slic es
        0.733 0.981 0.214 0.495 0.377 0.020 ... %disp
        0.754 0.921 0.737 0.448 0.571 0.810 ... %disp cont.
        0.792 0.529 0.418 0.337 0.683 0.499 ... %disp cont.
        0.586 0.326 0.442 0.826 0.036 0.042 ... %hbond
        0.496 0.779 0.868 0.658 0.186 0.653 ... %hbond 
        1 ... % combinatorial z
        0.25]; % cavity   
x0_7 = [0.771 -71.754 0.612 0.621 0.582 ... %slic es
        0.733 0.981 0.214 0.495 0.377 0.020 ... %disp
        0.754 0.921 0.737 0.448 0.571 0.810 ... %disp cont.
        0.792 0.529 0.418 0.337 0.683 0.499 ... %disp cont.
        0.586 0.326 0.442 0.826 0.036 0.042 ... %hbond
        0.496 0.779 0.868 0.658 0.186 0.653 ... %hbond 
        2 ... % combinatorial z
        1]; % cavity   
x0_8 = [0.102 -36.661 0.310 0.022 0.203 ... %slic es
        1.216 1.138 1.140 1.438 0.842 0.164 ... %disp
        0.000 1.090 0.860 0.641 1.095 0.160 ... %disp cont.
        0.327 0.886 0.442 0.828 1.153 0.358 ... %disp cont.
        0.796 0.210 0.489 0.145 0.121 0.088 ... %hbond
        0.444 0.465 0.109 0.978 2.708 0.088 ... %hbond 
        3 ... % combinatorial z
        1]; % cavity   
x0_9 = [0.179 -41.754 0.987 0.876 0.392 ... %slic es
        0.823 0.891 0.891 0.669 0.466 0.676 ... %disp
        0.161 0.789 0.365 0.452 0.298 0.877 ... %disp cont.
        0.387 0.112 0.814 0.774 0.373 0.113 ... %disp cont.
        0.188 0.249 0.717 0.437 0.448 0.450 ... %hbond
        0.888 0.225 0.522 0.506 0.215 0.646 ... %hbond 
        4 ... % combinatorial z
        0.25]; % cavity   
x0_10 = [0.454 -48.813 -0.541 -0.548 -0.062 ... %slic es
        0.463 0.787 0.777 0.785 0.807 0.552 ... %disp
        0.455 0.354 0.365 0.917 0.835 0.281 ... %disp cont.
        0.756 0.031 0.622 0.484 0.992 0.417 ... %disp cont.
        0.642 0.355 0.401 0.036 0.666 0.506 ... %hbond
        0.319 0.539 0.426 0.652 0.314 0.668 ... %hbond 
        2 ... % combinatorial z
        1]; % cavity   
x0_11 = [0.453 -48.813	-0.541 -0.548	-0.062	... %slic es
        1.216	1.138	1.140	1.438	0.842	0.164 ... %disp
        0.000	1.090	0.860	0.641	1.095	0.160 ... %disp cont.
        0.327	0.886	0.442	0.828	1.153	0.358 ... %disp cont.
        0.796	0.210	0.489	0.145	0.121	0.088 ... %hbond
        0.444	0.465	0.109	0.978	2.708	0.088 ... %hbond 
        2 ... % combinatorial z
        1]; % cavity 			

x0 = x0_1;

opt_file_name = 'OptSlicCdc_1.mat';

% alpha : x(1)                      % O-sp3-H dispersion coeff : x(20)     
% beta : x(2)                       % P dispersion coeff : x(21)    
% gamma : x(3)                      % S dispersion coeff : x(22)     
% mu : x(4)                         % q_s H-bond coeff : x(23)  
% phi_static : x(5)                 % n_amine H-bond coeff : x(24)          
% Br dispersion coeff : x(6)        % n_amide H-bond coeff : x(25)
% C-sp dispersion coeff : x(7)      % n_nitro H-bond coeff : x(26)  
% C-sp2 dispersion coeff : x(8)     % n_other H-bond coeff : x(27)   
% C-sp3 dispersion coeff : x(9)     % o_carbonyl H-bond coeff : x(28)   
% Cl dispersion coeff : x(10)       % o_ester H-bond coeff : x(29) 
% F dispersion coeff : x(11)        % o_nitro H-bond coeff : x(30)
% H dispersion coeff : x(12)        % o_hydroxyl H-bond coeff : x(31)
% I dispersion coeff : x(13)        % fluorine H-bond coeff : x(32)
% N-sp dispersion coeff : x(14)     % h_oh H-bond coeff : x(33)   
% N-sp2 dispersion coeff : x(15)    % h_nh H-bond coeff : x(34)    
% N-sp3 dispersion coeff : x(16)    % h_other H-bond coeff : x(35)    
% O-sp2 dispersion coeff : x(17)    % z combinatorial coeff : x(36)    
% O-sp2-N dispersion coeff : x(18)  % cavity rescaling coeff : x(37)      
% O-sp3 dispersion coeff : x(19)

% upper bound
ub = [+2 +200 +100 +20  +0.1 ...
      20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 ...
      20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 4];

% lower bound
lb = [-2 -200 -100 -20  -0.1 ...
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% optimization options
options = optimoptions('lsqnonlin', 'MaxIter', 30);
options = optimoptions(options,'Display', 'iter');

% objective function (SLIC_es + CDC_np + hb)
y = @(x)ObjectiveFromBEMCosmo(x);

[x, resnorm, residual, exitflag, output,] = lsqnonlin(y, x0, lb, ub, ...
    options);
[err, calc, ref, es, np, hb, disp, disp_slsl, disp_svsl, disp_svsv, cav, ...
    comb] = ObjectiveFromBEMCosmo(x);

[err0, calc0, ref0, es0, np0, hb0, disp0, disp_slsl0, disp_svsl0, ...
    disp_svsv0, cav0, comb0] = ObjectiveFromBEMCosmo(x0);

[~, id]=ismember(training_set, mol_list);
disp_mob = allData.disp_mobley(id); 
cav_mob = allData.cav_mobley(id); 
np_mob = allData.np_mobley(id); 
es_mob = allData.es_mobley(id); 
np_SLIC = allData.np_SLIC(id); 
es_SLIC= allData.es_SLIC(id);
rmse = rms(calc-ref);
rmse_np = rms(np_mob-np);
rmse_disp = rms(disp_mob-disp);
rmse_cav = rms(cav_mob-cav);
rmse_es = rms(es_mob-es);
rmse_eshb = rms(es_mob-es-hb);

% save the results
save(opt_file_name, 'x', 'training_set', 'mol_list', 'ref', ...
     'calc', 'es', 'np', 'hb', 'disp', ...
     'disp_slsl', 'disp_svsl', 'disp_svsv', 'comb', 'cav', ...
     'disp_mob', 'cav_mob', 'np_mob', 'es_mob', 'np_SLIC', ...
     'rmse', 'rmse_np', 'rmse_disp', 'rmse_cav', 'rmse_es', 'rmse_eshb', ...
     'x0', 'calc0', 'es0', 'np0', 'hb0', 'disp0', 'disp_slsl0', ...
     'disp_svsl0', 'disp_svsv0', 'comb0', 'cav0', 'epsOut');
