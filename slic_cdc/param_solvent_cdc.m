% Path information 
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

% allData includes atom parameters for 495 neutral small molecules that are REQUIRED
% for parameterization and prediction runs. This includes dispersion-atom-types, 
% Hbond-atom-types, surface-area fractions etc.
allData = readtable('all_data.csv');

% COSMO-SAC Dispersion atom types
% all_atom_types = {'br','c-sp','c-sp2','c-sp3','cl',...
%                   'f','h','i','n-sp','n-sp3',...
%                   'n-sp2','o-sp2','o-sp3-h',...
%                   'o-sp2-n','o-sp3','p','s'};
              
% COSMO-SAC H-bond atom types
% allHbondTypes = {'n_amine','n_amide','n_nitro',...
%                  'n_other','o_carbonyl','o_ester',...
%                  'o_nitro','o_hydroxyl','fluorine',...
%                  'h_oh','h_nh','h_other'};
             
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

% slic-cdc calculations require a training-set with many compounds to ensure that
% all different atom-type-specific fitting parameters have enough representative
% compounds.

training_set = {'4_bromophenol', 'ethanamide', 'teflurane', '4_chloroaniline',...
                '2_methylpropane', '222_trifluoroethanol', '2_fluorophenol',...
                '2_iodopropane', 'iodobenzene', '1_nitropentane', '3_cyanophenol',...
                'pyridine','4_nitroaniline','14_dioxane','acetic_acid','butan_1_ol',...
                'methyl_acetate','propanone','triethyl_phosphate','trimethyl_phosphate',...
                'methanethiol','dimethyl_sulfate','piperidine','ethylamine','N_methylacetamide',...
                'nitromethane','nonanal','benzaldehyde','methanol',...
                '3_methyl_1h_indole','anthracene','124_trimethylbenzene','2_naphthylamine'...
                ,'4_formylpyridine','cyclohexylamine','dimethyl_sulfide','hex_1_ene',...
                'n_butanethiol','naphthalene','33_dimethylbutan_2_one','333_trimethoxypropionitrile',...
                'chloroethane','diethyl_sulfide','ethene','imidazole',...
                'methyl_octanoate','n_octane','n_propylbenzene',...
                'p_cresol','propanoic_acid','tetrahydropyran','trichloroethene',...
                '2_methoxyaniline','2_methylhexane','2_nitropropane','26_dimethylpyridine',...
                'benzene','but_1_ene','but_1_yne','m_xylene','methane',...
                'n_pentylamine','p_dibromobenzene'};
        
dG_list = allData.dG_expt; 
dG_disp_mob = allData.disp_mobley; 
dG_cav_mob = allData.cav_mobley; 
dG_np_mob = allData.np_mobley; 
dG_es_mob = allData.es_mobley; 
dG_np_SLIC = allData.np_SLIC; 
dG_es_SLIC = allData.es_SLIC; 
all_solutes = allData.solute;
mol_list = all_solutes;
solventAreas = allData{495,9:79};
solventATypes = allData{495,80:144};
solventHbData = allData{495,145:end};
solventVdWA = allData{495,11};
solventVdWV = allData{495,12};
soluteAreas = allData{:,9:79};
soluteATypes = allData{:,80:144};
soluteHbData = allData{:,145:end};
soluteVdWA = allData{:,11};
soluteVdWV = allData{:,12};
temperature = 24.85 + KelvinOffset;
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
  soluteAtomAreas{i} = allData{index,9:79};
  soluteAtomTypes{i} = {allData{index,80:144}};
  soluteHbondData{i} = allData{index,145:end};
  referenceData{i} = dG_list(index);
  solute_VdWA{i} = allData{index,11};
  solute_VdWV{i} = allData{index,12};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solventHbondData{i} = solventHbData;
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  atom_vols{i} = allData{index,14};
  temp{i} = temperature;
  chdir(curdir);
  addProblemCosmo(training_set{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},...
                  soluteAtomAreas{i},soluteAtomTypes{i},soluteHbondData{i},...
                  solute_VdWV{i},solute_VdWA{i},...
                  solventAtomAreas{i},solventAtomTypes{i},solventHbondData{i},...
                  solvent_VdWV{i},solvent_VdWA{i},...
                  atom_vols{i},temp{i});
end

% optimization

% initial guesses
x0  =  [0.453   -48.813	-0.541	-0.548	-0.062	... %slic es
        1.216	1.138	1.140	1.438	0.842	0.164 ... %disp
        0.000	1.090	0.860	0.641	1.095	0.160 ... %disp cont.
        0.327	0.886	0.442	0.828	1.153	0.358 ... %disp cont.
        0.796	0.210	0.489	0.145	0.121	0.088 ... %hbond
        0.444	0.465	0.109	0.978	2.708	0.088 ... %hbond 
        2.430 ... % combinatorial z
        2.646]; % cavity 		

% upper bound
ub = [+2 +200 +100 +20  +0.1 ...
      3*x0(6:end)];

% lower bound
lb = [-2 -200 -100 -20  -0.1 ...
      0.4*x0(6:end)];

% optimization options
options = optimoptions('lsqnonlin','MaxIter',6);
options = optimoptions(options,'Display', 'iter');

% objective function (SLIC_es + CDC_np + hb)
y = @(x)ObjectiveFromBEMCosmo(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np,hb,disp,disp_slsl,disp_svsl,disp_svsv,cav,comb]=ObjectiveFromBEMCosmo(x);
[err0,calc0,ref0,es0,np0,hb0,disp0,disp_slsl0,disp_svsl0,disp_svsv0,cav0,comb0]=ObjectiveFromBEMCosmo(x0);
[~,id]=ismember(training_set,mol_list);
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
save('OptSlicCdcWater.mat','x','training_set','mol_list','ref','calc','es','np','hb','disp',...
     'disp_slsl','disp_svsl','disp_svsv','comb','cav',...
     'disp_mob','cav_mob','np_mob','es_mob','np_SLIC',...
     'rmse','rmse_np','rmse_disp','rmse_cav','rmse_es','rmse_eshb',...
     'x0','calc0','es0','np0','hb0','disp0', 'disp_slsl0','disp_svsl0',...
     'disp_svsv0','comb0','cav0','epsOut');
