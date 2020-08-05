% Note: the path information is set up assuming that the slic_matlab libraries is 
%       located at 'HOME' directory. If your library is located elswhere, you need
%       to change Home variable to point to the parent folder of the library 
clc
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ionflag=0; % ionflag=0 : ions data are not included in the testset 
           % ionflag=1 : ions data are included in the testset

calcflag=1; % calcflag=0 : we load (previously calculated) dGs and then calculate other thermodynamic properties
            % calcflag=1 : we calculate dGs using BEM first and then other thermodynamic properties
                 
dataset='mobley'; % 'mobley' or 'mnsol'. In both cases we use mobley syrface areas 

temp_min= 4.85; % lower bound of the temperature interval
temp_max=44.85; % upper bound in the temperature interval
tempdiv=5; % number of temperatures that we calculate solvation free energies 
temp_C=linspace(temp_min,temp_max,tempdiv); % temperatures in Celcius            
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if calcflag==1

    % Load the optimized parameters from file
    if ionflag == 1
        optimizedParamFile = 'OptSlicSasaWater_w_ion';
    elseif ionflag == 0
        optimizedParamFile = 'OptSlicSasaWater_wo_ion';
    end
    ParamWatInfo=load(optimizedParamFile);
    x = ParamWatInfo.xvec;
    
    for j=1:tempdiv
        clear global
        loadConstants
        convertKJtoKcal = 1/joulesPerCalorie;
        global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
        saveMemory = 1;
        writeLogfile = 1;
        logfileName = 'junklogfile';
        epsIn  =  1; % dielectric constant of the solutes
        KelvinOffset = 273.15;
        conv_factor = 332.112;
        staticpotential = 0.0; % this only affects charged molecules;
        kappa = 0.0;  % should be zero, meaning non-ionic solutions!
        epsOut = (-1.410e-6)*temp_C(j)^3+(9.398e-4)*temp_C(j)^2-0.40008*temp_C(j)+87.740;  
        % temperature dependence of the dielectric constant of Water T in C
        % from Mollerup15 (Modeling the permittivity of electrolyte solutions)
        
        UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa',kappa,...
                                 'conv_factor',conv_factor,'staticpotential',staticpotential);
        Data = readtable('../ref_data/thermo.csv');
        % Experimental data for 502 neutral molecules and 9 monovalent ions
        % from multiple references - See '../ref_data/thermo_expt.xls' for more details
        % Var1 = solutes | 2 = dG_expt | 3 = Surface_Area | 4 = Volume 
        %    5 = dG_np   | 6 = dG_es   | 7 = dG_Mobley
    
        Data_analogs = readtable('../ref_data/analogs.csv');
        % Experimental data for 12 Amino Acid side-chain analog neutral molecules 
        % See '../ref_data/thermo_expt.xls' for more details
        % Var1 = analogs_list | 2 = dG | 3 = dH | 4 = -TdS | 5 = Cp 
        
        Data_aa = readtable('../ref_data/aa.csv');
        % Experimental data for 12 Amino Acid side-chains from
        % Hess06 - Hydration Thermodynamic Properties of Amino Acid Analogues
        % Var1 = aa_list | 2 = dG | 3 = dH | 4 = -TdS | 5 = Cp 
        
        Data_full = readtable('../ref_data/thermo_full.csv');
        % Full thermodynamic experimental data (dG,dS,dH and Cp) for 111 neutral compounds
        % from multiple references - See '../ref_data/thermo_expt.xls' for more details
        % Var1 = molecules_list | 2 = dG | 3 = dH | 4 = -TdS | 5 = Cp 
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            
        mol_list = Data.Var1;
        if ionflag == 0
            mol_list = Data.Var1(1:502);
        end
        dG_list = Data.Var2;
        surfArea_list = Data.Var3;
        dG_Mobley = Data.Var7(1:502);

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
        [errfinal(j,:),calcE(j,:),refE(j,:),es(j,:),np(j, :)]=ObjectiveFromBEMSA(x(j,:));
    end  
    save('RunSlicSasaWaterThermo.mat','errfinal','calcE','refE','es','np','temp_C','x','mol_list','dG_Mobley');
end

RunThermoResults=load('RunSlicSasaWaterThermo.mat');
dGfunc=struct(); % structure that has the information of the solutes and the linear function that 
                 % fits the calculated values of dG at different temperatures
loadConstants
dG_Mobley=0;
errfinal=RunWater.errfinal;
calcE = RunWater.calcE; % calculated values for dG at different temperatures
refE = RunWater.refE; % calculated values for dG at different temperatures
if strcmp(dataset,'mobley')
    dG_Mobley=RunWater.dG_Mobley;
end
    
es=RunWater.es;
np=RunWater.np;
x = RunWater.x;
TEMP=RunWater.temp_C';
temp_K=TEMP+273.15;
[m,index]=ismember(24.85,TEMP);
mol_list=RunWater.mol_list;

for i=1:length(mol_list)
    f = @(R) (R(1)-R(2)*(temp_K-298)+R(3)*((temp_K-298)-temp_K.*log(temp_K./298)))-calcE(:,i);
    R0=[refE(index,i),1,1];
    options=optimoptions('lsqnonlin','StepTolerance',1e-6);
    options=optimoptions(options,'OptimalityTolerance',1e-6);
    options=optimoptions(options,'FunctionTolerance',1e-6);
    [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);
    dGfunc(i).name=mol_list(i); 
    dGfunc(i).dg=R(1);
    dGfunc(i).ds=R(2);
    dGfunc(i).cp=R(3); 
    dsvec(i)=dGfunc(i).ds*1000;
    cpvec(i)=dGfunc(i).cp*1000;
    resnorm(i)=resnorm;
    exitflag(i)=exitflag;
    output(i)=output;
end

 dg_rms_298_MD=0;
 dg_rms_298=rms(calcE-refE);

 if strcmp(dataset,'mobley')
    dg_rms_298_mol=rms(calcE(index,1:502)-refE(index,1:502));
    dg_rms_298_MD_mol=rms(dG_Mobley'-refE(index,1:502));
    dg_rms_298_ion=rms(calcE(index,503:end)-refE(index,503:end));
 end


if ionflag == 1
    output_name='RunSlicSasaWater_w_ions';
elseif ionflag == 0
    output_name='RunSlicSasaWater_wo_ions';
end

save(output_name,'errfinal','calcE','refE','es','np','TEMP',...
     'x','mol_list','dGfunc','dsvec','cpvec','dg_rms_298_ion',...
     'dg_rms_298_mol','dg_rms_298_MD_mol','index','dG_Mobley',...
     'resnorm','residual','output','exitflag');



    