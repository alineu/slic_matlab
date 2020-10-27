import os
import numpy as np
import pandas as pd
import subprocess
import random
from sklearn.decomposition import PCA
import warnings
import shutil
import csv
from pathlib import Path
from scipy.io import loadmat
import matplotlib.pyplot as plt


# Setting Path Info
homedir = os.environ['HOME']
robustness_repo_path = '%s/repos/slic-robustness-thermo' % homedir
setschenow_repo_path = '%s/repos/setschenow-data' % homedir
dropbox_path = '%s/Dropbox' % homedir
nlbc_path ='%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test' % dropbox_path
msms_path = '/Users/Ali/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1'

# Main Functions

def get_cmap(n, name='jet'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)
    # usage:
    # 
    # cmap = get_cmap(50)
    # ax.plot(x,y,':o',c=cmap(10),mfc='w',markersize=5,label='label1') 
    
def insert(file_name, string):
    
    """
    Inserts "string" ahead of "original_file" and saves as "new_file_name"
    """
    pwd = os.getcwd()
    if os.path.exists(os.path.join(pwd,file_name)):
        
        with open(file_name,'r') as f1:

            with open('temp.txt','w') as f2: 

                f2.write(string)
                f2.write(f1.read())

        f2.close()
        f1.close()

        os.rename('temp.txt',file_name)
    else:
        return "file %s does not exist!" % file_name
    pass

def get_Bondi_radius(atom):
    
    atoms = ["H" , "C" , "N" , "O" , "F" , "P" , "S" , "Cl",\
             "Ar", "As", "Br", "Cd", "Cu", "Ga", "Au", "He",\
             "In", "I" , "Kr", "Pb", "Li", "Mg", "Hg", "Ne",\
             "Ni", "Pd", "Pt", "K" , "Se", "Si", "Ag", "Na",\
             "Te", "Tl", "Sn", "U" , "Xe", "Zn"]

    radii = [1.20, 1.70, 1.55, 1.52, 1.47, 1.80, 1.80, 1.75,\
             1.88, 1.85, 1.85, 1.62, 1.40, 1.87, 1.66, 1.40,\
             1.93, 1.98, 2.02, 2.02, 1.82, 1.73, 1.70, 1.54,\
             1.64, 1.63, 1.80, 2.75, 1.90, 2.10, 1.90, 2.27,\
             2.06, 1.96, 2.17, 1.86, 2.16, 1.37]

    Bondi_radii_dict = dict(zip(atoms, radii))
    
    return Bondi_radii_dict.get(atom,"Atom Not Found!")


def get_Joung_radius(atom):
    
    atoms = ["LI", "NA", "K" , "RB", "CS", "F" , "CL", "BR", "I"] 

    radii = [0.9430, 1.2595, 1.5686, 1.6680, 1.8179, 2.1188, 2.3120, 2.3994, 2.6312]

    Joung_radii_dict = dict(zip(atoms, radii))  
    
    return Joung_radii_dict.get(atom,"Atom Not Found!")


def gen_stern_xyzr(stern_thickness,coeff=1,xyzr_file='test.xyzr',extension='.stern'):

    xyzr_stern = xyzr_file + extension + '%s' % coeff
    
    X = np.loadtxt(xyzr_file)
    X[:,3] += stern_thickness
    X[:,3] *= float(coeff)

    np.savetxt(xyzr_stern, X, fmt='%5.5f')
    
    return xyzr_stern


def write_file_info(file_name,info):
    tmp=open(file_name, "w")
    tmp.write(info)
    tmp.close()
    pass

def get_vert_face_info(filename):
    
    with open('%s.vert'%filename) as f1:

        vert = [line.rstrip('\n').split() for line in f1]

    f1.close()
    with open('%s.face'%filename) as f2:

        face = [line.rstrip('\n').split() for line in f2]

    f2.close()
    vert_data = np.reshape(vert,np.shape(vert)).astype(float)
    face_data = np.reshape(face,np.shape(face)).astype(int)-1
    
    return vert_data,face_data


def get_setschenow_data(salt,setschenowPath = setschenow_repo_path):
    
    a=[]
    salt_csv_path = os.path.join(setschenowPath,'%s.csv'%salt)
    with open(salt_csv_path, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            a.append(row)
    data = np.reshape(a,np.shape(a))
    mols = data[1:,0]
    setschenow_coeffs = data[1:,1]
    
    return mols, setschenow_coeffs

def msms_salt_mol(salt, probe_radius = 1.4, shrink_radius_factor = 1, concentration = 0,
    density = 2,molXYZR = "test.xyzr", path_to_nlbc_data = nlbc_path, msmsPath = msms_path):
    
    curDir = os.getcwd()
    probe_radi_dict = {'NaCl': (get_Joung_radius('NA') + get_Joung_radius('CL'))/2,
                       'KCl' : (get_Joung_radius('K')  + get_Joung_radius('CL'))/2,
                       'NaBr': (get_Joung_radius('NA') + get_Joung_radius('BR'))/2,
                       'LiCl': (get_Joung_radius('LI') + get_Joung_radius('CL'))/2}
    
    stern_thickness = dictionary.get(salt,'Not Found')
    
    for molecule in get_setschenow_data(salt)[0]:
        
        path_to_mol_data = os.path.join(path_to_nlbc_data,molecule) 
        os.chdir(path_to_mol_data)
        salt_diel = "diel_%s_%s_%s"%(molecule,salt,concentration)
        msms_diel = "%s%s%s%s%s%s%s%s%s" % (msmsPath,' -if ',molXYZR,' -no_header -de ',density,
                                            ' -prob ',probe_radius,' -of ',salt_diel)
        molXYZR_stern = gen_stern_xyzr(stern_thickness)
        salt_stern = "stern_%s_%s_%s"%(molecule,salt,concentration)
        msms_stern = "%s%s%s%s%s%s%s%s%s" % (msmsPath,' -if ',molXYZR_stern,' -no_header -de ',density,
                                             ' -prob ',probe_radius,' -of ',salt_stern)
        subprocess.Popen(msms_diel, shell=True).wait()
        subprocess.Popen(msms_stern,shell=True).wait()
        mol_salt_srf = "test_setschenow_%s.srf" % salt
        Path('%s/%s' % (path_to_mol_data,salt_diel)).touch()
        Path('%s/%s' % (path_to_mol_data,salt_stern)).touch()
        write_file_info(mol_salt_srf,"f\nf\n./%s\n./%s\n\n\n0\n" % (salt_stern,salt_diel))
        os.chdir(path_to_nlbc_data)
    os.chdir(curDir)
    pass

def get_SASA_SES_from_msms(
    *,molecule=None,mol_list=None,coef=1,probe_radius=1.4,st_thickness=0,density = 2,
    molXYZR = "test.xyzr",path_to_nlbc_data = nlbc_path,msmsPath=msms_path):
    
    if molecule and mol_list:
        raise TypeError('Not both!')
    if not mol_list: 
        mol_list = [molecule]
    pwd = os.getcwd()
    cp_cmd = 'cp -rf %s/nonpolar/get_surf_info.sh %s' % (setschenow_repo_path,pwd)
    subprocess.Popen(cp_cmd, shell=True).wait()
    chmod_cmd = 'chmod +x ./get_surf_info.sh'
    subprocess.Popen(chmod_cmd, shell=True).wait()
    arr = np.zeros((4,1)).ravel()
    molXYZR_grown = gen_stern_xyzr(st_thickness,coef,molXYZR)
    molXYZR_surf = "surf_%s_%s"%(molecule,coef)
    msms_grow = "%s%s%s%s%s%s%s > msmslog" % (msmsPath,' -if ',molXYZR_grown,' -no_header -de ',density,
                                              ' -prob ',probe_radius,' -of ',molXYZR_surf)
    subprocess.Popen(msms_grow, shell=True).wait()
    info_cmd="./get_surf_info.sh msmslog"
    subprocess.Popen(info_cmd, shell=True).wait()
    if os.stat('surf_data').st_size == 0:
        return arr
    arr = np.loadtxt('surf_data', unpack=True).ravel()
    
    return arr #SES_analytical,SAS_area,SES_volume,SES_area


def get_pqr_from_file(*,file=None,files=None):
    
    if file and files:
        raise TypeError('Not both!')
    if not files: 
        files = [file]
    net_pqr=[]
    
    for i in range(len(files)):
        
        with open(files[i]) as f:

            pqr = [line.rstrip('\n').split() for line in f]

        f.close()
        pqr_new = np.reshape(pqr,np.shape(pqr))[:,-5:].astype(float)
        net_pqr.append(pqr_new)
    pqr=pqr_new
    print(len(net_pqr))
    print(len(files))
    
    if len(net_pqr)==1:
        return pqr
    for j in range(len(net_pqr)-1):
        pqr=np.vstack((pqr,net_pqr[-j-2]))
    return pqr

def get_molecule_geometry(*,molecule=None,mol_list=None):
    
    if molecule and mol_list:
        raise TypeError('Not both!')
    if not mol_list: 
        mol_list = [molecule]
    
    coefs = [1.0]
    df = pd.read_csv('mol_data.csv')
    mollist = df['mol_list'].values
    path_to_nlbc_data ='/Users/Ali/repos/setschenow-data/nlbc_test'
    path_to_np='/Users/Ali/repos/setschenow-data/nonpolar'
    n=len(mollist)
    data = np.zeros((n,len(coefs),4))
    for i in range(n):
        path_to_mol_data = os.path.join(path_to_nlbc_data,mollist[i])
        os.chdir(path_to_mol_data)
        print(path_to_mol_data)
        for j in range(len(coefs)):
            try:
                data[i,j]=get_SASA_from_msms(mollist[i],coefs[j],probe_radius=0.1)
            except ValueError:
                data[i,j]=[-1,-1,-1,-1]
                pass
    os.chdir(path_to_np)
    np.savetxt(os.path.join(path_to_np,'data_vdw.txt'),data.reshape(n,4),fmt='%4.2f')

def mol_pca(fileName):
    
    line_count = pd.read_csv("%s/test.pqr" % fileName,header=None).values.shape[0]
    pos=[]
    pos.append([pd.read_csv("%s/test.pqr" % fileName,header=None).values[i][0].split()[5:8] for i in range(line_count)])
    rad=[]
    rad.append([pd.read_csv("%s/test.pqr" % fileName,header=None).values[i][0].split()[-1] for i in range(line_count)])
    if float(rad[0][0])>0.5:
        mean_r = np.mean(np.array(rad).astype(float))
    else:
        mean_r=1.5
    pos = np.array(pos).astype(float).reshape((line_count,3))
    pca = PCA()
    pca.fit(pos)
    u=pca.explained_variance_+2*mean_r
    norm_u = u/np.sqrt(u@u)
    return pos, pca.explained_variance_,sorted(u),sorted(norm_u)

    
def grow_surf(coeff,xyzr_diel='test.xyzr'):

    xyzr_stern = xyzr_diel + '%3.2f' % coeff

    X = np.loadtxt(xyzr_diel)

    X[:,3] *= coeff

    np.savetxt(xyzr_stern, X, fmt='%5.5f')
    return xyzr_stern

def write_matlab_setschenow_file(
    solvent,salt,concentration,ref_number=4,training_set=['ethanol'],
    x0=[1.23,-15.34,0.63,0.49,2.820,0.003,1.640],run_type='param',
    max_iter=8,csv_path='%s/data' % setschenow_repo_path):
    
    run_dict = {'param':[0,0,'testset'],'run':[1,1,'mol_list']}

    info = """%% Path information
close all;
clear all;
Home = getenv('HOME');
addpath(sprintf('%%s/repos/pointbem',Home));
addpath(sprintf('%%s/repos/panelbem',Home));
addpath(sprintf('%%s/repos/testasymmetry',Home));
addpath(sprintf('%%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%%s/repos/testasymmetry/born',Home));
addpath(sprintf('%%s/repos/setschenow-data',Home));

salt = '%s';
concentration = %s; %% in molar
kappa = sqrt(3*concentration)/3.047; 
epsOut = calculateDielectricConstant(salt, concentration);
fid = fopen('setschenow.csv','r');
Data = textscan(fid,'%%s %%f  %%f','delimiter',',');
fclose(fid);
mol_list = Data{1};
dG_list = Data{2};
surfArea_list = Data{3};
setschenowDatafile = sprintf('%%s.csv',salt);

loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = %d;
writeLogfile = %d;
logfileName = 'junklogfile';
epsIn  =  1;
Tbase = 300; 
KelvinOffset = 273.15;
conv_factor = 332.112;
activityDatafile = sprintf('activity_coefficients/%%s%%i.csv',salt,%d);
activityFile = fopen(activityDatafile,'r');
C = textscan(activityFile,'%%f %%f %%f %%s %%f %%f','HeaderLines',1,'Delimiter',',');
fclose(activityFile);


densityDatafile = sprintf('salt_densities/%%s.csv',salt);
densityFile = fopen(densityDatafile,'r');
D = textscan(densityFile,'%%f %%f %%f %%f %%*[^\\n]','HeaderLines',1,'Delimiter',',');
fclose(densityFile);


[setschenow_mol_list, setschenow_coefficients] = textread(setschenowDatafile,'%%s %%f','delimiter',',','headerlines',1);
epsOut = calculateDielectricConstant(salt, concentration);

staticpotential = 0.0; 
molal_conc = interp1(D{1, 2},D{1, 4},concentration);

ddg_K = interp1(C{1, 1},C{1, 5},molal_conc);
ddg_Cl = interp1(C{1, 1},C{1, 6},molal_conc);
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
kappa,'conv_factor',conv_factor,...
'staticpotential',staticpotential);

testset  = {'%s'};
curdir=pwd;

for i=1:length(%s)
  dir=sprintf('%%s/repos/setschenow-data/nlbc_test/%%s',getenv('HOME'),%s{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%%s/test_setschenow_%%s.srf',dir,salt);
  chargeDist{i} = pqrData.q;\n""" % (salt,concentration,run_dict.get(run_type,"not Found!")[0],\
                                     run_dict.get(run_type,"not Found!")[1],ref_number,\
                                     "','".join(training_set),\
                                     run_dict.get(run_type,"not Found!")[2],\
                                     run_dict.get(run_type,"not Found!")[2])
    
    info_end_param="""  foo = strcmp(mol_list,%s{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\\n');
    keyboard
  end
  referenceData{i} = dG_list(index);
  surfArea{i} = surfArea_list(index);
  chdir(curdir);
  addProblemSA(%s{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end
x0 = [%s];
lb = [0 -200 -100 -1  -0.1  -0.1  -4];
ub = [+2 +200 +100 +1  +0.1  +0.1  +4];

options = optimoptions('lsqnonlin','MaxIter',%d);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);
rms_train = rms(ref-calc);
solvent = '%s';
Optfile = sprintf('Opt%%s_%%s_%%s.mat',solvent,salt,concentration);
save(Optfile,'x','ref','calc','es','np','x0','calc0','es0','np0','rms_train');""" % (run_dict.get(run_type,"not Found!")[2],\
                                                                                     run_dict.get(run_type,"not Found!")[2],\
                                                                                     ' '.join(str("{0:.3f}".format(e)) for e in x0),\
                                                                                     max_iter,solvent)

    info_end_run = """  referenceData{i} = dG_list(i);
  surfArea{i} = surfArea_list(i);
  chdir(curdir);
  addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end
solvent = '%s';
Optfile = sprintf('Opt%%s_%%s_%%s.mat',solvent,salt,concentration);
ParamInfo = load(Optfile);
x = ParamInfo.x;
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x);
rms_test = rms(refE-calcE);
Runfile = sprintf('Run%%s_%%s_%%s.mat',solvent,salt,concentration);
save(Runfile,'x','errfinal','calcE','refE','es','np','rms_test');""" % (solvent)

    if run_type == 'run':
        info += info_end_run
    else:
        info += info_end_param
    conc = "%3.1f" % concentration
    conc_str = "%s" % conc
    conc_s = "%s%s"%(conc_str[-3],conc_str[-1])
    tmp=open("%s%s_%s_%s.m" % (run_type,solvent,salt,conc_s), "w")
    tmp.write(info)
    tmp.close()
    pass


def get_solute_family_data(solvent,ref_csv = '%s/data/mobley_all.csv' % robustness_repo_path,family_limit=10):
    
    ref_df = pd.read_csv(ref_csv,names=["solute","surf_area","volume","dG_expt","dG_Mob","dG_es","dG_np","family"])
    solvent_df = read_solvent(solvent)
    solvent_solutes = solvent_df['solute'].values
    solute_families = ref_df.family[solvent_solutes].values
    solute_families_unique = set(solute_families)
    len_families = len(solute_families_unique)
    family_dict = dict(zip(solute_families_unique,np.arange(0,len_families)))
    color_index = [family_dict.get(x) for x in solute_families]
    solvent_df['solute_col_id'] = pd.Series(color_index, index=solvent_df.index)
    solvent_df['family'] = pd.Series(ref_df.family[list(solvent_solutes)].values,index=solvent_df.index)
    
    return solvent_df,len_families


def plot_solvent_results(
    solvent_family,solvent,training_data='test_consistent',tr_set_len=14,
    fig_size=(14, 8),plot_tr = True,classify_solutes=False,min_len_solute_class=3):
    
    solvent_fam_data_path = os.path.join(robustness_repo_path,"data/%s" % solvent_family)
    if training_data=="test_all":
        plot_data_path = os.path.join(solvent_fam_data_path,'test_all')
    elif training_data=="test_consistent":
        tr_set_folder_name = "test_%s_%d" % (solvent_family,tr_set_len)
        plot_data_path = os.path.join(solvent_fam_data_path,\
                                      tr_set_folder_name,\
                                      training_data)
    else:
        rand_set_num = int(training_data[-1])
        tr_set_folder_name = "test_%s_%d" % (solvent_family,tr_set_len)
        if rand_set_num!=0:
            training_data = "test_random_%d"% (rand_set_num-1)
        else:
            training_data = "test_random_9"
        plot_data_path = os.path.join(solvent_fam_data_path,\
                                      tr_set_folder_name,\
                                      training_data)
    trainig_MatFile = "%s/Opt%s.mat" % (plot_data_path,solvent)
    test_MatFile = "%s/Run%s.mat" % (plot_data_path,solvent)
    
    x_tr = loadmat(trainig_MatFile).get('x')[0]
    ref_tr_vec = loadmat(trainig_MatFile).get('ref').ravel()
    calc_tr_vec = loadmat(trainig_MatFile).get('calc').ravel()
    rms_tr = rmse(calc_tr_vec,ref_tr_vec)
    ref_test_vec = loadmat(test_MatFile).get('refE').ravel()
    calc_test_vec = loadmat(test_MatFile).get('calcE').ravel()
    rms_test = rmse(calc_test_vec,ref_test_vec)
    f, ax = plt.subplots(figsize=fig_size,dpi=250)
    if plot_tr:
        ax.scatter(ref_tr_vec,calc_tr_vec, s=50, marker='d',\
                   facecolors='r', edgecolors='b',linewidth='1',label="RMSE_train = %4.2f" % rms_tr,zorder=10)
    if not classify_solutes:
        ax.scatter(ref_test_vec,calc_test_vec, s=100, marker='o',\
                   facecolors='none', edgecolors='r',label="RMSE_test = %4.2f" % rms_test,zorder=9)
    else:
        df, len_families = get_solute_family_data(solvent)
        cm = get_cmap(len(np.unique(df['family'].values))*100,'jet')
        print(len(np.unique(df['family'].values)))
        df['refE'] = pd.Series(ref_test_vec, index=df.index)
        df['calcE'] = pd.Series(calc_test_vec, index=df.index)
        colored_df = pd.concat(g for _, g in df.groupby("solute_col_id") if len(g) > min_len_solute_class)
        family_ids = np.unique(colored_df['family'].values)
        uncolored_df = pd.concat(g for _, g in df.groupby("solute_col_id") if len(g) <= min_len_solute_class)
        colorless_id = uncolored_df.solute_col_id.iloc[0]
        for g in family_ids:
            ix = df['family'].values == g
            color_id = df[df.solute_col_id.where(ix).notnull()].solute_col_id.iloc[0]
            ax.scatter(df['refE'].where(ix),df['calcE'].where(ix),facecolors=cm(color_id*80+1),label = g, s = 100)
        ax.scatter(uncolored_df['refE'],uncolored_df['calcE'],facecolors=cm(len(np.unique(df['family'].values))*100),label = 'other', s = 100)
    ax.set_xlabel("$\Delta\,G_{\,expt}^{\,solv}$ " + solvent,fontsize=16)
    ax.set_ylabel("$\Delta\,G_{\,calc}^{\,solv}$ " + solvent,fontsize=16)
    ax.set_title("Experimental vs. predicted solvation free energies of "+\
                 solvent+" using "+ training_data + " with training-set length "+\
                 "%d"%tr_set_len,fontsize=16)
    ax.plot(ref_test_vec,ref_test_vec,'k-') # refline
    ax.legend()
    
    return ax

def plot_water_es(solvent_family,training_data='test_consistent',ref_csv = '%s/data/mobley_all.csv' % robustness_repo_path,
    tr_set_len=14,fig_size=(14, 8),plot_tr = True,
    classify_solutes=False,min_len_solute_class=20,hide_other_solutes=False):
    
    solvent='water'
    ref_df = pd.read_csv(ref_csv,names=["solute", "surf_area","volume","dG_expt","dG_Mob","dG_es","dG_np","family"])
    ref_df.set_index("solute",inplace=True)
    solvent_fam_data_path = os.path.join(robustness_repo_path,"data/%s" % solvent_family)
    if training_data=="test_all":
        plot_data_path = os.path.join(solvent_fam_data_path,'test_all')
    elif training_data=="test_consistent":
        tr_set_folder_name = "test_%s_%d" % (solvent_family,tr_set_len)
        plot_data_path = os.path.join(solvent_fam_data_path,\
                                      tr_set_folder_name,\
                                      training_data)
    else:
        rand_set_num = int(training_data[-1])
        tr_set_folder_name = "test_%s_%d" % (solvent_family,tr_set_len)
        if rand_set_num!=0:
            training_data = "test_random_%d"% (rand_set_num-1)
        else:
            training_data = "test_random_9"
        plot_data_path = os.path.join(solvent_fam_data_path,\
                                      tr_set_folder_name,\
                                      training_data)
    trainig_MatFile = "%s/Opt%s.mat" % (plot_data_path,'water')
    test_MatFile = "%s/Run%s.mat" % (plot_data_path,'water')
    
    ref_test_es = ref_df.dG_es.values[:502]
    calc_test_es = loadmat(test_MatFile).get('es').ravel()
    rms_test = rmse(calc_test_es,ref_test_es)
    f, ax = plt.subplots(figsize=fig_size,dpi=250)
    if not classify_solutes:
        ax.scatter(ref_test_es,calc_test_es, s=100, marker='o',\
               facecolors='none', edgecolors='r',label="RMSE_test = %4.2f" % rms_test)
        ax.set_title("MD vs. predicted electrostatic solvation free energies of "+\
                 "water using "+ training_data + " with training-set length "+\
                 "%d"%tr_set_len,fontsize=16)
    else:
        df, len_families = get_solute_family_data('water')
        cm = get_cmap(len(np.unique(df['family'].values))*100,'jet')
        print(len(np.unique(df['family'].values)))
        df['refEs'] = pd.Series(ref_test_es, index=df.index)
        df['calcEs'] = pd.Series(calc_test_es, index=df.index)
        colored_df = pd.concat(g for _, g in df.groupby("solute_col_id") if len(g) > min_len_solute_class)
        family_ids = np.unique(colored_df['family'].values)
        uncolored_df = pd.concat(g for _, g in df.groupby("solute_col_id") if len(g) <= min_len_solute_class)
        colorless_id = uncolored_df.solute_col_id.iloc[0]
        for g in family_ids:
            ix = df['family'].values == g
            color_id = df[df.solute_col_id.where(ix).notnull()].solute_col_id.iloc[0]
            ax.scatter(df['refEs'].where(ix),df['calcEs'].where(ix),facecolors=cm(color_id*80+1),label = g, s = 100)
        if not hide_other_solutes:
            ax.scatter(uncolored_df['refEs'],uncolored_df['calcEs'],facecolors=cm(len(np.unique(df['family'].values))*100),label="other families", s = 40,alpha=0.2)
        ax.set_title("MD vs. predicted electrostatic solvation free energies of "+\
                 "water using "+ training_data + " with training-set length "+\
                 "%d - RMSE = %4.2f kcal/mol"%(tr_set_len,rms_test),fontsize=14)
    
    ax.set_xlabel("$\Delta\,G_{\,es,MD}^{\,solv}$ " + solvent,fontsize=16)
    ax.set_ylabel("$\Delta\,G_{\,es,calc}^{\,solv}$ " + solvent,fontsize=16)
    ax.plot(ref_test_es,ref_test_es,'k-') # refline
    ax.legend()
    
    return ax

def rmse(calc, ref, decimal_points=3):
    return np.round(np.sqrt(((calc - ref) ** 2).mean()),decimal_points)

def read_solvent(
    solvent,path_to_solvent = '%s/data' % robustness_repo_path):
    
    solvent_csv_file = "%s.csv" % solvent
    file_path = os.path.join(path_to_solvent,solvent_csv_file)
    df = pd.read_csv(file_path,names=["solute", "dG"])
    return df
    
def read_solvent_family(
    family,path_to_solvent= '%s/data' % robustness_repo_path):
    family_csv_file = "%s.csv" % family
    file_path = os.path.join(path_to_solvent,family_csv_file)
    df = pd.read_csv(file_path,names=["solute"])    
    return df

def generate_surf_area_csv(
    solvent,path_to_solvent= '%s/data' % robustness_repo_path,surf_area_file = 'mobley_all.csv'):
    
    file_path = os.path.join(path_to_solvent,surf_area_file)
    surf_df = pd.read_csv(file_path,names=["solute", "surf_area","volume","dG_expt","dG_Mob","dG_es","dG_np","family"])
    surf_df.set_index("solute",inplace=True)
    df_solvent = read_solvent(solvent)
    solvent_surf_df = surf_df.loc[list(df_solvent.solute.values)]
    solvent_surf_df.to_csv("%s_surf.csv" %solvent, header=False)#, index=False)
    pass
                           
def generate_random_training_set(
    solvent_family,training_set_length = 8, path_to_solvent= '%s/data' % robustness_repo_path):
    
    df = read_solvent_family(solvent_family)
    
    return random.sample(list(df.solute.values),training_set_length)


def write_matlab_file(
    solvent,eps_out,training_set=['ethanol'],x0=[0.5,-60,-0.5,-0.5*np.tanh(0.5),0,0,0],
    run_type='param',max_iter=8,csv_path='%s/data' % robustness_repo_path):
    
    run_dict = {'param':[0,0,'testset'],'run':[1,1,'mol_list']}

    info = """%% Path information
close all;
clear all;
Home = getenv('HOME');
addpath(sprintf('%%s/repos/pointbem',Home));
addpath(sprintf('%%s/repos/panelbem',Home));
addpath(sprintf('%%s/repos/testasymmetry',Home));
addpath(sprintf('%%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%%s/repos/testasymmetry/born',Home));
addpath(sprintf('%%s/repos/slic-robustness-thermo/data',Home));
addpath(sprintf('%%s/repos/setschenow-data',Home));

%% a bunch of useful variables and constants. also defining the global
%% variable \"ProblemSet\" which well use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = %d;
writeLogfile = %d;
logfileName = 'junklogfile';
epsIn  =  1;
Tbase = 300; 
epsOut = %5.2f; %% from MNSol
mytemp=Tbase;
KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; %% this only affects charged molecules;
kappa = 0.0;  %% should be zero, meaning non-ionic solutions!


%% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
kappa,'conv_factor',conv_factor,...
'staticpotential',staticpotential);
solvent = '%s';
solvent_csv=sprintf('data/%%s.csv',solvent);
solvent_surf_csv=sprintf('%%s_surf.csv',solvent);
fid = fopen(solvent_csv,'r'); 
Data = textscan(fid,'%%s %%f  %%f','delimiter',',');
fclose(fid);
mol_list = Data{1};
dG_list = Data{2};
fid_surf = fopen(solvent_surf_csv,'r'); 
Data_surf = textscan(fid_surf,'%%s %%f  %%f','delimiter',',');
fclose(fid_surf);
surfArea_list = Data_surf{2};
vol_list = Data_surf{3};
testset  = {'%s'};
curdir=pwd;

for i=1:length(%s)
  dir=sprintf('%%s/Dropbox/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%%s',getenv('HOME'),%s{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;%%chargeDistribution;\n""" % (run_dict.get(run_type,"not Found!")[0],\
                                                          run_dict.get(run_type,"not Found!")[1],\
                                                          eps_out,solvent,"','".join(training_set),\
                                                          run_dict.get(run_type,"not Found!")[2],\
                                                          run_dict.get(run_type,"not Found!")[2])

    info_end_param="""  foo = strcmp(mol_list,testset{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\\n');
    keyboard
  end
  referenceData{i} = dG_list(index);
  surfArea{i} = surfArea_list(index);
  chdir(curdir);
  addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end
x0 = [%s];
lb = [0 -200 -100 -1  -0.1  -0.1  -4];
ub = [+2 +200 +100 +1  +0.1  +0.1  +4];

options = optimoptions('lsqnonlin','MaxIter',%d);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);
rms_train = rms(ref-calc);
save('Opt%s.mat','x','ref','calc','es','np','x0','calc0','es0','np0','rms_train');""" % (' '.join(str("{0:.3f}".format(e)) for e in x0),max_iter,solvent)

    info_end_run = """  referenceData{i} = dG_list(i);
  surfArea{i} = surfArea_list(i);
  chdir(curdir);
  addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end
ParamInfo = load('Opt%s.mat');
x = ParamInfo.x;
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x);
rms_test = rms(refE-calcE);
save('Run%s.mat','x','errfinal','calcE','refE','es','np','rms_test');""" % (solvent,solvent)

    if run_type == 'run':
        info += info_end_run
    else:
        info += info_end_param
        
    tmp=open("%s%s.m" % (run_type,solvent), "w")
    tmp.write(info)
    tmp.close()
    pass

def run_matlab_python(solvents,run_type='param',repo_path=robustness_repo_path):
    solvent_list = []
    for solvent in solvents:
        solvent_list.append("%s%s"%(run_type,solvent))
    py_script = """import matlab.engine
import os
Home = os.environ['HOME']
Solvent_List = ['%s']

for solvent in Solvent_List:
    eng = matlab.engine.start_matlab()
    eng.eval(solvent,nargout=0)
    eng.quit()""" % ("','".join(solvent_list))
    
    # map the inputs to the function blocks
    tmp=open("%s_solvents.py" % run_type, "w")
    tmp.write(py_script)
    tmp.close() 
    run_cmd = '/usr/bin/python3 %s_solvents.py' % run_type
    subprocess.Popen(run_cmd, shell=True).wait()
    pass

def get_solute_family_data(solvent,ref_csv = '%s/data/mobley_all.csv' % robustness_repo_path,family_limit=10):
    
    ref_df = pd.read_csv(ref_csv,names=["solute","surf_area","volume","dG_expt","es","np","family"])
    solvent_df = read_solvent(solvent)
    solvent_solutes = solvent_df['solute'].values
    solute_families = ref_df.family[solvent_solutes].values
    solute_families_unique = set(solute_families)
    len_families = len(solute_families_unique)
    family_dict = dict(zip(solute_families_unique,np.arange(0,len_families)))
    color_index = [family_dict.get(x) for x in solute_families]
    solvent_df['solute_col_id'] = pd.Series(color_index, index=solvent_df.index)
    solvent_df['family'] = pd.Series(ref_df.family[list(solvent_solutes)].values,index=solvent_df.index)
    
    return solvent_df,len_families

def get_mol_data_from_pqr(file='test.pqr'):
    
    with open(file) as f:

        pqr = [line.rstrip('\n').split() for line in f]

    f.close()
    pqr_new = np.reshape(pqr,np.shape(pqr))
    pqr_str = pqr_new[:,:5]
    pqr_num = pqr_new[:,-5:].astype(float)
    pqr_f = np.hstack((pqr_str,pqr_num))
    return pqr_f

def write_param_files(solvent_family,tr_set_lens = [8,11,14],\
                      num_rand_tests=10,read_tr_set = True):


    for training_set_length in tr_set_lens:
        solvent_list, all_set, consistent_training_set = solvent_data_dict.get((training_set_length,solvent_family),"Not Found!!")
        test_folder_name = "test_%s_%d" % (solvent_family,training_set_length)
        common_molecules = "molecules_%s" % solvent_family

        test_folder_path = os.path.join(repo_path,test_folder_name)
        if os.path.exists(test_folder_path):

            shutil.rmtree(test_folder_path) 

        os.mkdir(test_folder_path)
        os.chdir(test_folder_path)

        # Param consistent training-set
        consistent_path = os.path.join(test_folder_path,'test_consistent')
        os.mkdir(consistent_path)
        os.chdir(consistent_path)
        np.savetxt('consistent_training_set.txt',consistent_training_set, fmt='%s')
        for solvent in solvent_list:
            eps_index = list(eps_df.solvent.values).index(solvent)
            eps_solvent = eps_df.eps.values[eps_index]
            generate_surf_area_csv(solvent)
            write_matlab_file(solvent,eps_solvent,consistent_training_set)
        # run_matlab_python(solvent_list)
        os.chdir(test_folder_path)

        # Param All (all common solutes in the training-set)
        all_set_path = os.path.join(test_folder_path,'test_all')
        os.mkdir(all_set_path)
        os.chdir(all_set_path)
        np.savetxt('all_training_set.txt',all_set, fmt='%s')
        for solvent in solvent_list:
            eps_index = list(eps_df.solvent.values).index(solvent)
            eps_solvent = eps_df.eps.values[eps_index]
            generate_surf_area_csv(solvent)
            write_matlab_file(solvent,eps_solvent,all_set)
        # run_matlab_python(solvent_list)
        os.chdir(test_folder_path)

        # Param random training-sets
        training_sets = []
        for i in range(num_rand_tests):
            test_path = os.path.join(test_folder_path,'test_random_%d' % i)
            os.mkdir(test_path)
            os.chdir(test_path)
            if not read_tr_set:
                tr_set = generate_random_training_set(common_molecules,\
                         training_set_length=len(consistent_training_set))
            else:
                tr_set_file = '%s/training_set_%d.txt' % (test_path,i)
                tr_set = [line.rstrip('\n') for line in open(tr_set_file)]
            training_sets.append(tr_set)
            np.savetxt('training_set_%d.txt'%i,training_sets[i], fmt='%s')

            for solvent in solvent_list[:]:

                eps_index = list(eps_df.solvent.values).index(solvent)
                eps_solvent = eps_df.eps.values[eps_index]
                generate_surf_area_csv(solvent)

                write_matlab_file(solvent,eps_solvent,training_sets[i],run_type='param')

            # run_matlab_python(solvent_list)
        os.chdir(test_folder_path)
    pass



# Data for dielectric constants of different solvents (from MNsol)

eps_csv_file = "eps.csv"
file_path = os.path.join('%s/data'%robustness_repo_path,eps_csv_file)
eps_df = pd.read_csv(file_path,names=["solvent","eps"]) 

# Solvent Families

alkanes_solvent_list = ['pentane','hexane','heptane','octane','nonane','decane']
alcohols_solvent_list = ['methanol','ethanol','propanol','butanol','octanol','water']
multiple_solvent_list = ['benzene','hexadecane','chloroform','cyclohexane','carbontet','hexane','toluene','xylene','water']
solvent_families = ['alkanes','alcohols','multiple']
all_solvents = ['2methylpyridine','4methyl2pentanone','aceticacid','acetonitrile','acetonitrile_ions',
                'acetophenone','aniline','benzene','benzonitrile','benzylalcohol','bromobenzene',
                'bromoethane','bromoform','bromooctane','butanol','butanol_ions','butanone','butylacetate',
                'butylbenzene','carbondisulfide','carbontet','chlorobenzene','chloroform','chlorohexane',
                'cyclohexane','cyclohexanone','decalin','decane','decanol','dibromoethane','dibutylether',
                'dichloroethane','dichloroethane_ions','dichloromethane_only_ions','diethylether',
                'diisopropylether','dimethylacetamide_ions','dimethylformamide','dimethylformamide_ions',
                'dimethylpyridine','dimethylsulfoxide','dimethylsulfoxide_ions','dodecane','ethanol',
                'ethanol_ions','ethoxybenzene','ethylacetate','ethylbenzene','fluorobenzene','fluoroctane',
                'heptane','heptanol','hexadecane','hexadecyliodide','hexamethylphosphoramide_only_ions',
                'hexane','hexanol','iodobenzene','ions_alcohols','ions_polar','isobutanol','isooctane',
                'isopropanol','isopropylbenzene','mcresol','mesitylene','methanol','methanol_ions',
                'methoxyethanol','methylenechloride','methylformamide','molecules_alcohols',
                'molecules_alkanes','molecules_multi','nitrobenzene_ions','nitroethane','nitromethane_ions',
                'nonane','nonanol','octane','octanol','octanol_ions','odichlorobenzene','onitrotoluene',
                'pentadecane','pentane','pentanol','perfluorobenzene','phenylether','propanol','propanol_ions',
                'propanone','propanone_ions','propylenecarbonate_only_ions','secbutanol','secbutylbenzene',
                'tbutylbenzene','tetrachloroethene','tetrahydrofuran','tetrahydrothiophenedioxide','tetralin',
                'tetramethylenesulphone_only_ions','toluene','tributylphosphate','triethylamine',
                'trimethylbenzene','undecane','water','water_ions','xylene']

# All common solutes in each family

alkanes_all_set = list(read_solvent_family("molecules_alkanes").solute.values)
alcohols_all_set = list(read_solvent_family("molecules_alcohols").solute.values)
multiple_all_set = list(read_solvent_family("molecules_multiple").solute.values)

# Training sets of different lengths


alkanes_consistent_training_set_8 = ["33_dimethylbutan_2_one","ethanol","heptan_2_one",
                                      "heptan_1_ol","hexan_2_one","n_propylamine",
                                      "methyl_hexanoate","pentan_1_ol"]

alcohols_consistent_training_set_8 = ["anthracene","4_nitroaniline","pyrene",
                                      "n_hexane","n_heptane","toluene",
                                      "ethane","benzene"]

multiple_consistent_training_set_8 = ["acetic_acid","ethanol","ethylamine",
                                      "n_butyl_acetate","n_octane","nitromethane",
                                      "p_cresol","pyridine"]

alkanes_consistent_training_set_11 = ["33_dimethylbutan_2_one","ethanol","heptan_2_one","heptan_1_ol",
                                      "n_butyl_acetate","methanol","ethylamine","hexan_2_one",
                                      "n_propylamine","methyl_hexanoate","pentan_1_ol"]

alcohols_consistent_training_set_11 = ["ethanol","propan_1_ol","butanone","anthracene",
                                       "4_nitroaniline","pyrene","n_hexane","n_heptane",
                                       "toluene","ethane","benzene"]

multiple_consistent_training_set_11 = ["acetic_acid","ethanol","ethylamine","n_butyl_acetate",
                                       "n_octane","nitromethane","p_cresol","pyridine",
                                       "methyl_pentanoate","butan_1_ol","propanone"]

alkanes_consistent_training_set_14 = ["33_dimethylbutan_2_one","butan_1_ol","ethanol","ethyl_acetate","ethylamine",
                                      "heptan_1_ol","heptan_2_one","hexan_2_one","methanol","methyl_hexanoate",
                                      "methyl_propanoate","n_butyl_acetate","n_propylamine","pentan_1_ol"]

alcohols_consistent_training_set_14 =["4_nitroaniline","anthracene","benzene","butanone","cyclohexane",
                                      "methanol","n_butane","n_hexane","n_octane","nitromethane",
                                      "propan_1_ol","pyrene","toluene","triethylamine"]

multiple_consistent_training_set_14 =["14_dioxane","ethanol","ethyl_acetate","ethylamine","hexan_1_ol",
                                      "methyl_pentanoate","n_octane","nitromethane","o_cresol","phenol",
                                      "propan_1_ol","propanoic_acid","propanone","toluene"]


solvent_data_dict = {}
solvent_data_dict[(8,"alkanes")] = [alkanes_solvent_list,
                                    alkanes_all_set,
                                    alkanes_consistent_training_set_8]

solvent_data_dict[(11,"alkanes")] = [alkanes_solvent_list,
                                     alkanes_all_set,
                                     alkanes_consistent_training_set_11]

solvent_data_dict[(14,"alkanes")] = [alkanes_solvent_list,
                                     alkanes_all_set,
                                     alkanes_consistent_training_set_14]

solvent_data_dict[(8,"alcohols")] = [alcohols_solvent_list,
                                     alcohols_all_set,
                                     alcohols_consistent_training_set_8]

solvent_data_dict[(11,"alcohols")] = [alcohols_solvent_list,
                                      alcohols_all_set,
                                      alcohols_consistent_training_set_11]

solvent_data_dict[(14,"alcohols")] = [alcohols_solvent_list,
                                      alcohols_all_set,
                                      alcohols_consistent_training_set_14]

solvent_data_dict[(8,"multiple")] = [multiple_solvent_list,
                                     multiple_all_set,
                                     multiple_consistent_training_set_8]

solvent_data_dict[(11,"multiple")] = [multiple_solvent_list,
                                      multiple_all_set,
                                      multiple_consistent_training_set_11]

solvent_data_dict[(14,"multiple")] = [multiple_solvent_list,
                                      multiple_all_set,
                                      multiple_consistent_training_set_14]
