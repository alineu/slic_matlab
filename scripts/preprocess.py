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

# Set path
homedir = os.environ['HOME']
repo_path = '%s/repos/slic_matlab' % homedir
nlbc_path ='%s/ref_data/nlbc_test' % repo_path
msms_path = '%s/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1' % homedir

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
    if os.path.exists(os.path.join(pwd, file_name)):

        with open(file_name, 'r') as f1:

            with open('temp.txt', 'w') as f2:

                f2.write(string)
                f2.write(f1.read())

        f2.close()
        f1.close()

        os.rename('temp.txt', file_name)
    else:
        return "file %s does not exist!" % file_name
    pass


def get_Bondi_radius(atom):
    """
    Bondi Radii for different atoms
    """
    atoms = ["H", "C", "N", "O", "F", "P", "S", "Cl",
             "Ar", "As", "Br", "Cd", "Cu", "Ga", "Au", "He",
             "In", "I", "Kr", "Pb", "Li", "Mg", "Hg", "Ne",
             "Ni", "Pd", "Pt", "K", "Se", "Si", "Ag", "Na",
             "Te", "Tl", "Sn", "U", "Xe", "Zn"]

    radii = [1.20, 1.70, 1.55, 1.52, 1.47, 1.80, 1.80, 1.75,
             1.88, 1.85, 1.85, 1.62, 1.40, 1.87, 1.66, 1.40,
             1.93, 1.98, 2.02, 2.02, 1.82, 1.73, 1.70, 1.54,
             1.64, 1.63, 1.80, 2.75, 1.90, 2.10, 1.90, 2.27,
             2.06, 1.96, 2.17, 1.86, 2.16, 1.37]

    Bondi_radii_dict = dict(zip(atoms, radii))

    return Bondi_radii_dict.get(atom, "Atom Not Found!")


def get_Joung_radius(atom):
    """
    Joung-Cheatham Radii for different atoms
    """
    atoms = ["LI", "NA", "K", "RB", "CS", "F", "CL", "BR", "I"]

    radii = [0.9430, 1.2595, 1.5686, 1.6680,
             1.8179, 2.1188, 2.3120, 2.3994, 2.6312]

    Joung_radii_dict = dict(zip(atoms, radii))

    return Joung_radii_dict.get(atom, "Atom Not Found!")


def gen_stern_xyzr(stern_thickness, coeff=1, xyzr_file='test.xyzr', extension='.stern'):

    xyzr_stern = xyzr_file + extension + '%s' % coeff

    X = np.loadtxt(xyzr_file)
    X[:, 3] += stern_thickness
    X[:, 3] *= float(coeff)

    np.savetxt(xyzr_stern, X, fmt='%5.5f')

    return xyzr_stern


def write_file_info(file_name, info):
    tmp = open(file_name, "w")
    tmp.write(info)
    tmp.close()
    pass


def get_vert_face_info(filename):

    with open('%s.vert' % filename) as f1:

        vert = [line.rstrip('\n').split() for line in f1]

    f1.close()
    with open('%s.face' % filename) as f2:

        face = [line.rstrip('\n').split() for line in f2]

    f2.close()
    vert_data = np.reshape(vert, np.shape(vert)).astype(float)
    face_data = np.reshape(face, np.shape(face)).astype(int) - 1

    return vert_data, face_data


def get_SASA_SES_from_msms(
        *, molecule=None, mol_list=None, coef=1, probe_radius=1.4, st_thickness=0, density=2,
        molXYZR="test.xyzr", path_to_nlbc_data=nlbc_path, msmsPath=msms_path):

    if molecule and mol_list:
        raise TypeError('Not both!')
    if not mol_list:
        mol_list = [molecule]
    pwd = os.getcwd()
    cp_cmd = 'cp -rf %s/nonpolar/get_surf_info.sh %s' % (
        setschenow_repo_path, pwd)
    subprocess.Popen(cp_cmd, shell=True).wait()
    chmod_cmd = 'chmod +x ./get_surf_info.sh'
    subprocess.Popen(chmod_cmd, shell=True).wait()
    arr = np.zeros((4, 1)).ravel()
    molXYZR_grown = gen_stern_xyzr(st_thickness, coef, molXYZR)
    molXYZR_surf = "surf_%s_%s" % (molecule, coef)
    msms_grow = "%s%s%s%s%s%s%s > msmslog" % (msmsPath, ' -if ', molXYZR_grown, ' -no_header -de ', density,
                                              ' -prob ', probe_radius, ' -of ', molXYZR_surf)
    subprocess.Popen(msms_grow, shell=True).wait()
    info_cmd = "./get_surf_info.sh msmslog"
    subprocess.Popen(info_cmd, shell=True).wait()
    if os.stat('surf_data').st_size == 0:
        return arr
    arr = np.loadtxt('surf_data', unpack=True).ravel()

    return arr  # SES_analytical,SAS_area,SES_volume,SES_area


def get_pqr_from_file(*, file=None, files=None):

    if file and files:
        raise TypeError('Not both!')
    if not files:
        files = [file]
    net_pqr = []

    for i in range(len(files)):

        with open(files[i]) as f:

            pqr = [line.rstrip('\n').split() for line in f]

        f.close()
        pqr_new = np.reshape(pqr, np.shape(pqr))[:, -5:].astype(float)
        net_pqr.append(pqr_new)
    pqr = pqr_new
    print(len(net_pqr))
    print(len(files))

    if len(net_pqr) == 1:
        return pqr
    for j in range(len(net_pqr) - 1):
        pqr = np.vstack((pqr, net_pqr[-j - 2]))
    return pqr

def rmse(calc, ref, decimal_points=3):
    return np.round(np.sqrt(((calc - ref) ** 2).mean()), decimal_points)


def get_mol_data_from_pqr(file='test.pqr'):

    with open(file) as f:

        pqr = [line.rstrip('\n').split() for line in f]

    f.close()
    pqr_new = np.reshape(pqr, np.shape(pqr))
    pqr_str = pqr_new[:, :5]
    pqr_num = pqr_new[:, -5:].astype(float)
    pqr_f = np.hstack((pqr_str, pqr_num))
    return pqr_f


def grow_surf(coeff, offset, suffix='', xyzr_diel='test.xyzr'):

    xyzr_stern = xyzr_diel + '_' + suffix
    X = np.loadtxt(xyzr_diel)
    X[:, 3] *= coeff
    X[:, 3] += offset
    np.savetxt(xyzr_stern, X, fmt='%5.5f')
    return xyzr_stern


def get_SASA_from_msms(molecule, coef, offset=0, probe_radius=1.4, density=2.0, 
                       molXYZR="test.xyzr", 
                       path_to_nlbc_data=nlbc_path, 
                       msmsPath=msms_path):
    
    #dictionary = {'NaCl': 1.415,'NH4Cl': 1.590,'NaBr': 1.490,'Na2SO4': 1.587,\
    #              'LiCl': 1.200,'KCl': 1.595,'K2SO4': 1.827,'CaCl2': 1.493,'NH42SO4': 1.820}
    #salt='NaCl'
    pwd = os.getcwd()
    cp_cmd = 'cp -rf /Users/Ali/repos/setschenow-data/nonpolar/get_surf_info.sh %s' % pwd
    subprocess.Popen(cp_cmd, shell=True).wait()
    chmod_cmd = 'chmod +x ./get_surf_info.sh'
    subprocess.Popen(chmod_cmd, shell=True).wait()
    arr = np.zeros((4,1)).ravel()
    #os.system("cp %s ." % get_surf)
    molXYZR_grown = grow_surf(coef,offset)
    molXYZR_surf = "surf_%s_%s"%(molecule,coef)
    msms_grow = "%s%s%s%s%3.1f%s%3.1f%s%s > msmslog" % (msmsPath,' -if ',molXYZR_grown,' -no_header -de ',density,
                                                 ' -prob ',probe_radius,' -of ',molXYZR_surf)
    #msms_log = "msmslog"
    subprocess.Popen(msms_grow, shell=True).wait()
    #os.system("%s > %s" % (msms_grow,msms_log))
    info_cmd="./get_surf_info.sh msmslog"
    subprocess.Popen(info_cmd, shell=True).wait()
    if os.stat('surf_data').st_size == 0:
        return arr
    arr = np.loadtxt('surf_data', unpack=True).ravel()
    return arr
  

def get_element(string):
    return (re.sub("\d", "", string)).lower()

def bare_area(radius,SASA=True,probe_rad=1.4):
    if SASA:
        radius += probe_rad
    return np.round(4*np.pi*radius**2 ,3)

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + ids + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Failed'

def get_orbital_info(solute,path=nlbc_path):
    file_name = os.path.join(path,solute,'test.prmtop')
    data=[]
    buff=[]
    with open(file_name, "r") as ifile:
        for line in ifile:
            if line.startswith("%FORMAT(20a4)"):
                data.append(next(ifile).split())
                buff.append(next(ifile).split())
    l=[]
    if buff[-2][0]!="%FLAG":
        l.append(data[-2])
        l.append(buff[-2])
        data_out = list(itertools.chain(*l))
    else:
        data_out=data[-2]
        
    if ('no' in data_out and 'o' in data_out):
        indices = [i for i, x in enumerate(data_out) if x == "o"]
        for index in indices:
            data_out[index]="on"
    return list(np.asarray(data_out))


def get_info_msms(molecule,probe_radius = 1.4,density=2.0,molXYZR = "test.xyzr",
                  path_to_nlbc_data = nlbc_path,
                  msmsPath = '/Users/Ali/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1'):
    """Molecular surface data from msms"""    
    pwd = os.getcwd()
    cp_cmd = 'cp -rf /Users/Ali/repos/setschenow-data/nonpolar/get_surf_info.sh %s' % pwd
    subprocess.Popen(cp_cmd, shell=True).wait()
    chmod_cmd = 'chmod +x ./get_surf_info.sh'
    subprocess.Popen(chmod_cmd, shell=True).wait()
    arr = np.zeros((4,1)).ravel()
    molXYZR_surf = "surf_%s"%(molXYZR)
    msms_cmd = "%s%s%s%s%3.1f%s%3.1f%s%s > msmslog" % (msmsPath,' -if ',molXYZR,' -no_header -de ',density,
                                                 ' -prob ',probe_radius,' -of ',molXYZR_surf)
    subprocess.Popen(msms_cmd, shell=True).wait()
    info_cmd="./get_surf_info.sh msmslog"
    subprocess.Popen(info_cmd, shell=True).wait()
    if os.stat('surf_data').st_size == 0:
        return arr
    arr = np.loadtxt('surf_data', unpack=True).ravel()
    return arr

def get_area_centroid_from_points(points):
    vec1 = points[0] - points[1]
    vec2 = points[0] - points[2]
    centroid = np.sum(points, axis=0)/3
    cross = 0.5*np.cross(vec1,vec2)
    area = np.sqrt(cross@cross)
    return area,centroid,cross/np.linalg.norm(cross)
    # return area,centroid
def get_dist(p1,p2):
    dist_vec = p1-p2
    return np.sqrt(dist_vec@dist_vec)