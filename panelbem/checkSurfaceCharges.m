addpath('../pointbem');
addpath('../testasymmetry');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
kappa = 0.125;

res = 'phe';
pdbFile = 'residues/phe/phe.pdb';
crgFile = 'residues/phe/phe.crg';
pqrData = loadPdbAndCrg(pdbFile,crgFile);

srfNoSternFile = 'residues/phe/phe_scaledcharmm_2.srf'
srfData = loadSrfIntoPanels(srfNoSternFile);
numPanels =  length(surfData.areas);

srfSternFile   = 'residues/phe/phe_scaledcharmm_stern_2.srf'
srfSternData = loadSternSrfIntoPanels(srfSternFile);

% part 1: no-salt : PCM formulation
bemPcm = makePanelBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);

% part 2: no-salt : YL formulation
bemYoonDiel = makePanelBemYoonDielMatrices(srfData,pqrData,epsIn,epsOut);

surfaceData = bemYoonDiel.A \ (bemYoonDiel.B * pqrData.q);
potentials = surfaceData(1:numPanels);
surfaceChargeDensity = bemYoonDiel.dielDielOp.V \ potentials;
pointCharges = surfaceChargeDensity .* srfData.areas';

for i=1:length(pointCharges)
  for j=1:length(pqrData.q)
    coulMatrix(j,i) = 1/norm(srfData.centroids(i,:) - pqrData.xyz(j,:));
  end
end


%% part 3: salt with Stern, YL/NLBC formulation
%bemYoonStern = makePanelBemSternMatrices(srfSternData, ...
%					 pqrData, epsIn, ...
%					 epsOut, kappa);



