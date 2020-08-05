function [E,electrostatic,nonpolar,hb,disp,...
          disp_sl_sl,disp_sv_sl,disp_sv_sv,cav,comb] = runTestCosmo(params, problem, chargeDistribution)
global UsefulConstants ProblemSet

epsIn  = UsefulConstants.epsIn;    % solute dielectric constant
epsOut = UsefulConstants.epsOut;   % solvent dielectric constant
kappa  = UsefulConstants.kappa;    % 0.0 for non-ionic
conv_factor = UsefulConstants.conv_factor; % 332.112
kB = 0.001987;
%% make sure all of our BEM stuff is loaded
index = initializeProblem(problem);
%% assign current test's charge distribution to the PQR
pqrData = problem.pqrData;
pqrData.q = chargeDistribution;  
				 
[phiReac, phiBndy, dphiDnBndy] = ...
                  solvePanelConsistentSternAsym(ProblemSet(index).srfSternData.dielBndy(1), ...
				  ProblemSet(index).srfSternData.sternBndy(1), ...
				  pqrData, ProblemSet(index).bemYoonStern, ...
				  epsIn, epsOut, kappa, ...
				  conv_factor, params, ...
				  ProblemSet(index).asymBemPcm);

% phiReac holds the vector of reaction potentials at the charge
% locations (in kcal/mol/e).  This is the response due to
% electrostatic polarization 
dG_asym = 0.5 * pqrData.q' * phiReac;

% the additional term is due to the work done against the static
% potential field that arises due to water structure around even
% uncharged solutes.  we model it as a constant field so the extra
% free energy = potential * totalCharge

% electrostatic
electrostatic = dG_asym + params.phiStatic*sum(pqrData.q);

% now account for the nonpolar solvation and hbonding terms;
% dG_np = dG_disp + dG_cav + dG_comb
solute_data = problem.soluteAtomAreas;
solute_atom_types = problem.soluteAtomTypes;
solute_vdw_v = problem.solute_VdWV;
solute_vdw_a = problem.solute_VdWA;
solute_hbond_data = problem.soluteHbondData;
solvent_vdw_v = problem.solvent_VdWV;
solvent_vdw_a = problem.solvent_VdWA;
solvent_data = problem.solventAtomAreas;
solvent_atom_types = problem.solventAtomTypes;
solvent_hbond_data = problem.solventHbondData;
temp = problem.temperature;
disp_coeffs = params.dispCoeffs;
hbond_coeffs = params.hbondCoeffs;
z_comb = params.zComb;
q_s = params.q_s;
cavity_coeff = params.cavity_coeff;
atom_vols = problem.atom_vols;
area_solute = solute_data(3);
area_solvent = solvent_data(3);
vol_solute = solute_data(4);
vol_solvent = solvent_data(4);

% hbond
hb = CalcHBondE(solute_hbond_data,solvent_hbond_data,hbond_coeffs,temp);

% cavity
cav = kB*temp*CalcCavityE(solute_vdw_a,solute_vdw_v,atom_vols,cavity_coeff,temp);

% dispersion
lambda = 2;
disp_sv_sv = (solute_vdw_v/solvent_vdw_v)*...
              CalcDispersionE(solvent_data,solvent_atom_types,...
                             solvent_data,solvent_atom_types,...
                             disp_coeffs,q_s,lambda,temp);

disp_sv_sl = CalcDispersionE(solute_data,solute_atom_types,...
                             solvent_data,solvent_atom_types,...
                             disp_coeffs,q_s,lambda,temp);
                
disp_sl_sl = (CalcDispersionE(solute_data,solute_atom_types,...
                               solute_data,solute_atom_types,...
                               disp_coeffs,q_s,lambda,temp));
        
disp = 2*disp_sv_sl - disp_sv_sv;  

% combinatorial
comb = -kB*temp*CalcCombinatorialE(solute_vdw_a,solute_vdw_v,...
                    solvent_vdw_a,solvent_vdw_v,z_comb);

% nonpolar
nonpolar = comb + disp + cav;

% total
E = electrostatic + nonpolar + hb;
