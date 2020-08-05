function [E,electrostatic,nonpolar] = runTestSA(params, problem, chargeDistribution)
global UsefulConstants ProblemSet

%alpha = params.alpha;
%beta  = params.beta;
%gamma = params.EfieldOffset;
%mu    =  -alpha * tanh(-gamma);

epsIn  = UsefulConstants.epsIn;    % 1
epsOut = UsefulConstants.epsOut;   % 80
kappa  = UsefulConstants.kappa;    % 0.0 for non-ionic
conv_factor = UsefulConstants.conv_factor; % 332.112
staticpotential = UsefulConstants.staticpotential;  % 10.7 kcal/mol/e,
                                                    % according to our
                                                    % earlier work

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
electrostatic = dG_asym + params.phiStatic*sum(pqrData.q);

% now account for the nonpolar solvation term
nonpolar = params.surfAreaConstant * problem.surfArea + params.NPoffset;
E = electrostatic + nonpolar;
