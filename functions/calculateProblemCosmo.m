function [calculatedE, referenceE, electrostatic, nonpolar, ...
          dG_hb, dG_disp, dG_disp_sl_sl, dG_disp_sv_sl, ...
          dG_disp_sv_sv, dG_cav, dG_comb] = ...
          calculateProblemCosmo(problem, params)

% This function is intended to be called by
% CalculateEnergiesFromBEM.  you are supposed to pass in a problem
% (one element from ProblemSet) and a param structure.
numTestsInProblem = problem.numTestsInProblem;
calculatedE = zeros(numTestsInProblem, 1);

% Run a BEM calculation for each test charge distribution.  Note,
% importantly, i have a subtle but nice innovation here.  runTest
% will initialize the BEM matrices if it needs to, but not
% otherwise! this will save a ton of time
%keyboard
for i=1:numTestsInProblem
  [calculatedE(i),electrostatic(i),nonpolar(i),dG_hb(i),dG_disp(i),...
   dG_disp_sl_sl(i),dG_disp_sv_sl(i),dG_disp_sv_sv(i),dG_cav(i),dG_comb(i)] = ...
   runTestCosmo_2(params, problem, problem.chargeDistribution(:,i));
end

% Here's where the information about the reference result (whether
% experiment, MD, or other, e.g. MSA) comes in.  the info is
% encapsulated "inside" the problem structure, which improves modularity.
referenceE = problem.reference;

% Double check the reference results and calculated results are the
% same length.  this should be impossible to fail because of the
% check in addProblem.
if length(calculatedE) ~= length(referenceE)
  fprintf('Error in calculate problem %s!\n',problem.name);
  fprintf('length(calculatedE)=%d\nlength(referenceE)=%d\n',length(calculatedE),length(referenceE));
  keyboard
end
