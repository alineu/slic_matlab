function problemIndex =  initializeProblem(problem)
global UsefulConstants ProblemSet

% this use of problemIndex is SUBTLE! please ask if the motivation/use
% is not clear, after you think about it for a while.
problemIndex = -1;
for i=1:length(ProblemSet)
  if ProblemSet(i).index == problem.index
    problemIndex = i;
  end
end
if problemIndex==-1
  fprintf('Error in initializeProblem: Could not identify problem in ProblemSet!\n');
  keyboard
end

% Key detail: notice that we use ProblemSet(problemIndex) throughout,
% INSTEAD OF "problem" itself.

if ProblemSet(problemIndex).uninitialized

  % useful constants
  epsIn  = UsefulConstants.epsIn;    % e.g., 1
  epsOut = UsefulConstants.epsOut;   % e.g., 80
  kappa  = UsefulConstants.kappa;    % e.g., 0.0 for non-ionic
  
  % loadSternSrfIntoPanels goes through our custom, BEM for
  % molecules .srf file.  It's worth reading to understand about
  % the file format, and the complications of dealing with meshes.
  ProblemSet(problemIndex).srfSternData = ...
      loadSternSrfIntoPanels(ProblemSet(problemIndex).srfFile);

  % the operators from makePanelBemEcfQualMatrices are used for calculating the
  % electric field just inside the interface
  ProblemSet(problemIndex).bemPcm = ...
      makePanelBemEcfQualMatrices(ProblemSet(problemIndex).srfSternData.dielBndy(1), ...
				  ProblemSet(problemIndex).pqrData, epsIn, epsOut);

  % here is the main integral equation operator, preconditioner, etc.
  ProblemSet(problemIndex).bemYoonStern = makePanelBemSternMatrices(ProblemSet(problemIndex).srfSternData, ...
						  ProblemSet(problemIndex).pqrData, epsIn, ...
						  epsOut, kappa);

  % here are the asymmetry-related perturbations to the problem
  % that are independent of the electric field.  have a look.
  ProblemSet(problemIndex).asymBemPcm = makePanelAsymEcfCollocMatrices(ProblemSet(problemIndex).srfSternData.dielBndy(1), ...
						  ProblemSet(problemIndex).bemPcm, ...
						  ProblemSet(problemIndex).pqrData);

  % unset uninitialized field, so it's understood that the next
  % "test" calculation won't repeat these expensive (and memory
  % hungry!!) steps
  ProblemSet(problemIndex).uninitialized = 0;
end

