function addProblemSA(name, pqrData, srfFile, chargeDistribution, ...
		      reference, surfArea)
global ProblemSet

% check to make sure that the test charge distributions are the
% right length
if size(chargeDistribution,1) ~= length(pqrData.q)
  fprintf('Error in createProblem: number of charges in chargeDistribution is not equal to the number of atoms in the problem!\n');
  keyboard
end

% then define a variable numTests that says, how many tests (charge
% distributions) do we have.  
numTests = size(chargeDistribution,2);

true = 1;  %  just to keep the boolean nature of "uninitialized"
           %  clear in this function
	   
% when we first say "global ProblemSet" at the uppermost level of
% Matlab, it defines it as a zero-length vector of floats, which
% means I can't say "ProblemSet(1) = struct( ... ).  Try this for
% yourself if it's not clear what the issue is.

newIndex = length(ProblemSet)+1;
newproblem = struct('index',newIndex,...
		    'name',name,...
		    'pqrData',pqrData,...
		    'reference', reference,...
		    'srfFile', srfFile,...
		    'uninitialized',true,...
		    'numTestsInProblem',numTests,...
		    'chargeDistribution',chargeDistribution,...
		    'bemYoonStern',[],...
		    'asymBemPcm',[],...
		    'bemPcm',[],...
		    'srfSternData',[],...
		    'surfArea',surfArea);

if newIndex == 1
  ProblemSet = newproblem;
else
  ProblemSet(newIndex) = newproblem;
end

