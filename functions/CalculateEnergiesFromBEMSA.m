function [calculatedE, referenceE,electrostatic,nonpolar ] = CalculateEnergiesFromBEMSA(Params)
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName

% define empty vectors so we can extend them easily
calculatedE  = [];
referenceE = [];
electrostatic  =[];
nonpolar = [];
% how many problem geometries do we have ?  note that this can be
% DIFFERENT FROM the number of overall "BEM calculations" we do.

numProblems = length(ProblemSet);

for i=1:numProblems
    i
  curProblem = ProblemSet(i);
  [newCalculatedE, newReferenceE,newElectrostatic,newNonpolar]  = calculateProblemSA(curProblem, Params);

if saveMemory
  ProblemSet(i).bemPcm = [];
  ProblemSet(i).bemYoonStern = [];
  ProblemSet(i).asymBemPcm = [];
  ProblemSet(i).uninitialized = 1;
end

if writeLogfile
  fid = fopen(logfileName,'a');
  fprintf(fid,'%s,%f,%f,%f,%f\n',curProblem.name,newReferenceE(1), ...
	  newCalculatedE(1),newElectrostatic(1),newNonpolar(1));
  fclose(fid);
end

  % append calculated and reference results to total list  
  calculatedE = [calculatedE; newCalculatedE];
  referenceE = [referenceE; newReferenceE];
  electrostatic= [electrostatic; newElectrostatic];
  nonpolar = [nonpolar; newNonpolar];
end

% there's no need for a trailing "end" at the end of a Matlab
% function, when it's in a .m file by itself!