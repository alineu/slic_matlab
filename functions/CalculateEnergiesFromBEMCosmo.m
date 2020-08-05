function [calculatedE, referenceE,electrostatic,nonpolar,...
          hb,disp,disp_sl_sl,disp_sv_sl,disp_sv_sv,cav,comb] = CalculateEnergiesFromBEMCosmo(Params)
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName

% define empty vectors so we can extend them easily
calculatedE  = [];
referenceE = [];
electrostatic  =[];
nonpolar = [];
hb = [];
disp = [];
disp_sl_sl= [];
disp_sv_sl= [];
disp_sv_sv= [];
cav = [];
comb = [];

% how many problem geometries do we have ?  note that this can be
% DIFFERENT FROM the number of overall "BEM calculations" we do.

numProblems = length(ProblemSet);

for i=1:numProblems
    i
  curProblem = ProblemSet(i);
  [newCalculatedE, newReferenceE,newElectrostatic,newNonpolar,newHb,newDisp,...
      newDisp_sl_sl,newDisp_sv_sl,newDisp_sv_sv,newCav,newComb] = ...
      calculateProblemCosmo(curProblem, Params);

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
  hb = [hb; newHb];
  disp = [disp; newDisp];
  disp_sl_sl = [disp_sl_sl; newDisp_sl_sl];
  disp_sv_sl = [disp_sv_sl; newDisp_sv_sl];
  disp_sv_sv = [disp_sv_sv; newDisp_sv_sv];
  cav = [cav; newCav];
  comb = [comb; newComb];
  
end

% there's no need for a trailing "end" at the end of a Matlab
% function, when it's in a .m file by itself!