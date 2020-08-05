function [Error,calculatedE,referenceE,electrostatic,nonpolar,...
          hb,disp,disp_sl_sl,disp_sv_sl,disp_sv_sv,cav,comb] = ...
          ObjectiveFromBEMCosmo(x)
%ERROR          Returns the deviance of the experimental MD FEP energy results from the
%               calculated FEP free energy.  OBJECTIVEFUNCTION(Params)
%               takes in i situations and calculates the difference between the MD,
%               experimental energy of solvation, 
%               MD(i) and the calculated energy of solvation E(i).
%               
%               

% "wrap" opt variables into physical simulation variables:
%  convert x vector to our usual NLBC structure. 
Params = MakeParamsStructCosmo(x);
   
% physical simulations and obtain reference results
[calculatedE,referenceE,electrostatic,nonpolar,hb,disp,...
    disp_sl_sl,disp_sv_sl,disp_sv_sv,cav,comb] = CalculateEnergiesFromBEMCosmo(Params); 

% "unwrap" physical simulation results into optimization problem
Error = referenceE - calculatedE;