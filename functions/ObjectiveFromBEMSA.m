function [Error,calculatedE,referenceE,electrostatic,nonpolar] = ObjectiveFromBEMSA(x)
%ERROR          Returns the deviance of the experimental MD FEP energy results from the
%               calculated FEP free energy.  OBJECTIVEFUNCTION(Params)
%               takes in i situations and calculates the difference between the MD,
%               experimental energy of solvation, 
%               MD(i) and the calculated energy of solvation E(i).
%               
%               

% "wrap" opt variables into physical simulation variables:
%  convert x vector to our usual NLBC structure. 
Params = MakeParamsStructSA(x);
   
% physical simulations and obtain reference results
[calculatedE,referenceE,electrostatic,nonpolar] = CalculateEnergiesFromBEMSA(Params); 

% "unwrap" physical simulation results into optimization problem
Error = referenceE - calculatedE;

