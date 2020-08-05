

function [phiReac, phiDielBndy, dphiDnDielBndy] = ...
    solveConsistentSternAsym(dielSurfData, sternSurfData, pqr, bem, ...
			     epsIn, epsOut, kappa, convFactor, asymParams)

picardIterLimit = 10;
maxGMRESIter = min(100, size(bem.A,1));
numDielPanels = length(dielSurfData.weights);
numSternPanels = length(sternSurfData.weights);
x = zeros(size(bem.A,1),1); 
rhs = bem.B * pqr.q;

% initial guesses for Cauchy data.  We could use Coulomb field 
phiDielBndy = x(1:numDielPanels);
dphiDnDielBndy = x(numDielPanels+1:2*numDielPanels);
phiSternBndy = x(2*numDielPanels+1:2*numDielPanels+numSternPanels);
dphiDnSternBndy = ones(numSternPanels,1)*(-sum(pqr.q)/ epsOut)/ ...
    sum(sternSurfData.weights);

for picardIter = 1:picardIterLimit
  [curA, curP] = calcYoonSternAsymPert(asymParams, dielSurfData, sternSurfData, bem, ...
				      pqr, epsIn, epsOut, kappa, ...
				       phiDielBndy, dphiDnDielBndy,...
				       phiSternBndy,dphiDnSternBndy);
  [L,U]=lu(curP);
  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, L,U);

  phiDielBndy = x(1:numDielPanels);
  dphiDnDielBndy = x(numDielPanels+1:2*numDielPanels);
  phiSternBndy  = x(2*numDielPanels+1:2*numDielPanels+numSternPanels);
  dphiDnSternBndy  = x(2*numDielPanels+numSternPanels+1:2*numDielPanels+2*numSternPanels);
end

phiReac = convFactor * bem.C * x;