function [phiReac, phiDielBndy, dphiDnDielBndy,x,curA] = ...
    solvePanelConsistentSternAsym(dielSurfData, sternSurfData, pqr, bem, ...
			     epsIn, epsOut, kappa, convFactor, ...
				  asymParams, asymBem)

picardIterLimit = 10;
maxGMRESIter = min(100, size(bem.A,1));
numDielPanels = length(dielSurfData.areas);
numSternPanels = length(sternSurfData.areas);
x = zeros(size(bem.A,1),1); 
rhs = bem.B * pqr.q;

% initial guesses for Cauchy data.  We could use Coulomb field 
phiDielBndy = x(1:numDielPanels);
dphiDnDielBndy = x(numDielPanels+1:2*numDielPanels);
phiSternBndy = x(2*numDielPanels+1:2*numDielPanels+numSternPanels);
dphiDnSternBndy = ones(numSternPanels,1)*(-sum(pqr.q)/ epsOut)/ ...
    sum(sternSurfData.areas);

for picardIter = 1:picardIterLimit
  if picardIter ==1
    [curA, curP] = calcPanelYoonSternAsymPert(asymParams, asymBem, dielSurfData, sternSurfData, bem, ...
					      pqr, epsIn, epsOut, kappa, ...
					      phiDielBndy, dphiDnDielBndy, ...
					      phiSternBndy, dphiDnSternBndy);
    [L,U]=lu(curP);
  else
    curA = calcPanelYoonSternAsymPert(asymParams, asymBem, dielSurfData, sternSurfData, bem, ...
				      pqr, epsIn, epsOut, kappa, ...
				      phiDielBndy, dphiDnDielBndy, ...
				      phiSternBndy, dphiDnSternBndy);
  end
  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, L,U,x);
  fprintf('took %d GMRES iters at picard iteration %d\n', max(iter), ...
	  picardIter)
  %  keyboard
  phiDielBndy = x(1:numDielPanels);
  dphiDnDielBndy = x(numDielPanels+1:2*numDielPanels);
  phiSternBndy  = x(2*numDielPanels+1:2*numDielPanels+numSternPanels);
  dphiDnSternBndy  = x(2*numDielPanels+numSternPanels+1:2*numDielPanels+2*numSternPanels);
end

phiReac = convFactor * bem.C * x;