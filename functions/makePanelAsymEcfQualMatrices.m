function asymBem = makePanelAsymEcfQualMatrices(surf, bem, pqr, method, ...
						nquad)

if strcmpi(method,'colloc')
  asymBem = makePanelAsymEcfCollocMatrices(surf, bem, pqr);
elseif strcmpi(method,'mean')
  asymBem = makePanelAsymEcfMeanMatrices(surf, bem, pqr);
elseif strcmpi(method,'quad')
  asymBem = makePanelAsymEcfQuadMatrices(surf, bem, pqr, nquad);
else 
  fprintf('Error: makePanelAsymEcfQualMatrices\n');
  fprintf('arg "method" must be one of ["colloc", "mean", "quad"]\n');
  asymBem = 0;
  return
end

