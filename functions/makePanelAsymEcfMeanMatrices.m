function asymBem = makePanelAsymEcfMeanMatrices(surf, bem, pqr)

Kp_NLBC = diag(1./surf.areas) * bem.surfsurfop.K'; %Kp * sigma = mean normal field
			 
dphidnCoul_NLBC = bem.chargesurfop.dphidnCoul;

h_quadrature = diag(surf.areas);

asymBem = struct('Kp_NLBC',Kp_NLBC,...
		 'dphidnCoul_NLBC',dphidnCoul_NLBC,...
		 'h_quadrature',h_quadrature);
