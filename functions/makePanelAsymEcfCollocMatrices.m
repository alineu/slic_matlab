function asymBem = makePanelAsymEcfCollocMatrices(surf, bem, pqr)

Kp_NLBC = bem.surfsurfop.Kp; % i.e. Kp * sigma will be the electric
                         % field at the centroid.
			 
dphidnCoul_NLBC = bem.chargesurfop.dphidnCoul;

h_quadrature = diag(surf.areas);

asymBem = struct('Kp_NLBC',Kp_NLBC,...
		 'dphidnCoul_NLBC',dphidnCoul_NLBC,...
		 'h_quadrature',h_quadrature);
