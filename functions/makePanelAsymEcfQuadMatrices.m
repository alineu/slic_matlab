function asymBem = makePanelAsymEcfQuadMatrices(surf, bem, pqr, nquad)

[points, normals, weights] = makePanelQuadraturePoints(surf, ...
						  nquad);

[V, K, Kp] = colloc_Laplace(surf.meshData, points, normals, weights);

Kp_NLBC = Kp; % i.e. Kp * sigma will be the electric field at the
              % quadrature points
			 
[phiCoul_NLBC, dphidnCoul_NLBC] = coulombMatrixCollocation(surf.meshData, ...
						  points, normals, pqr.xyz);


nq = nquad;
np = length(surf.areas);
iVec = reshape(kron(ones(1,nq),(1:np)')',[1 nq*np]);
jVec = repmat(1:nq,[1,np]) + (iVec-1)*nq;
h_quadrature = sparse(iVec,jVec,weights,np,np*nq,length(weights));

asymBem = struct('Kp_NLBC',Kp_NLBC,...
		 'dphidnCoul_NLBC',dphidnCoul_NLBC,...
		 'h_quadrature',h_quadrature);
