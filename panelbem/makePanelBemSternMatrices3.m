function bem = makePanelBemSternMatrices3(surf1, surf2, surf3, pqr1, pqr2, epsIn1, epsIn2, epsOut, ...
					 kappaArg)
global kappa;
kappa = kappaArg;

npDiel1 = length(surf1.dielBndy(1).areas);
npDiel2 = length(surf2.dielBndy(1).areas);

%Surface3 is the merged surface of the two molecules (atoms) and only 
% the stern surface will be used
npStern3 = length(surf3.dielBndy(1).areas);

%npTotal = npDiel1 + npStern1 + npDiel2 + npStern2;

nCharges1 = length(pqr1.q);
nCharges2 = length(pqr2.q);

diel1ChargeOp = makeSurfaceToChargePanelOperators(surf1.dielBndy(1), ...
						 pqr1);
diel2ChargeOp = makeSurfaceToChargePanelOperators(surf2.dielBndy(1), ...
						 pqr2);
Idiel1 = eye(npDiel1);
Idiel2 = eye(npDiel2);
Istern3 = eye(npStern3);



% src = diel1, target = diel1 
diel1Diel1Op = makeSurfaceToSurfaceLaplacePanelOperators(surf1.dielBndy(1));

% src = diel1, target = diel2 
diel1Diel2Op = makeSurfaceToSurfaceLaplacePanelOperators(surf1.dielBndy(1), ...
						  surf2.dielBndy(1));

% src = diel1, target = stern3
diel1ToStern3Op = makeSurfaceToSurfaceLaplacePanelOperators(surf1.dielBndy(1), ...
						  surf3.sternBndy(1));

% src = stern3, target = diel1
stern3ToDiel1Op = makeSurfaceToSurfaceLaplacePanelOperators(surf3.sternBndy(1), ...
					      surf1.dielBndy(1));
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
% src = diel2, target = diel2 
diel2Diel2Op = makeSurfaceToSurfaceLaplacePanelOperators(surf2.dielBndy(1));

% src = diel2, target = diel1 
diel2Diel1Op = makeSurfaceToSurfaceLaplacePanelOperators(surf2.dielBndy(1), ...
						  surf1.dielBndy(1));

% src = diel2, target = stern3
diel2ToStern3Op = makeSurfaceToSurfaceLaplacePanelOperators(surf2.dielBndy(1), ...
						  surf3.sternBndy(1));

% src = stern3, target = diel2
stern3ToDiel2Op = makeSurfaceToSurfaceLaplacePanelOperators(surf3.sternBndy(1), ...
					      surf2.dielBndy(1));
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAPLACE kernel for src = stern3, target = stern3                     

stern3Stern3Laplace = makeSurfaceToSurfaceLaplacePanelOperators(surf3.sternBndy(1));

% YUKAWA kernel for src = stern3, target = stern3

if kappa > 1e-6
  stern3Stern3Yukawa = ...
      makeSurfaceToSurfaceYukawaPanelOperators(surf3.sternBndy(1), kappa);
else
  stern3Stern3Yukawa = ...
      makeSurfaceToSurfaceLaplacePanelOperators(surf3.sternBndy(1));
end

ZD1D2 = zeros(npDiel1, npDiel2);
ZD1S3 = zeros(npDiel1, npStern3);
ZD2D1 = zeros(npDiel2, npDiel1);
ZD2S3 = zeros(npDiel2, npStern3);
ZS3D1 = zeros(npStern3, npDiel1);
ZS3D2 = zeros(npStern3, npDiel2);
% could just use Z' but this is clearer

%A
A11 = (Idiel1/2 + diel1Diel1Op.K);
A12 = -diel1Diel1Op.V;
A13 = ZD1D2;
A14 = ZD1D2;
A15 = ZD1S3;
A16 = ZD1S3;

A21 = (Idiel1/2 - diel1Diel1Op.K);
A22_base = diel1Diel1Op.V;
A23 = -diel2Diel1Op.K;
A24_base = diel2Diel1Op.V;
A25 = stern3ToDiel1Op.K;
A26 = -stern3ToDiel1Op.V;

A31 = ZD2D1;
A32 = ZD2D1;
A33 = (Idiel2/2 + diel2Diel2Op.K);
A34 = -diel2Diel2Op.V;
A35 = ZD2S3;
A36 = ZD2S3;

A41 = -diel1Diel2Op.K;
A42_base = diel1Diel2Op.V; 
A43 = (Idiel2/2 - diel2Diel2Op.K);
A44_base = diel2Diel2Op.V;
A45 = stern3ToDiel2Op.K;
A46 = -stern3ToDiel2Op.V;

A51 = -diel1ToStern3Op.K;
A52_base = diel1ToStern3Op.V; 
A53 = -diel2ToStern3Op.K;
A54_base = diel2ToStern3Op.V; 
A55 = (Istern3/2+stern3Stern3Laplace.K);
A56 = -stern3Stern3Laplace.V;

A61 = ZS3D1;
A62 = ZS3D1;
A63 = ZS3D2;
A64 = ZS3D2;
A65 = (Istern3/2-stern3Stern3Yukawa.K);
A66_base = stern3Stern3Yukawa.V;


A = [A11 A12 A13 A14 A15 A16;
     A21 (epsIn1/epsOut)*A22_base A23 (epsIn2/epsOut)*A24_base A25 A26;
     A31 A32 A33 A34 A35 A36;
     A41 (epsIn1/epsOut)*A42_base A43 (epsIn2/epsOut)*A44_base A45 A46;
     A51 (epsIn1/epsOut)*A52_base A53 (epsIn2/epsOut)*A54_base A55 A56;
     A61 A62 A63 A64 A65 (epsOut/epsOut)*A66_base;];
     
B1 = [diel1ChargeOp.phiCoul/epsIn1;
     zeros(npDiel1,nCharges1);
     zeros(npDiel2,nCharges1);
     zeros(npDiel2,nCharges1);
     zeros(npStern3,nCharges1);
     zeros(npStern3,nCharges1);];
     
B2 = [zeros(npDiel1,nCharges2);
     zeros(npDiel1,nCharges2);
     diel2ChargeOp.phiCoul/epsIn2;
     zeros(npDiel2,nCharges2);
     zeros(npStern3,nCharges2);
     zeros(npStern3,nCharges2);];


C1 = 4*pi*[-diel1ChargeOp.dlpToCharges diel1ChargeOp.slpToCharges ...
     zeros(nCharges1,npDiel2) zeros(nCharges1, npDiel2) ...
     zeros(nCharges1,npStern3) zeros(nCharges1, npStern3);];

C2 = 4*pi*[zeros(nCharges2,npDiel1) zeros(nCharges2, npDiel1) ...
     -diel2ChargeOp.dlpToCharges diel2ChargeOp.slpToCharges ...
     zeros(nCharges2,npStern3) zeros(nCharges2, npStern3);];
 
bem = struct('A', A,'B1', B1, 'B2', B2 ,...
         'C1',C1 ,'C2', C2 , ...
	     'diel1ChargeOp',diel1ChargeOp,...
	     'diel2ChargeOp',diel2ChargeOp,...
	     'diel1Diel1Op',diel1Diel1Op,...
	     'diel1Diel2Op',diel1Diel2Op,...
	     'diel2Diel2Op',diel2Diel2Op,...
	     'diel2Diel1Op',diel2Diel1Op,...
	     'diel1ToStern3Op',diel1ToStern3Op,...
	     'diel2ToStern3Op',diel2ToStern3Op,...
	     'stern3ToDiel1Op',stern3ToDiel1Op,...
	     'stern3ToDiel2Op',stern3ToDiel2Op,...
	     'stern3Stern3Laplace',stern3Stern3Laplace,...
	     'stern3Stern3Yukawa',stern3Stern3Yukawa,...
	     'A11',A11,'A12',A12,'A13',A13,'A14',A14,...
	     'A15',A15,'A16',A16,...
	     'A21',A21,'A22_base',A22_base,'A23',A23,'A24_base',A24_base,...
	     'A25',A25,'A26',A26,...
         'A31',A31,'A32',A32,'A33',A33,'A34',A34,...
	     'A35',A35,'A36',A36,...
	     'A41',A41,'A42_base',A42_base,'A43',A43,'A44_base',A44_base,...
	     'A45',A45,'A46',A46,...
	     'A51',A51,'A52_base',A52_base,'A53',A53,'A54_base',A54_base,...
	     'A55',A55,'A56',A56,...
         'A61',A61,'A62',A62,'A63',A63,'A64',A64,...
	     'A65',A65,'A66_base',A66_base);
