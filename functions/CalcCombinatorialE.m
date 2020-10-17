function combinatorialE = CalcCombinatorialE(solute_vdw_a, ...
	solute_vdw_v, solvent_vdw_a, solvent_vdw_v,z)

q = solute_vdw_a/solvent_vdw_a;
r = solute_vdw_v/solvent_vdw_v;
combinatorialE = (1-r+log(r)-(z/2)*(q)*(1-(r/q)+log(r/q)));