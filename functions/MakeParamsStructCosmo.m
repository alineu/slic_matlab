function Params = MakeParamsStructCosmo(x)
Params = struct('alpha', x(1),'beta', x(2),'EfieldOffset', x(3), ...
		'mu', x(4), 'phiStatic', x(5), ...
        'dispCoeffs', x(6:22),'q_s', x(23), ...
        'hbondCoeffs', x(24:33), ...
        'zComb', x(34),'cavity_coeff', x(35));
        
fprintf('Params: \n\n');
fprintf('alpha = %0.3f \n', x(1));
fprintf('beta = %0.3f \n', x(2));
fprintf('gamma = %0.3f \n', x(3));
fprintf('mu = %0.3f \n', x(4));
fprintf('phi_static = %0.3f \n', x(5));
fprintf('Br dispersion coeff = %0.3f \n', x(6));
fprintf('C-sp dispersion coeff = %0.3f \n', x(7));
fprintf('C-sp2 dispersion coeff = %0.3f \n', x(8));
fprintf('C-sp3 dispersion coeff = %0.3f \n', x(9));
fprintf('Cl dispersion coeff = %0.3f \n', x(10));
fprintf('F dispersion coeff = %0.3f \n', x(11));
fprintf('H dispersion coeff = %0.3f \n', x(12));
fprintf('I dispersion coeff = %0.3f \n', x(13));
fprintf('N-sp dispersion coeff = %0.3f \n', x(14));
fprintf('N-sp2 dispersion coeff = %0.3f \n', x(15));
fprintf('N-sp3 dispersion coeff = %0.3f \n', x(16));
fprintf('O-sp2 dispersion coeff = %0.3f \n', x(17));
fprintf('O-sp2-N dispersion coeff = %0.3f \n', x(18));
fprintf('O-sp3 dispersion coeff = %0.3f \n', x(19));
fprintf('O-sp3-H dispersion coeff = %0.3f \n', x(20));
fprintf('P dispersion coeff = %0.3f \n', x(21));
fprintf('S dispersion coeff = %0.3f \n', x(22));
fprintf('q_s dispersion coeff = %0.3f \n', x(23));
fprintf('n_amn_hoh hb coeff = %0.3f \n', x(24));
fprintf('n_amd_hoh hb coeff = %0.3f \n', x(25));
fprintf('n_no2_hoh hb coeff = %0.3f \n', x(26));
fprintf('n_other_hoh hb coeff = %0.3f \n', x(27));
fprintf('o_crbnl_hoh hb coeff = %0.3f \n', x(28));
fprintf('o_estr_hoh hb coeff = %0.3f \n', x(29));
fprintf('o_no2_hoh hb coeff = %0.3f \n', x(30));
fprintf('o_oh_hoh hb coeff = %0.3f \n', x(31));
fprintf('fl_hoh hb coeff = %0.3f \n', x(32));
fprintf('o_oh_hnh hb coeff = %0.3f \n', x(33));
fprintf('z combinatorial coeff = %0.3f \n', x(34));
fprintf('cavity rescaling coeff = %0.3f \n\n', x(35));


