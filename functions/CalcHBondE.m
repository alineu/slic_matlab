function hbondE = CalcHBondE(hbond_data, hbond_coeffs, temp)
kB = 0.001987;
hb_ids = find(hbond_data);
num_hbs = hbond_data(hb_ids);
hb_coeffs = hbond_coeffs(hb_ids);
hbondE = sum(num_hbs.*hb_coeffs.*exp(-hb_coeffs/(kB*temp))./(exp(-hb_coeffs/(kB*temp))+1));
