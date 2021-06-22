function hbondE = CalcHBondE(hbond_data, hbond_coeffs, temp)
kB = 0.001987;
if isempty(hbond_data{2})
    
    hb_ids = find(hbond_data{1});
    num_hbs = hbond_data{1}(hb_ids);
    hb_coeffs = hbond_coeffs(hb_ids);
    hbondE = sum(num_hbs.*hb_coeffs.*exp(-hb_coeffs/(kB*temp))./(exp(-hb_coeffs/(kB*temp))+1));
else
    hb_ids_solute = find(hbond_data{1});
    hb_ids_solvent = find(hbond_data{2});
    vol_ratio = hbond_data{3};
    num_hbs_solute = hbond_data{1}(hb_ids_solute);
    num_hbs_solvent = hbond_data{2}(hb_ids_solvent);
    hb_coeffs_solute = hbond_coeffs(hb_ids_solute);
    hb_coeffs_solvent = hbond_coeffs(hb_ids_solvent);
    hbondE = sum(num_hbs_solute.*hb_coeffs_solute.*exp(-hb_coeffs_solute/(kB*temp))./(exp(-hb_coeffs_solute/(kB*temp))+1)) - ...
             vol_ratio*sum(num_hbs_solvent.*hb_coeffs_solvent.*exp(-hb_coeffs_solvent/(kB*temp))./(exp(-hb_coeffs_solvent/(kB*temp))+1));
             
end