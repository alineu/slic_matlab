function dispE = CalcDispersionE(solute_data, solute_atom_types, ...
                                 solvent_data, solvent_atom_types, ...
                                 disp_coeffs, q_s, lambda, temperature)
                                      
density_ratio = density_water(temperature)/density_water(298);
all_atom_types = {'br','c-sp','c-sp2','c-sp3','cl',...
                  'f','h','i','n-sp','n-sp3',...
                  'n-sp2','o-sp2','o-sp3-h',...
                  'o-sp2-n','o-sp3','p','s'};
              
all_atom_radii = {2.16, 2,2, 2,2.05,...
                  1.72, 1.3, 2.32, 1.83, 1.83,...
                  1.83, 1.72, 1.72,...
                  1.72, 1.72, 2.12, 2.16};
              
kB = 0.001987;
beta=1/kB/temperature; % Boltzmann factor
solute_areas = solute_data(7:end);
unique_atom_types_solute = unique(string(solute_atom_types));
unique_atom_types_solute = unique_atom_types_solute(2:end);
vol_mix = solvent_data(4);
solvent_areas = solvent_data(7:end);
unique_atom_types_solvent = unique(string(solvent_atom_types));
unique_atom_types_solvent = unique_atom_types_solvent(2:end);
solute_types = cell2table(solute_atom_types);
solvent_types = cell2table(solvent_atom_types);

dispE=0;

for i=1:length(unique_atom_types_solvent)
    
    solvent_atom_type = unique_atom_types_solvent{i};
    [~,id_coeff_solvent] = ismember(solvent_atom_type,all_atom_types);
    coeff_solvent = disp_coeffs(id_coeff_solvent);
    radii_solvent = all_atom_radii(id_coeff_solvent);
    radii_solvent_value = cell2table(radii_solvent);
    radii_solvent = radii_solvent_value.Variables;
    [m, ~] = ismember(solvent_types.Variables,{solvent_atom_type});
    areas_solvent = solvent_areas(m);
    m_tau_solvent = sum(areas_solvent.^q_s);
    dispSoluteSolvent=0;
    
    for j=1:length(unique_atom_types_solute)
        
        solute_atom_type = unique_atom_types_solute{j};
        [~,id_coeff_solute] = ismember(solute_atom_type,all_atom_types);
        coeff_solute = disp_coeffs(id_coeff_solute);
        radii_solute = all_atom_radii(id_coeff_solute);
        radii_solute_value = cell2table(radii_solute);
        radii_solute = radii_solute_value.Variables;
        radii_ij = (radii_solvent+radii_solute)/2;
        disp_coeff_ij = sqrt(coeff_solute*coeff_solvent);
        square_well_exp = -radii_ij^3*(lambda^3-1)*...
                           disp_coeff_ij*(1+disp_coeff_ij*beta/2);
        [m, ~] = ismember(solute_types.Variables,{solute_atom_type});
        areas_solute = solute_areas(m);
        m_tau_solute = sum(areas_solute.^q_s);
        dispSoluteSolvent = dispSoluteSolvent + m_tau_solute*square_well_exp;
        
    end
    
    dispE = dispE + m_tau_solvent*dispSoluteSolvent;
    
end

dispE = density_ratio*dispE/vol_mix; %in kcal/mol
