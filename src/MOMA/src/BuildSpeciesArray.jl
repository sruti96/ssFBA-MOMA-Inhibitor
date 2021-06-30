function BuildSpeciesArray(TXTL, tEND, number_of_species)

    # initialize Species Vector
    species_array = zeros(tEND+1,number_of_species)
    initial_species = readdlm("config/SpeciesDict/met_mM_ss.dat") # concentration at steady state from literature
    species_array[1,:] = initial_species

    # look up experimental values
    ConcData = TXTL["ConcDictionary"]
    m_idx = ConcData["Species_index_metabolite"]
    aa_idx = ConcData["Species_index_aa"]
    met_conc = ConcData["Conc_met"]
    aa_conc = ConcData["Conc_aa"]

    # correct Species concentration vector with data
    species_array[1,m_idx] .= met_conc[1,:]
    species_array[1,aa_idx] .= aa_conc[1,:]
    species_array[1,95] = 30.0 # maltodextrin (20-40mM from ACS Paper)
    species_array[1,152] = 1e-7 # initial mRNA level
    species_array[1,148] = 0.0 # GFP
    species_array[1,154] = 0.0 # CO2
    species_array[1,157] = 0.0 # GFP_e

    return species_array
end
