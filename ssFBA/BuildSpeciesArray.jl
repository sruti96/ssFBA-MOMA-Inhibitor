function BuildSpeciesArray(TXTL, tEND, number_of_species)
#Initialize Species Vector
    Species_array = zeros(tEND+1,number_of_species)
    initial_species = readdlm("config/SpeciesDict/met_mM_ss.dat") #Concentration at steady state from literature
    Species_array[1,:] = initial_species

    #Look up experimental values
    ConcData = TXTL[:ConcDictionary]
    m_idx = ConcData[:Species_index_metabolite]
    aa_idx = ConcData[:Species_index_aa]
    met_conc = ConcData[:Conc_met]
    aa_conc = ConcData[:Conc_aa]

    #Correct Species concentration vector with data
    Species_array[1,m_idx] .= met_conc[1,:]
    Species_array[1,aa_idx] .= aa_conc[1,:]
    Species_array[1,95] = 30.0 #Maltodextrin (20-40mM from ACS Paper)
    Species_array[1,152] = 1e-7 #Initial mRNA level
    Species_array[1,148] = 0.0 #GFP
    Species_array[1,154] = 0.0 #CO2
    Species_array[1,157] = 0.0 #GFP_e
    return Species_array
end
