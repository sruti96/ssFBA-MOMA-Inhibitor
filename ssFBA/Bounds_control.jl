include("Kinetics.jl")
function Bounds_control(DF,TXTL,Species_vector,idx)

FB = DF["default_flux_bounds_array"]
SBA = DF["species_bounds_array"]
  Species = Species_vector[idx,:]
  Species = max.(0.0,Species) #Make sure no negative species
  #=======================Look up flux and enzyme activity data======================#
  ConcDict = TXTL[:ConcDictionary]
  aa_idx = ConcDict[:Species_index_aa]
  aa_flux = ConcDict[:Flux_aa]
  energy_flux = ConcDict[:Flux_energy]
  gly_flux = ConcDict[:Flux_gly]
  reducing_flux = ConcDict[:Flux_reducing]
  tca_flux = ConcDict[:Flux_tca]
  energy_idx = ConcDict[:index_energy]
  gly_idx = ConcDict[:index_gly]
  reducing_idx = ConcDict[:index_reducing]
  tca_idx = ConcDict[:index_tca]
  rxn_enzyme_idx = ConcDict[:index_enzyme_activity]
  vmax_2h = ConcDict[:enzyme_vmax_2h]
  vmax_8h = ConcDict[:enzyme_vmax_8h]

  #==================================Flux RATES===========================================#
  KineticFlux = Kinetics(TXTL,Species)
  FB[1:200,2] .= DF["Enzyme_sample"] .* KineticFlux[1:200]
  FB[176:180,2] .= KineticFlux[176:180]
  FB[64,2] = 0.0 #DNP LEAK

  #Define kinetic rates for enzyme activity assays
  for i = 1:length(rxn_enzyme_idx)
    if idx <= 60
      FB[rxn_enzyme_idx[i],2] = DF["Enzyme_sample"] * vmax_2h[i]
    elseif idx > 60
      FB[rxn_enzyme_idx[i],2] = DF["Enzyme_sample"] * vmax_8h[i]
    end
  end

  #===============================Update Species Bound Array========================================#
  #Match FB (flux bounds) constraints to SBA (species bounds)
  network = DF["stoichiometric_matrix"];
  SBA[:,1] = max.(-Species,network*FB[:,2])
  SBA[:,2] = network*FB[:,2]

  #Check Bounds
  SBA[:,1] = min.(SBA[:,1],0.0)
  SBA[:,2] = max.(SBA[:,2],0.0)

  #=======================Constrain to experimental data======================================#
  delta = 0.1 #Percent of range to sample on between lower and upper bound
  #AminoAcids
  aa_id_pos = findall(aa_flux[idx,:].>0.0)
  aa_id_neg = findall(aa_flux[idx,:].<0.0)
  aa_tmp = findall(-0.012.<aa_flux[idx,:].<0.0)
  SBA[aa_idx[aa_id_pos],2] .= DF["AA_pos"][aa_id_pos] .* aa_flux[idx,aa_id_pos]
  SBA[aa_idx[aa_id_neg],1] .= max.(-Species[aa_idx[aa_id_neg]],DF["AA_neg"][aa_id_neg] .* aa_flux[idx,aa_id_neg])
  SBA[aa_idx[aa_id_neg],2] .= max.(-Species[aa_idx[aa_id_neg]],(DF["AA_neg"][aa_id_neg] .- delta).*aa_flux[idx,aa_id_neg])
  #Allow for expiremental noise
  SBA[aa_idx[aa_tmp],1] .= max.(-Species[aa_idx[aa_tmp]],DF["AA_pos"][aa_tmp] .* (-0.01))
  SBA[aa_idx[aa_tmp],2] .= 0.0.*DF["AA_pos"][aa_tmp] .* aa_flux[idx,aa_tmp]

  #TCA metabolites
  tca_id_pos = findall(tca_flux[idx,:].>0.0)
  tca_id_neg = findall(tca_flux[idx,:].<0.0)
  tca_tmp = findall(-0.05.<tca_flux[idx,:].<0.05)
  SBA[tca_idx[tca_id_pos],1] .= (0.0) .* tca_flux[idx,tca_id_pos]
  SBA[tca_idx[tca_id_pos],2] .= DF["TCA_pos"][tca_id_pos] .* tca_flux[idx,tca_id_pos]
  SBA[tca_idx[tca_id_neg],1] .= max.(-Species[tca_idx[tca_id_neg]],DF["TCA_neg"][tca_id_neg] .* tca_flux[idx,tca_id_neg])
  SBA[tca_idx[tca_id_neg],2] .= max.(-Species[tca_idx[tca_id_neg]],(DF["TCA_neg"][tca_id_neg] .- delta).*tca_flux[idx,tca_id_neg])

  #glycolysis + PPP metabolites
  gly_id_pos = findall(gly_flux[idx,:].>0.0)
  gly_id_neg = findall(gly_flux[idx,:].<0.0)
  SBA[gly_idx[gly_id_pos],1] .= 0.0*gly_flux[idx,gly_id_pos]
  SBA[gly_idx[gly_id_pos],2] .= DF["GLY_pos"][gly_id_pos] .* gly_flux[idx,gly_id_pos]
  SBA[gly_idx[gly_id_neg],1] .= max.(-Species[gly_idx[gly_id_neg]],DF["GLY_neg"][gly_id_neg] .* gly_flux[idx,gly_id_neg])
  SBA[gly_idx[gly_id_neg],2] .= max.(-Species[gly_idx[gly_id_neg]],(DF["GLY_neg"][gly_id_neg] .- delta).*gly_flux[idx,gly_id_neg])

  #energy metabolites
  energy_id_pos = findall(energy_flux[idx,:].>0.0)
  energy_id_neg = findall(energy_flux[idx,:].<0.0)
  energy_tmp = findall(-0.05.<energy_flux[idx,:].<0.01)
  SBA[energy_idx[energy_id_pos],1] .= 0.0.*energy_flux[idx,energy_id_pos]
  SBA[energy_idx[energy_id_pos],2] .= DF["ENERGY_pos"][energy_id_pos] .* energy_flux[idx,energy_id_pos]
  SBA[energy_idx[energy_id_neg],1] .= max.(-Species[energy_idx[energy_id_neg]],DF["ENERGY_neg"][energy_id_neg] .* energy_flux[idx,energy_id_neg])
  SBA[energy_idx[energy_id_neg],2] .= max.(-Species[energy_idx[energy_id_neg]],0.0.*energy_flux[idx,energy_id_neg])

  #reducing metabolites
  reducing_id_pos = findall(reducing_flux[idx,:].>0.0)
  reducing_id_neg = findall(reducing_flux[idx,:].<0.0)
  SBA[reducing_idx[reducing_id_pos],1] .= 0.0*reducing_flux[idx,reducing_id_pos]
  SBA[reducing_idx[reducing_id_pos],2] .= DF["REDUCING_pos"][reducing_id_pos] .* reducing_flux[idx,reducing_id_pos]
  SBA[reducing_idx[reducing_id_neg],1] .= max.(-Species[reducing_idx[reducing_id_neg]], DF["REDUCING_neg"][reducing_id_neg] .* reducing_flux[idx,reducing_id_neg])
  SBA[reducing_idx[reducing_id_neg],2] .= max.(-Species[reducing_idx[reducing_id_neg]],0.0.*reducing_flux[idx,reducing_id_neg])

  #Check Bounds (Make sure consumption isn't greater than species available)
  SBA[:,1] = max.(SBA[:,1],-Species)

  #Allow nucleotide mono and di phosphates to accumulate (allows accumulation of mRNA)
  SBA[24,2] = .05; #AMP
  SBA[38,2] = .05; #CMP
  SBA[71,2] = .05; #GMP
  SBA[141,2] =.05; #UMP
  SBA[18,2] = .05; #ADP
  SBA[34,2] = .05; #CDP
  SBA[60,2] = .05; #GDP
  SBA[140,2] =.05; #UDP

  SBA[150,2] = 0.0 #RIBOSOME
  SBA[157,1:2] .= 0.0 #gfp_e
  SBA[148,1] = 0.0# gfp
  SBA[148,2] = 100.0# gfp
  FB[213:214,2] .= 0.0 #Ethanol and mglx export
  SBA[53,2] = 0.0 #Formate
  SBA[47,2] = 0.0 #Ethanol

  DF["species_bounds_array"] = SBA
  DF["default_flux_bounds_array"] = FB
  return DF
end
