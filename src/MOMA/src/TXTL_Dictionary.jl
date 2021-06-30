function generate_TXTL_dictionary(TXTL_parameters)

  # Protein of Interest Sequences
  gene_concentration = 0.005 #μM from Experiment
  volume_in_L = 14e-6

  length_of_gene_in_nt = 717 #nt
  characteristic_gene_length = 1000 #nt
  length_factor_transcription = (characteristic_gene_length/length_of_gene_in_nt)

  length_of_prot_in_aa = 239 #aa
  characteristic_prot_length = 330 #aa
  length_factor_translation = (characteristic_prot_length/length_of_prot_in_aa)

  # Setup the mRNA elongation rate, and global translation
  RNAP_concentration = TXTL_parameters[1]; #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
  max_transcription_rate = TXTL_parameters[2];  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
  RIBOSOME_concentration = TXTL_parameters[3]; #RIBOSOME_concentration #uM  <0.0023mM (ACS SynBio Garamella 2016)
  max_translation_rate = TXTL_parameters[4]; #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

  # calculate the time time constant for transcription -
  transcription_initiation_time_constant = 1 #s" Varner has 45s
  kcat_transcription = max_transcription_rate*(3600/length_of_gene_in_nt)                         # hr^-1
  kcat_transcription_initiation = (1/transcription_initiation_time_constant)*(3600)              # hr^-1
  X_tau_factor = (kcat_transcription)/(kcat_transcription_initiation)

  # calculate the time constant for translation -
  translation_initiation_time_constant=15 #s #bionumbers_reference":"BIND:109525 (60sec Underwood, Swartz, Puglisi 2005 Biotech Bioeng)
  kcat_translation = max_translation_rate*(3600/length_of_prot_in_aa)                             # hr^-1
  kcat_translation_initiation = (1/translation_initiation_time_constant)*(3600)                  # hr^-1
  L_tau_factor = (kcat_translation)/(kcat_translation_initiation)

  # calculate the degradation constants for mRNA and protein -
  mRNA_half_life = 0.290833 #h (0.290833 hr ACS SynBio Garamella 2016)
  kdX = -(1/mRNA_half_life)*log(0.5)

  protein_half_life =70 #h
  kdL = -(1/protein_half_life)*log(0.5)

  # compute VX and VL -
  polysome_gain = 10;
  kEX = max_transcription_rate*(3600/characteristic_gene_length)  # hr^-1
  kEL = max_translation_rate*(3600/characteristic_prot_length)    # hr^-1
  VX = kEX*RNAP_concentration                                     # μM/hr
  VL = kEL*RIBOSOME_concentration*polysome_gain                   # μM/hr

  # compute saturation coefficient for X and L -
  elongation_slope = 0.3 # μM "McClure 1980= 0.3"
  KX = 0.3 # elongation_slope  # μM
  KL = 0.6 # mM (600uM)
  km_tx = 0.03 #Saturation constant for nucleotides involved in TX

  # setup promoter model -
  Promoter = Dict{String,Any}()
  Promoter["K1"] = 0.014
  Promoter["K2"] = 10
  Promoter["coop_parameter"] = 1
  Promoter["k_binding"] = 130
  Promoter["inducer_concentration"] = 35.0

  # =============================== DO NOT EDIT BELOW THIS LINE ============================== #
  # stuff for bounds -
  TXTL_dictionary  = Dict{String,Any}()
  TXTL_dictionary["gene_concentration"] = gene_concentration
  TXTL_dictionary["length_factor_transcription"] = length_factor_transcription
  TXTL_dictionary["length_factor_translation"] = length_factor_translation
  TXTL_dictionary["RX_concentration"] = RNAP_concentration
  TXTL_dictionary["RL_concentration"] = RIBOSOME_concentration
  TXTL_dictionary["X_tau_factor"] = X_tau_factor
  TXTL_dictionary["L_tau_factor"] = L_tau_factor
  TXTL_dictionary["kdX"] = kdX
  TXTL_dictionary["kdL"] = kdL
  TXTL_dictionary["KX"] = KX
  TXTL_dictionary["KL"] = KL
  TXTL_dictionary["kdL"] = kdL
  TXTL_dictionary["VX"] = VX
  TXTL_dictionary["VL"] = VL
  TXTL_dictionary["km_tx"] = km_tx
  TXTL_dictionary["Promoter"] = Promoter
  # ============================== DO NOT EDIT ABOVE THIS LINE ============================== #

  return TXTL_dictionary
end
