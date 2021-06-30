function Kinetics(TXTL,Species)

    S = max.(0.0,Species) #Concentration of species (mM)

    #Load enzyme levels (Garenne, 2019)
    E = readdlm("config/KineticDict/enzyme_level_nM.dat")./1e6
    E_idx = findall(E.<2.4e-5) #Find all enzymes that were not reported
    E[E_idx] .= 5e-5 #Set all unknown enzymes to 50nM

    #Convert Kcat vector to (1/hr)
    kcat = readdlm("config/KineticDict/kcat_1over_sec.dat").*3600;
    Km = readdlm("config/KineticDict/Km_mM_matrix.dat");

    #Calculate TXTL rates
    TX = compute_transcription_rate(TXTL)
    TL = compute_translation_rate(S[152],TXTL)
    mRNA_degradation_rate = TXTL[:kdX]
    km_tx = TXTL[:km_tx]

    #Calculate rxn rates
    rxn = zeros(200,1)
    for i = 1:200
        rxn[i] = kcat[i]*E[i]
    end

    #Redefine maltodextrin consumption rate and TXTL rates
    rxn[1] = kcat[1]*E[1] * S[95]/(Km[95,1]+S[95])
    rxn[176] = TX
    rxn[177] = TX * S[31]/(km_tx+S[31]) * S[41]/(km_tx+S[41]) * S[72]/(km_tx+S[72]) * S[142]/(km_tx+S[142])
    rxn[178] = S[152] * mRNA_degradation_rate
    rxn[179] = TL
    rxn[180] = TL

    #Check for negative rates
    rxn = max.(0.0,rxn)
    return rxn
end
