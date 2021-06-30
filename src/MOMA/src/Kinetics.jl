# ----------------------------------------------------------------------------------- #
# Copyright (c) 2021 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

function compute_transcription_rate(TXTL_dictionary::Dict{String,Any})

  # compute rX_hat bound -
  VX = TXTL_dictionary["VX"]
  length_factor_transcription = TXTL_dictionary["length_factor_transcription"]
  gene_concentration = TXTL_dictionary["gene_concentration"]
  X_tau_factor = TXTL_dictionary["X_tau_factor"]
  KX = TXTL_dictionary["KX"]

  # compute the control variable value -
  u_variable = calculate_transcription_control(TXTL_dictionary)

  # kinetic limit -
  transcription_rate = VX*length_factor_transcription*(gene_concentration/(KX*X_tau_factor+(1+X_tau_factor)*gene_concentration))/1000 #mM/hr

  # return transcription rate -
  return transcription_rate*u_variable
end

function compute_translation_rate(mRNA, TX, TXTL_dictionary::Dict{String,Any})

  # Get some stuff from the parameter dictionary -
  VL = TXTL_dictionary["VL"]
  length_factor_transcription = TXTL_dictionary["length_factor_translation"]
  L_tau_factor = TXTL_dictionary["L_tau_factor"]
  KL = TXTL_dictionary["KL"]

  # compute the control variable value -
  w_variable = calculate_translation_control(TXTL_dictionary)

  # kinetic limit -
  translation_rate = VL*length_factor_transcription*(mRNA/(KL*L_tau_factor+(1+L_tau_factor)*mRNA))/1000 #mM/hr

  # return translation_rate-
  return translation_rate*w_variable
end

function calculate_kinetic_limits(TXTL,Species)

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
    TL = compute_translation_rate(S[152],TX, TXTL)
    mRNA_degradation_rate = TXTL["kdX"]
    km_tx = TXTL["km_tx"]

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
