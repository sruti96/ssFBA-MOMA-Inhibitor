# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
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
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2019-06-08T16:39:01.235
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("config/Network.dat");

	# Setup default flux bounds array -
	default_bounds_array = [
		0	100.0	;	# 1 M_maltodextrin6_c+M_pi_c --> M_maltodextrin5_c+M_g1p_c
		0	100.0	;	# 2 M_maltodextrin5_c+M_pi_c --> M_maltodextrin4_c+M_g1p_c
		0	100.0	;	# 3 M_maltodextrin4_c+M_pi_c --> M_maltodextrin3_c+M_g1p_c
		0	100.0	;	# 4 M_maltodextrin3_c+M_pi_c --> M_maltose_c+M_g1p_c
		0	100.0	;	# 5 M_maltose_c+M_pi_c --> M_glc_D_c+M_g1p_c
		0	100.0	;	# 6 M_g1p_c --> M_g6p_c
		0	100.0	;	# 7 M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c
		0	100.0	;	# 8 M_g6p_c --> M_f6p_c
		0	100.0	;	# 9 M_f6p_c --> M_g6p_c
		0	100.0	;	# 10 M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c
		0	100.0	;	# 11 M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0	100.0	;	# 12 M_fdp_c --> M_dhap_c+M_g3p_c
		0	100.0	;	# 13 M_dhap_c+M_g3p_c --> M_fdp_c
		0	100.0	;	# 14 M_dhap_c --> M_g3p_c
		0	100.0	;	# 15 M_g3p_c --> M_dhap_c
		0	100.0	;	# 16 M_dhap_c+M_nadph_c+M_h_c --> M_gly3p_c+M_nadp_c
		0	100.0	;	# 17 M_gly3p_c+M_nadp_c --> M_dhap_c+M_nadph_c+M_h_c
		0	100.0	;	# 18 M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0	100.0	;	# 19 M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0	100.0	;	# 20 M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0	100.0	;	# 21 M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0	100.0	;	# 22 M_3pg_c --> M_2pg_c
		0	100.0	;	# 23 M_2pg_c --> M_3pg_c
		0	100.0	;	# 24 M_2pg_c --> M_h2o_c+M_pep_c
		0	100.0	;	# 25 M_h2o_c+M_pep_c --> M_2pg_c
		0	100.0	;	# 26 M_adp_c+M_pep_c --> M_atp_c+M_pyr_c
		0	100.0	;	# 27 M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0	100.0	;	# 28 M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c
		0	100.0	;	# 29 M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c
		0	100.0	;	# 30 M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c
		0	100.0	;	# 31 M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0	100.0	;	# 32 M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0	100.0	;	# 33 M_6pgl_c+M_h2o_c --> M_6pgc_c
		0	100.0	;	# 34 M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
		0	100.0	;	# 35 M_ru5p_D_c --> M_xu5p_D_c
		0	100.0	;	# 36 M_xu5p_D_c --> M_ru5p_D_c
		0	100.0	;	# 37 M_r5p_c --> M_ru5p_D_c
		0	100.0	;	# 38 M_ru5p_D_c --> M_r5p_c
		0	100.0	;	# 39 M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0	100.0	;	# 40 M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 41 M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 42 M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0	100.0	;	# 43 M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0	100.0	;	# 44 M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0	100.0	;	# 45 M_6pgc_c --> M_2ddg6p_c+M_h2o_c
		0	100.0	;	# 46 M_2ddg6p_c --> M_g3p_c+M_pyr_c
		0	100.0	;	# 47 M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c
		0	100.0	;	# 48 M_cit_c --> M_icit_c
		0	100.0	;	# 49 M_icit_c --> M_cit_c
		0	100.0	;	# 50 M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c
		0	100.0	;	# 51 M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c
		0	100.0	;	# 52 M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c
		0	100.0	;	# 53 M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0	100.0	;	# 54 M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0	100.0	;	# 55 M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c
		0	100.0	;	# 56 M_fum_c+M_h2o_c --> M_mal_L_c
		0	100.0	;	# 57 M_mal_L_c --> M_fum_c+M_h2o_c
		0	100.0	;	# 58 M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c
		0	100.0	;	# 59 M_oaa_c+M_h_c+M_nadh_c --> M_mal_L_c+M_nad_c
		0	100.0	;	# 60 2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_hx_c
		0	100.0	;	# 61 4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_hx_c
		0	100.0	;	# 62 2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_hx_c
		0	100.0	;	# 63 M_adp_c+M_pi_c+4.0*M_hx_c --> M_atp_c+4.0*M_h_c+M_h2o_c
		0	100.0	;	# 64 M_hx_c --> M_h_c
		0	100.0	;	# 65 3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_hx_c
		0	100.0	;	# 66 M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0	100.0	;	# 67 M_nadh_c+M_nadp_c+2.0*M_hx_c --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0	100.0	;	# 68 M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c
		0	100.0	;	# 69 M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c
		0	100.0	;	# 70 M_ppi_c+M_h2o_c --> 2.0*M_pi_c
		0	100.0	;	# 71 M_icit_c --> M_glx_c+M_succ_c
		0	100.0	;	# 72 M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c
		0	100.0	;	# 73 M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c
		0	100.0	;	# 74 M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c
		0	100.0	;	# 75 M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0	100.0	;	# 76 M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0	100.0	;	# 77 M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0	100.0	;	# 78 M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0	100.0	;	# 79 M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c
		0	100.0	;	# 80 M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c
		0	100.0	;	# 81 M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c
		0	100.0	;	# 82 M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c
		0	100.0	;	# 83 M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c
		0	100.0	;	# 84 M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0	100.0	;	# 85 M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c
		0	100.0	;	# 86 M_ala_L_c+M_akg_c --> M_pyr_c+M_glu_L_c
		0	100.0	;	# 87 M_accoa_c+2.0*M_glu_L_c+3.0*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c --> M_coa_c+2.0*M_adp_c+2.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c
		0	100.0	;	# 88 M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c
		0	100.0	;	# 89 M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c
		0	100.0	;	# 90 M_asp_L_c+M_atp_c+M_nh3_c --> M_asn_L_c+M_ppi_c+M_amp_c
		0	100.0	;	# 91 M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_ac_c
		0	100.0	;	# 92 M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c
		0	100.0	;	# 93 M_akg_c+M_nadph_c+M_nh3_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0	100.0	;	# 94 M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh3_c+M_h_c
		0	100.0	;	# 95 M_glu_L_c+M_atp_c+M_nh3_c --> M_gln_L_c+M_adp_c+M_pi_c
		0	100.0	;	# 96 M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c
		0	100.0	;	# 97 M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c
		0	100.0	;	# 98 M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c
		0	100.0	;	# 99 2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
		0	100.0	;	# 100 M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
		0	100.0	;	# 101 M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c
		0	100.0	;	# 102 M_chor_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
		0	100.0	;	# 103 M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
		0	100.0	;	# 104 M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
		0	100.0	;	# 105 M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
		0	100.0	;	# 106 M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c
		0	100.0	;	# 107 M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c
		0	100.0	;	# 108 2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
		0	100.0	;	# 109 M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c --> M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c
		0	100.0	;	# 110 M_asp_L_c --> M_fum_c+M_nh3_c
		0	100.0	;	# 111 M_asn_L_c+M_amp_c+M_ppi_c --> M_nh3_c+M_asp_L_c+M_atp_c
		0	100.0	;	# 112 M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c
		0	100.0	;	# 113 M_mglx_c+M_nad_c+M_h2o_c --> M_pyr_c+M_nadh_c+M_h_c
		0	100.0	;	# 114 M_ser_L_c --> M_nh3_c+M_pyr_c
		0	100.0	;	# 115 M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
		0	100.0	;	# 116 M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
		0	100.0	;	# 117 M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c
		0	100.0	;	# 118 M_thr_L_c+M_pi_c+M_adp_c --> M_nh3_c+M_for_c+M_atp_c+M_prop_c
		0	100.0	;	# 119 M_trp_L_c+M_h2o_c --> M_indole_c+M_nh3_c+M_pyr_c
		0	100.0	;	# 120 M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh3_c+M_pyr_c
		0	100.0	;	# 121 M_lys_L_c --> M_co2_c+M_cadav_c
		0	100.0	;	# 122 M_gln_L_c+M_h2o_c --> M_nh3_c+M_glu_L_c
		0	100.0	;	# 123 M_glu_L_c --> M_co2_c+M_gaba_c
		0	100.0	;	# 124 M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadh_c
		0	100.0	;	# 125 M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadph_c
		0	100.0	;	# 126 M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
		0	100.0	;	# 127 M_gtp_c+4.0*M_h2o_c --> M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c
		0	100.0	;	# 128 M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c
		0	100.0	;	# 129 M_4adochor_c --> M_4abz_c+M_pyr_c
		0	100.0	;	# 130 M_4abz_c+M_78mdp_c --> M_78dhf_c+M_h2o_c
		0	100.0	;	# 131 M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c
		0	100.0	;	# 132 M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c
		0	100.0	;	# 133 M_gly_L_c+M_thf_c+M_nad_c --> M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c
		0	100.0	;	# 134 M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c --> M_gly_L_c+M_thf_c+M_nad_c
		0	100.0	;	# 135 M_mlthf_c+M_nadp_c --> M_methf_c+M_nadph_c
		0	100.0	;	# 136 M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c
		0	100.0	;	# 137 M_h2o_c+M_methf_c --> M_10fthf_c+M_h_c
		0	100.0	;	# 138 M_10fthf_c+M_h_c --> M_h2o_c+M_methf_c
		0	100.0	;	# 139 M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c
		0	100.0	;	# 140 M_mlthf_c+M_h_c+M_nadph_c --> M_5mthf_c+M_nadp_c
		0	100.0	;	# 141 M_r5p_c+M_atp_c --> M_prpp_c+M_amp_c
		0	100.0	;	# 142 2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c+2.0*M_h_c
		0	100.0	;	# 143 M_clasp_c+M_asp_L_c+M_q8_c --> M_or_c+M_q8h2_c+M_h2o_c+M_pi_c
		0	100.0	;	# 144 M_prpp_c+M_or_c --> M_omp_c+M_ppi_c
		0	100.0	;	# 145 M_omp_c --> M_ump_c+M_co2_c
		0	100.0	;	# 146 M_utp_c+M_atp_c+M_nh3_c --> M_ctp_c+M_adp_c+M_pi_c
		0	100.0	;	# 147 M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c
		0	100.0	;	# 148 M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c
		0	100.0	;	# 149 M_atp_c+M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c
		0	100.0	;	# 150 M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c
		0	100.0	;	# 151 M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
		0	100.0	;	# 152 M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c
		0	100.0	;	# 153 M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c
		0	100.0	;	# 154 M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c
		0	100.0	;	# 155 M_saicar_c --> M_fum_c+M_aicar_c
		0	100.0	;	# 156 M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c
		0	100.0	;	# 157 M_faicar_c --> M_imp_c+M_h2o_c
		0	100.0	;	# 158 M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c
		0	100.0	;	# 159 M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c
		0	100.0	;	# 160 M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
		0	100.0	;	# 161 M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c
		0	100.0	;	# 162 M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c
		0	100.0	;	# 163 M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c
		0	100.0	;	# 164 M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c
		0	100.0	;	# 165 M_atp_c+M_h2o_c --> M_adp_c+M_pi_c
		0	100.0	;	# 166 M_utp_c+M_h2o_c --> M_udp_c+M_pi_c
		0	100.0	;	# 167 M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c
		0	100.0	;	# 168 M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c
		0	100.0	;	# 169 M_udp_c+M_atp_c --> M_utp_c+M_adp_c
		0	100.0	;	# 170 M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c
		0	100.0	;	# 171 M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c
		0	100.0	;	# 172 M_atp_c+M_ump_c --> M_adp_c+M_udp_c
		0	100.0	;	# 173 M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c
		0	100.0	;	# 174 M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c
		0	100.0	;	# 175 M_amp_c+M_atp_c --> 2.0*M_adp_c
		0	100.0	;	# 176 GENE_gfp+RNAP --> OPEN_GENE_gfp
		0	100.0	;	# 177 OPEN_GENE_gfp+195.0*M_gtp_c+157.0*M_ctp_c+157.0*M_utp_c+208.0*M_atp_c+717.0*M_h2o_c --> mRNA_gfp+GENE_gfp+RNAP+717.0*M_ppi_c
		0	100.0	;	# 178 mRNA_gfp --> 195.0*M_gmp_c+157.0*M_cmp_c+157.0*M_ump_c+208.0*M_amp_c
		0	100.0	;	# 179 mRNA_gfp+RIBOSOME --> RIBOSOME_START_gfp
		0	100.0	;	# 180 RIBOSOME_START_gfp+478.0*M_gtp_c+478.0*M_h2o_c+8.0*M_ala_L_c_tRNA_c+6.0*M_arg_L_c_tRNA_c+13.0*M_asn_L_c_tRNA_c+18.0*M_asp_L_c_tRNA_c+2.0*M_cys_L_c_tRNA_c+16.0*M_glu_L_c_tRNA_c+8.0*M_gln_L_c_tRNA_c+22.0*M_gly_L_c_tRNA_c+9.0*M_his_L_c_tRNA_c+12.0*M_ile_L_c_tRNA_c+21.0*M_leu_L_c_tRNA_c+20.0*M_lys_L_c_tRNA_c+6.0*M_met_L_c_tRNA_c+12.0*M_phe_L_c_tRNA_c+10.0*M_pro_L_c_tRNA_c+10.0*M_ser_L_c_tRNA_c+16.0*M_thr_L_c_tRNA_c+M_trp_L_c_tRNA_c+11.0*M_tyr_L_c_tRNA_c+18.0*M_val_L_c_tRNA_c --> RIBOSOME+mRNA_gfp+PROTEIN_gfp+478.0*M_gdp_c+478.0*M_pi_c+239.0*tRNA_c
		0	100.0	;	# 181 8.0*M_ala_L_c+8.0*M_atp_c+8.0*tRNA_c+8.0*M_h2o_c --> 8.0*M_ala_L_c_tRNA_c+8.0*M_amp_c+8.0*M_ppi_c
		0	100.0	;	# 182 6.0*M_arg_L_c+6.0*M_atp_c+6.0*tRNA_c+6.0*M_h2o_c --> 6.0*M_arg_L_c_tRNA_c+6.0*M_amp_c+6.0*M_ppi_c
		0	100.0	;	# 183 13.0*M_asn_L_c+13.0*M_atp_c+13.0*tRNA_c+13.0*M_h2o_c --> 13.0*M_asn_L_c_tRNA_c+13.0*M_amp_c+13.0*M_ppi_c
		0	100.0	;	# 184 18.0*M_asp_L_c+18.0*M_atp_c+18.0*tRNA_c+18.0*M_h2o_c --> 18.0*M_asp_L_c_tRNA_c+18.0*M_amp_c+18.0*M_ppi_c
		0	100.0	;	# 185 2.0*M_cys_L_c+2.0*M_atp_c+2.0*tRNA_c+2.0*M_h2o_c --> 2.0*M_cys_L_c_tRNA_c+2.0*M_amp_c+2.0*M_ppi_c
		0	100.0	;	# 186 16.0*M_glu_L_c+16.0*M_atp_c+16.0*tRNA_c+16.0*M_h2o_c --> 16.0*M_glu_L_c_tRNA_c+16.0*M_amp_c+16.0*M_ppi_c
		0	100.0	;	# 187 8.0*M_gln_L_c+8.0*M_atp_c+8.0*tRNA_c+8.0*M_h2o_c --> 8.0*M_gln_L_c_tRNA_c+8.0*M_amp_c+8.0*M_ppi_c
		0	100.0	;	# 188 22.0*M_gly_L_c+22.0*M_atp_c+22.0*tRNA_c+22.0*M_h2o_c --> 22.0*M_gly_L_c_tRNA_c+22.0*M_amp_c+22.0*M_ppi_c
		0	100.0	;	# 189 9.0*M_his_L_c+9.0*M_atp_c+9.0*tRNA_c+9.0*M_h2o_c --> 9.0*M_his_L_c_tRNA_c+9.0*M_amp_c+9.0*M_ppi_c
		0	100.0	;	# 190 12.0*M_ile_L_c+12.0*M_atp_c+12.0*tRNA_c+12.0*M_h2o_c --> 12.0*M_ile_L_c_tRNA_c+12.0*M_amp_c+12.0*M_ppi_c
		0	100.0	;	# 191 21.0*M_leu_L_c+21.0*M_atp_c+21.0*tRNA_c+21.0*M_h2o_c --> 21.0*M_leu_L_c_tRNA_c+21.0*M_amp_c+21.0*M_ppi_c
		0	100.0	;	# 192 20.0*M_lys_L_c+20.0*M_atp_c+20.0*tRNA_c+20.0*M_h2o_c --> 20.0*M_lys_L_c_tRNA_c+20.0*M_amp_c+20.0*M_ppi_c
		0	100.0	;	# 193 6.0*M_met_L_c+6.0*M_atp_c+6.0*tRNA_c+6.0*M_h2o_c --> 6.0*M_met_L_c_tRNA_c+6.0*M_amp_c+6.0*M_ppi_c
		0	100.0	;	# 194 12.0*M_phe_L_c+12.0*M_atp_c+12.0*tRNA_c+12.0*M_h2o_c --> 12.0*M_phe_L_c_tRNA_c+12.0*M_amp_c+12.0*M_ppi_c
		0	100.0	;	# 195 10.0*M_pro_L_c+10.0*M_atp_c+10.0*tRNA_c+10.0*M_h2o_c --> 10.0*M_pro_L_c_tRNA_c+10.0*M_amp_c+10.0*M_ppi_c
		0	100.0	;	# 196 10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA_c+10.0*M_h2o_c --> 10.0*M_ser_L_c_tRNA_c+10.0*M_amp_c+10.0*M_ppi_c
		0	100.0	;	# 197 16.0*M_thr_L_c+16.0*M_atp_c+16.0*tRNA_c+16.0*M_h2o_c --> 16.0*M_thr_L_c_tRNA_c+16.0*M_amp_c+16.0*M_ppi_c
		0	100.0	;	# 198 M_trp_L_c+M_atp_c+tRNA_c+M_h2o_c --> M_trp_L_c_tRNA_c+M_amp_c+M_ppi_c
		0	100.0	;	# 199 11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA_c+11.0*M_h2o_c --> 11.0*M_tyr_L_c_tRNA_c+11.0*M_amp_c+11.0*M_ppi_c
		0	100.0	;	# 200 18.0*M_val_L_c+18.0*M_atp_c+18.0*tRNA_c+18.0*M_h2o_c --> 18.0*M_val_L_c_tRNA_c+18.0*M_amp_c+18.0*M_ppi_c
		0	100.0	;	# 201 [] --> tRNA_c
		0	100.0	;	# 202 tRNA_c --> []
		0	100.0	;	# 203 PROTEIN_gfp --> gfp_e
		0	100.0	;	# 204 [] --> M_h2s_c
		0	100.0	;	# 205 M_h2s_c --> []
		0	100.0	;	# 206 [] --> M_h2o_c
		0	100.0	;	# 207 M_h2o_c --> []
		0	100.0	;	# 208 [] --> M_pi_c
		0	100.0	;	# 209 M_pi_c --> []
		0	100.0	;	# 210 [] --> M_nh3_c
		0	100.0	;	# 211 M_nh3_c --> []
		0	100.0	;	# 212 [] --> M_hco3_c
		0	100.0	;	# 213 M_etoh_c --> []
		0	100.0	;	# 214 M_mglx_c --> []
		0	100.0	;	# 215 M_prop_c --> []
		0	100.0	;	# 216 M_indole_c --> []
		0	100.0	;	# 217 M_cadav_c --> []
		0	100.0	;	# 218 M_gaba_c --> []
		0	100.0	;	# 219 M_glycoA_c --> []
		0	100.0	;	# 220 M_o2_e --> M_o2_c
		0	100.0	;	# 221 M_co2_c --> M_co2_e
		0	100.0	;	# 222 M_h_c --> M_h_e
		0	100.0	;	# 223 M_h_e --> M_h_c
	];

	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1 GENE_gfp
		0.0	0.0	;	# 2 M_10fthf_c
		0.0	0.0	;	# 3 M_13dpg_c
		0.0	0.0	;	# 4 M_2ddg6p_c
		0.0	0.0	;	# 5 M_2pg_c
		0.0	0.0	;	# 6 M_3pg_c
		0.0	0.0	;	# 7 M_4abz_c
		0.0	0.0	;	# 8 M_4adochor_c
		0.0	0.0	;	# 9 M_5mthf_c
		0.0	0.0	;	# 10 M_5pbdra
		0.0	0.0	;	# 11 M_6pgc_c
		0.0	0.0	;	# 12 M_6pgl_c
		0.0	0.0	;	# 13 M_78dhf_c
		0.0	0.0	;	# 14 M_78mdp_c
		0.0	0.0	;	# 15 M_ac_c
		0.0	0.0	;	# 16 M_accoa_c
		0.0	0.0	;	# 17 M_actp_c
		0.0	0.0	;	# 18 M_adp_c
		0.0	0.0	;	# 19 M_aicar_c
		0.0	0.0	;	# 20 M_air_c
		0.0	0.0	;	# 21 M_akg_c
		0.0	0.0	;	# 22 M_ala_L_c
		0.0	0.0	;	# 23 M_ala_L_c_tRNA_c
		0.0	0.0	;	# 24 M_amp_c
		0.0	0.0	;	# 25 M_arg_L_c
		0.0	0.0	;	# 26 M_arg_L_c_tRNA_c
		0.0	0.0	;	# 27 M_asn_L_c
		0.0	0.0	;	# 28 M_asn_L_c_tRNA_c
		0.0	0.0	;	# 29 M_asp_L_c
		0.0	0.0	;	# 30 M_asp_L_c_tRNA_c
		0.0	0.0	;	# 31 M_atp_c
		0.0	0.0	;	# 32 M_cadav_c
		0.0	0.0	;	# 33 M_cair_c
		0.0	0.0	;	# 34 M_cdp_c
		0.0	0.0	;	# 35 M_chor_c
		0.0	0.0	;	# 36 M_cit_c
		0.0	0.0	;	# 37 M_clasp_c
		0.0	0.0	;	# 38 M_cmp_c
		0.0	0.0	;	# 39 M_co2_c
		0.0	0.0	;	# 40 M_coa_c
		0.0	0.0	;	# 41 M_ctp_c
		0.0	0.0	;	# 42 M_cys_L_c
		0.0	0.0	;	# 43 M_cys_L_c_tRNA_c
		0.0	0.0	;	# 44 M_dhap_c
		0.0	0.0	;	# 45 M_dhf_c
		0.0	0.0	;	# 46 M_e4p_c
		0.0	0.0	;	# 47 M_etoh_c
		0.0	0.0	;	# 48 M_f6p_c
		0.0	0.0	;	# 49 M_faicar_c
		0.0	0.0	;	# 50 M_fdp_c
		0.0	0.0	;	# 51 M_fgam_c
		0.0	0.0	;	# 52 M_fgar_c
		0.0	0.0	;	# 53 M_for_c
		0.0	0.0	;	# 54 M_fum_c
		0.0	0.0	;	# 55 M_g1p_c
		0.0	0.0	;	# 56 M_g3p_c
		0.0	0.0	;	# 57 M_g6p_c
		0.0	0.0	;	# 58 M_gaba_c
		0.0	0.0	;	# 59 M_gar_c
		0.0	0.0	;	# 60 M_gdp_c
		0.0	0.0	;	# 61 M_glc_D_c
		0.0	0.0	;	# 62 M_gln_L_c
		0.0	0.0	;	# 63 M_gln_L_c_tRNA_c
		0.0	0.0	;	# 64 M_glu_L_c
		0.0	0.0	;	# 65 M_glu_L_c_tRNA_c
		0.0	0.0	;	# 66 M_glx_c
		0.0	0.0	;	# 67 M_gly3p_c
		0.0	0.0	;	# 68 M_gly_L_c
		0.0	0.0	;	# 69 M_gly_L_c_tRNA_c
		0.0	0.0	;	# 70 M_glycoA_c
		0.0	0.0	;	# 71 M_gmp_c
		0.0	0.0	;	# 72 M_gtp_c
		0.0	0.0	;	# 73 M_h2o2_c
		0.0	0.0	;	# 74 M_h2o_c
		0.0	0.0	;	# 75 M_h2s_c
		0.0	0.0	;	# 76 M_h_c
		0.0	0.0	;	# 77 M_hco3_c
		0.0	0.0	;	# 78 M_his_L_c
		0.0	0.0	;	# 79 M_his_L_c_tRNA_c
		0.0	0.0	;	# 80 M_hx_c
		0.0	0.0	;	# 81 M_icit_c
		0.0	0.0	;	# 82 M_ile_L_c
		0.0	0.0	;	# 83 M_ile_L_c_tRNA_c
		0.0	0.0	;	# 84 M_imp_c
		0.0	0.0	;	# 85 M_indole_c
		0.0	0.0	;	# 86 M_lac_D_c
		0.0	0.0	;	# 87 M_leu_L_c
		0.0	0.0	;	# 88 M_leu_L_c_tRNA_c
		0.0	0.0	;	# 89 M_lys_L_c
		0.0	0.0	;	# 90 M_lys_L_c_tRNA_c
		0.0	0.0	;	# 91 M_mal_L_c
		0.0	0.0	;	# 92 M_maltodextrin3_c
		0.0	0.0	;	# 93 M_maltodextrin4_c
		0.0	0.0	;	# 94 M_maltodextrin5_c
		0.0	0.0	;	# 95 M_maltodextrin6_c
		0.0	0.0	;	# 96 M_maltose_c
		0.0	0.0	;	# 97 M_met_L_c
		0.0	0.0	;	# 98 M_met_L_c_tRNA_c
		0.0	0.0	;	# 99 M_methf_c
		0.0	0.0	;	# 100 M_mglx_c
		0.0	0.0	;	# 101 M_mlthf_c
		0.0	0.0	;	# 102 M_mql8_c
		0.0	0.0	;	# 103 M_mqn8_c
		0.0	0.0	;	# 104 M_nad_c
		0.0	0.0	;	# 105 M_nadh_c
		0.0	0.0	;	# 106 M_nadp_c
		0.0	0.0	;	# 107 M_nadph_c
		0.0	0.0	;	# 108 M_nh3_c
		0.0	0.0	;	# 109 M_o2_c
		0.0	0.0	;	# 110 M_oaa_c
		0.0	0.0	;	# 111 M_omp_c
		0.0	0.0	;	# 112 M_or_c
		0.0	0.0	;	# 113 M_pep_c
		0.0	0.0	;	# 114 M_phe_L_c
		0.0	0.0	;	# 115 M_phe_L_c_tRNA_c
		0.0	0.0	;	# 116 M_pi_c
		0.0	0.0	;	# 117 M_ppi_c
		0.0	0.0	;	# 118 M_pro_L_c
		0.0	0.0	;	# 119 M_pro_L_c_tRNA_c
		0.0	0.0	;	# 120 M_prop_c
		0.0	0.0	;	# 121 M_prpp_c
		0.0	0.0	;	# 122 M_pyr_c
		0.0	0.0	;	# 123 M_q8_c
		0.0	0.0	;	# 124 M_q8h2_c
		0.0	0.0	;	# 125 M_r5p_c
		0.0	0.0	;	# 126 M_ru5p_D_c
		0.0	0.0	;	# 127 M_s7p_c
		0.0	0.0	;	# 128 M_saicar_c
		0.0	0.0	;	# 129 M_ser_L_c
		0.0	0.0	;	# 130 M_ser_L_c_tRNA_c
		0.0	0.0	;	# 131 M_succ_c
		0.0	0.0	;	# 132 M_succoa_c
		0.0	0.0	;	# 133 M_thf_c
		0.0	0.0	;	# 134 M_thr_L_c
		0.0	0.0	;	# 135 M_thr_L_c_tRNA_c
		0.0	0.0	;	# 136 M_trp_L_c
		0.0	0.0	;	# 137 M_trp_L_c_tRNA_c
		0.0	0.0	;	# 138 M_tyr_L_c
		0.0	0.0	;	# 139 M_tyr_L_c_tRNA_c
		0.0	0.0	;	# 140 M_udp_c
		0.0	0.0	;	# 141 M_ump_c
		0.0	0.0	;	# 142 M_utp_c
		0.0	0.0	;	# 143 M_val_L_c
		0.0	0.0	;	# 144 M_val_L_c_tRNA_c
		0.0	0.0	;	# 145 M_xmp_c
		0.0	0.0	;	# 146 M_xu5p_D_c
		0.0	0.0	;	# 147 OPEN_GENE_gfp
		0.0	0.0	;	# 148 PROTEIN_gfp
		0.0	0.0	;	# 149 RIBOSOME
		0.0	0.0	;	# 150 RIBOSOME_START_gfp
		0.0	0.0	;	# 151 RNAP
		0.0	0.0	;	# 152 mRNA_gfp
		0.0	0.0	;	# 153 tRNA_c
		-1.0	1.0	;	# 154 M_co2_e
		-1.0	1.0	;	# 155 M_h_e
		-1.0	1.0	;	# 156 M_o2_e
		-1.0	1.0	;	# 157 gfp_e
	];

	# Setup the objective coefficient•••••••••
	objective_coefficient_array = [
		0.0	;	# 1 R_gp_1::M_maltode•••••••••> M_maltodextrin5_c+M_g1p_c
		0.0	;	# 2 R_gp_2::M_maltodextrin5_c+M_pi_c --> M_maltodextrin4_c+M_g1p_c
		0.0	;	# 3 R_gp_3::M_maltodextrin4_c+M_pi_c --> M_maltodextrin3_c+M_g1p_c
		0.0	;	# 4 R_gp_4::M_maltodextrin3_c+M_pi_c --> M_maltose_c+M_g1p_c
		0.0	;	# 5 R_gp::M_maltose_c+M_pi_c --> M_glc_D_c+M_g1p_c
		0.0	;	# 6 R_pgm::M_g1p_c --> M_g6p_c
		0.0	;	# 7 R_glk_atp::M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c
		0.0	;	# 8 R_pgi::M_g6p_c --> M_f6p_c
		0.0	;	# 9 R_pgi_reverse::M_f6p_c --> M_g6p_c
		0.0	;	# 10 R_pfk::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c
		0.0	;	# 11 R_fdp::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0.0	;	# 12 R_fbaA::M_fdp_c --> M_dhap_c+M_g3p_c
		0.0	;	# 13 R_fbaA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c
		0.0	;	# 14 R_tpiA::M_dhap_c --> M_g3p_c
		0.0	;	# 15 R_tpiA_reverse::M_g3p_c --> M_dhap_c
		0.0	;	# 16 R_gpsA::M_dhap_c+M_nadph_c+M_h_c --> M_gly3p_c+M_nadp_c
		0.0	;	# 17 R_gpsA_reverse::M_gly3p_c+M_nadp_c --> M_dhap_c+M_nadph_c+M_h_c
		0.0	;	# 18 R_gapA::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0.0	;	# 19 R_gapA_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0.0	;	# 20 R_pgk::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0.0	;	# 21 R_pgk_reverse::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0.0	;	# 22 R_gpm::M_3pg_c --> M_2pg_c
		0.0	;	# 23 R_gpm_reverse::M_2pg_c --> M_3pg_c
		0.0	;	# 24 R_eno::M_2pg_c --> M_h2o_c+M_pep_c
		0.0	;	# 25 R_eno_reverse::M_h2o_c+M_pep_c --> M_2pg_c
		0.0	;	# 26 R_pyk::M_adp_c+M_pep_c --> M_atp_c+M_pyr_c
		0.0	;	# 27 R_pck::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0.0	;	# 28 R_ppc::M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c
		0.0	;	# 29 R_pdh::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c
		0.0	;	# 30 R_pps::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c
		0.0	;	# 31 R_zwf::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0.0	;	# 32 R_zwf_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0.0	;	# 33 R_pgl::M_6pgl_c+M_h2o_c --> M_6pgc_c
		0.0	;	# 34 R_gnd::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
		0.0	;	# 35 R_rpe::M_ru5p_D_c --> M_xu5p_D_c
		0.0	;	# 36 R_rpe_reverse::M_xu5p_D_c --> M_ru5p_D_c
		0.0	;	# 37 R_rpi::M_r5p_c --> M_ru5p_D_c
		0.0	;	# 38 R_rpi_reverse::M_ru5p_D_c --> M_r5p_c
		0.0	;	# 39 R_talAB::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0.0	;	# 40 R_talAB_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0.0	;	# 41 R_tkt1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0.0	;	# 42 R_tkt1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0.0	;	# 43 R_tkt2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0.0	;	# 44 R_tkt2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0.0	;	# 45 R_edd::M_6pgc_c --> M_2ddg6p_c+M_h2o_c
		0.0	;	# 46 R_eda::M_2ddg6p_c --> M_g3p_c+M_pyr_c
		0.0	;	# 47 R_gltA::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c
		0.0	;	# 48 R_acn::M_cit_c --> M_icit_c
		0.0	;	# 49 R_acn_reverse::M_icit_c --> M_cit_c
		0.0	;	# 50 R_icd::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c
		0.0	;	# 51 R_icd_reverse::M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c
		0.0	;	# 52 R_sucAB::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c
		0.0	;	# 53 R_sucCD::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0.0	;	# 54 R_sdh::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0.0	;	# 55 R_frd::M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c
		0.0	;	# 56 R_fum::M_fum_c+M_h2o_c --> M_mal_L_c
		0.0	;	# 57 R_fum_reverse::M_mal_L_c --> M_fum_c+M_h2o_c
		0.0	;	# 58 R_mdh::M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c
		0.0	;	# 59 R_mdh_reverse::M_oaa_c+M_h_c+M_nadh_c --> M_mal_L_c+M_nad_c
		0.0	;	# 60 R_cyd::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_hx_c
		0.0	;	# 61 R_cyo::4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_hx_c
		0.0	;	# 62 R_app::2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_hx_c
		0.0	;	# 63 R_atp::M_adp_c+M_pi_c+4.0*M_hx_c --> M_atp_c+4.0*M_h_c+M_h2o_c
		0.0	;	# 64 R_DNP_leak::M_hx_c --> M_h_c
		0.0	;	# 65 R_nuo::3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_hx_c
		0.0	;	# 66 R_pnt1::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0.0	;	# 67 R_pnt2::M_nadh_c+M_nadp_c+2.0*M_hx_c --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0.0	;	# 68 R_ndh1::M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c
		0.0	;	# 69 R_ndh2::M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c
		0.0	;	# 70 R_ppa::M_ppi_c+M_h2o_c --> 2.0*M_pi_c
		0.0	;	# 71 R_aceA::M_icit_c --> M_glx_c+M_succ_c
		0.0	;	# 72 R_aceB::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c
		0.0	;	# 73 R_maeA::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c
		0.0	;	# 74 R_maeB::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c
		0.0	;	# 75 R_pta::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0.0	;	# 76 R_pta_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0.0	;	# 77 R_ackA::M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0.0	;	# 78 R_ackA_reverse::M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0.0	;	# 79 R_acs::M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c
		0.0	;	# 80 R_adhE::M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c
		0.0	;	# 81 R_adhE_reverse::M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c
		0.0	;	# 82 R_ldh::M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c
		0.0	;	# 83 R_ldh_reverse::M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c
		0.0	;	# 84 R_pflAB::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0.0	;	# 85 R_alaAC::M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c
		0.0	;	# 86 R_alaAC_reverse::M_ala_L_c+M_akg_c --> M_pyr_c+M_glu_L_c
		0.0	;	# 87 R_arg::M_accoa_c+2.0*M_glu_L_c+3.0*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c --> M_coa_c+2.0*M_adp_c+2.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c
		0.0	;	# 88 R_aspC::M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c
		0.0	;	# 89 R_asnB::M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c
		0.0	;	# 90 R_asnA::M_asp_L_c+M_atp_c+M_nh3_c --> M_asn_L_c+M_ppi_c+M_amp_c
		0.0	;	# 91 R_cysEMK::M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_ac_c
		0.0	;	# 92 R_gltBD::M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c
		0.0	;	# 93 R_gdhA::M_akg_c+M_nadph_c+M_nh3_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0.0	;	# 94 R_gdhA_reverse::M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh3_c+M_h_c
		0.0	;	# 95 R_glnA::M_glu_L_c+M_atp_c+M_nh3_c --> M_gln_L_c+M_adp_c+M_pi_c
		0.0	;	# 96 R_glyA::M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c
		0.0	;	# 97 R_his::M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c
		0.0	;	# 98 R_ile::M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c
		0.0	;	# 99 R_leu::2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
		0.0	;	# 100 R_lys::M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
		0.0	;	# 101 R_met::M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c
		0.0	;	# 102 R_phe::M_chor_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
		0.0	;	# 103 R_pro::M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
		0.0	;	# 104 R_serABC::M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
		0.0	;	# 105 R_thr::M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
		0.0	;	# 106 R_trp::M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c
		0.0	;	# 107 R_tyr::M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c
		0.0	;	# 108 R_val::2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
		0.0	;	# 109 R_arg_deg::M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c --> M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c
		0.0	;	# 110 R_asp_deg::M_asp_L_c --> M_fum_c+M_nh3_c
		0.0	;	# 111 R_asn_deg::M_asn_L_c+M_amp_c+M_ppi_c --> M_nh3_c+M_asp_L_c+M_atp_c
		0.0	;	# 112 R_gly_deg::M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c
		0.0	;	# 113 R_mglx_deg::M_mglx_c+M_nad_c+M_h2o_c --> M_pyr_c+M_nadh_c+M_h_c
		0.0	;	# 114 R_ser_deg::M_ser_L_c --> M_nh3_c+M_pyr_c
		0.0	;	# 115 R_pro_deg::M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
		0.0	;	# 116 R_thr_deg1::M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
		0.0	;	# 117 R_thr_deg2::M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c
		0.0	;	# 118 R_thr_deg3::M_thr_L_c+M_pi_c+M_adp_c --> M_nh3_c+M_for_c+M_atp_c+M_prop_c
		0.0	;	# 119 R_trp_deg::M_trp_L_c+M_h2o_c --> M_indole_c+M_nh3_c+M_pyr_c
		0.0	;	# 120 R_cys_deg::M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh3_c+M_pyr_c
		0.0	;	# 121 R_lys_deg::M_lys_L_c --> M_co2_c+M_cadav_c
		0.0	;	# 122 R_gln_deg::M_gln_L_c+M_h2o_c --> M_nh3_c+M_glu_L_c
		0.0	;	# 123 R_glu_deg::M_glu_L_c --> M_co2_c+M_gaba_c
		0.0	;	# 124 R_gaba_deg1::M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadh_c
		0.0	;	# 125 R_gaba_deg2::M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadph_c
		0.0	;	# 126 R_chor::M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
		0.0	;	# 127 R_fol_e::M_gtp_c+4.0*M_h2o_c --> M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c
		0.0	;	# 128 R_fol_1::M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c
		0.0	;	# 129 R_fol_2a::M_4adochor_c --> M_4abz_c+M_pyr_c
		0.0	;	# 130 R_fol_2b::M_4abz_c+M_78mdp_c --> M_78dhf_c+M_h2o_c
		0.0	;	# 131 R_fol_3::M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c
		0.0	;	# 132 R_fol_4::M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c
		0.0	;	# 133 R_gly_fol::M_gly_L_c+M_thf_c+M_nad_c --> M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c
		0.0	;	# 134 R_gly_fol_reverse::M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c --> M_gly_L_c+M_thf_c+M_nad_c
		0.0	;	# 135 R_mthfd::M_mlthf_c+M_nadp_c --> M_methf_c+M_nadph_c
		0.0	;	# 136 R_mthfd_reverse::M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c
		0.0	;	# 137 R_mthfc::M_h2o_c+M_methf_c --> M_10fthf_c+M_h_c
		0.0	;	# 138 R_mthfc_reverse::M_10fthf_c+M_h_c --> M_h2o_c+M_methf_c
		0.0	;	# 139 R_mthfr2a::M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c
		0.0	;	# 140 R_mthfr2b::M_mlthf_c+M_h_c+M_nadph_c --> M_5mthf_c+M_nadp_c
		0.0	;	# 141 R_prpp_syn::M_r5p_c+M_atp_c --> M_prpp_c+M_amp_c
		0.0	;	# 142 R_or_syn_1::2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c+2.0*M_h_c
		0.0	;	# 143 R_or_syn_2::M_clasp_c+M_asp_L_c+M_q8_c --> M_or_c+M_q8h2_c+M_h2o_c+M_pi_c
		0.0	;	# 144 R_omp_syn::M_prpp_c+M_or_c --> M_omp_c+M_ppi_c
		0.0	;	# 145 R_ump_syn::M_omp_c --> M_ump_c+M_co2_c
		0.0	;	# 146 R_ctp_1::M_utp_c+M_atp_c+M_nh3_c --> M_ctp_c+M_adp_c+M_pi_c
		0.0	;	# 147 R_ctp_2::M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c
		0.0	;	# 148 R_A_syn_1::M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c
		0.0	;	# 149 R_A_syn_2::M_atp_c+M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c
		0.0	;	# 150 R_A_syn_3::M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c
		0.0	;	# 151 R_A_syn_4::M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
		0.0	;	# 152 R_A_syn_5::M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c
		0.0	;	# 153 R_A_syn_6::M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c
		0.0	;	# 154 R_A_syn_7::M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c
		0.0	;	# 155 R_A_syn_8::M_saicar_c --> M_fum_c+M_aicar_c
		0.0	;	# 156 R_A_syn_9::M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c
		0.0	;	# 157 R_A_syn_10::M_faicar_c --> M_imp_c+M_h2o_c
		0.0	;	# 158 R_A_syn_12::M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c
		0.0	;	# 159 R_xmp_syn::M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c
		0.0	;	# 160 R_gmp_syn::M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
		0.0	;	# 161 R_atp_amp::M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c
		0.0	;	# 162 R_utp_ump::M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c
		0.0	;	# 163 R_ctp_cmp::M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c
		0.0	;	# 164 R_gtp_gmp::M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c
		0.0	;	# 165 R_atp_adp::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c
		0.0	;	# 166 R_utp_adp::M_utp_c+M_h2o_c --> M_udp_c+M_pi_c
		0.0	;	# 167 R_ctp_adp::M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c
		0.0	;	# 168 R_gtp_adp::M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c
		0.0	;	# 169 R_udp_utp::M_udp_c+M_atp_c --> M_utp_c+M_adp_c
		0.0	;	# 170 R_cdp_ctp::M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c
		0.0	;	# 171 R_gdp_gtp::M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c
		0.0	;	# 172 R_atp_ump::M_atp_c+M_ump_c --> M_adp_c+M_udp_c
		0.0	;	# 173 R_atp_cmp::M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c
		0.0	;	# 174 R_atp_gmp::M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c
		0.0	;	# 175 R_adk_atp::M_amp_c+M_atp_c --> 2.0*M_adp_c
		0.0	;	# 176 transcriptional_initiation_gfp::GENE_gfp+RNAP --> OPEN_GENE_gfp
		0.0	;	# 177 transcription_gfp::OPEN_GENE_gfp+195.0*M_gtp_c+157.0*M_ctp_c+157.0*M_utp_c+208.0*M_atp_c+717.0*M_h2o_c --> mRNA_gfp+GENE_gfp+RNAP+717.0*M_ppi_c
		0.0	;	# 178 mRNA_degradation_gfp::mRNA_gfp --> 195.0*M_gmp_c+157.0*M_cmp_c+157.0*M_ump_c+208.0*M_amp_c
		0.0	;	# 179 translation_initiation_gfp::mRNA_gfp+RIBOSOME --> RIBOSOME_START_gfp
		0.0	;	# 180 translation_gfp::RIBOSOME_START_gfp+478.0*M_gtp_c+478.0*M_h2o_c+8.0*M_ala_L_c_tRNA_c+6.0*M_arg_L_c_tRNA_c+13.0*M_asn_L_c_tRNA_c+18.0*M_asp_L_c_tRNA_c+2.0*M_cys_L_c_tRNA_c+16.0*M_glu_L_c_tRNA_c+8.0*M_gln_L_c_tRNA_c+22.0*M_gly_L_c_tRNA_c+9.0*M_his_L_c_tRNA_c+12.0*M_ile_L_c_tRNA_c+21.0*M_leu_L_c_tRNA_c+20.0*M_lys_L_c_tRNA_c+6.0*M_met_L_c_tRNA_c+12.0*M_phe_L_c_tRNA_c+10.0*M_pro_L_c_tRNA_c+10.0*M_ser_L_c_tRNA_c+16.0*M_thr_L_c_tRNA_c+M_trp_L_c_tRNA_c+11.0*M_tyr_L_c_tRNA_c+18.0*M_val_L_c_tRNA_c --> RIBOSOME+mRNA_gfp+PROTEIN_gfp+478.0*M_gdp_c+478.0*M_pi_c+239.0*tRNA_c
		0.0	;	# 181 tRNA_c_charging_M_ala_L_c_gfp::8.0*M_ala_L_c+8.0*M_atp_c+8.0*tRNA_c+8.0*M_h2o_c --> 8.0*M_ala_L_c_tRNA_c+8.0*M_amp_c+8.0*M_ppi_c
		0.0	;	# 182 tRNA_c_charging_M_arg_L_c_gfp::6.0*M_arg_L_c+6.0*M_atp_c+6.0*tRNA_c+6.0*M_h2o_c --> 6.0*M_arg_L_c_tRNA_c+6.0*M_amp_c+6.0*M_ppi_c
		0.0	;	# 183 tRNA_c_charging_M_asn_L_c_gfp::13.0*M_asn_L_c+13.0*M_atp_c+13.0*tRNA_c+13.0*M_h2o_c --> 13.0*M_asn_L_c_tRNA_c+13.0*M_amp_c+13.0*M_ppi_c
		0.0	;	# 184 tRNA_c_charging_M_asp_L_c_gfp::18.0*M_asp_L_c+18.0*M_atp_c+18.0*tRNA_c+18.0*M_h2o_c --> 18.0*M_asp_L_c_tRNA_c+18.0*M_amp_c+18.0*M_ppi_c
		0.0	;	# 185 tRNA_c_charging_M_cys_L_c_gfp::2.0*M_cys_L_c+2.0*M_atp_c+2.0*tRNA_c+2.0*M_h2o_c --> 2.0*M_cys_L_c_tRNA_c+2.0*M_amp_c+2.0*M_ppi_c
		0.0	;	# 186 tRNA_c_charging_M_glu_L_c_gfp::16.0*M_glu_L_c+16.0*M_atp_c+16.0*tRNA_c+16.0*M_h2o_c --> 16.0*M_glu_L_c_tRNA_c+16.0*M_amp_c+16.0*M_ppi_c
		0.0	;	# 187 tRNA_c_charging_M_gln_L_c_gfp::8.0*M_gln_L_c+8.0*M_atp_c+8.0*tRNA_c+8.0*M_h2o_c --> 8.0*M_gln_L_c_tRNA_c+8.0*M_amp_c+8.0*M_ppi_c
		0.0	;	# 188 tRNA_c_charging_M_gly_L_c_gfp::22.0*M_gly_L_c+22.0*M_atp_c+22.0*tRNA_c+22.0*M_h2o_c --> 22.0*M_gly_L_c_tRNA_c+22.0*M_amp_c+22.0*M_ppi_c
		0.0	;	# 189 tRNA_c_charging_M_his_L_c_gfp::9.0*M_his_L_c+9.0*M_atp_c+9.0*tRNA_c+9.0*M_h2o_c --> 9.0*M_his_L_c_tRNA_c+9.0*M_amp_c+9.0*M_ppi_c
		0.0	;	# 190 tRNA_c_charging_M_ile_L_c_gfp::12.0*M_ile_L_c+12.0*M_atp_c+12.0*tRNA_c+12.0*M_h2o_c --> 12.0*M_ile_L_c_tRNA_c+12.0*M_amp_c+12.0*M_ppi_c
		0.0	;	# 191 tRNA_c_charging_M_leu_L_c_gfp::21.0*M_leu_L_c+21.0*M_atp_c+21.0*tRNA_c+21.0*M_h2o_c --> 21.0*M_leu_L_c_tRNA_c+21.0*M_amp_c+21.0*M_ppi_c
		0.0	;	# 192 tRNA_c_charging_M_lys_L_c_gfp::20.0*M_lys_L_c+20.0*M_atp_c+20.0*tRNA_c+20.0*M_h2o_c --> 20.0*M_lys_L_c_tRNA_c+20.0*M_amp_c+20.0*M_ppi_c
		0.0	;	# 193 tRNA_c_charging_M_met_L_c_gfp::6.0*M_met_L_c+6.0*M_atp_c+6.0*tRNA_c+6.0*M_h2o_c --> 6.0*M_met_L_c_tRNA_c+6.0*M_amp_c+6.0*M_ppi_c
		0.0	;	# 194 tRNA_c_charging_M_phe_L_c_gfp::12.0*M_phe_L_c+12.0*M_atp_c+12.0*tRNA_c+12.0*M_h2o_c --> 12.0*M_phe_L_c_tRNA_c+12.0*M_amp_c+12.0*M_ppi_c
		0.0	;	# 195 tRNA_c_charging_M_pro_L_c_gfp::10.0*M_pro_L_c+10.0*M_atp_c+10.0*tRNA_c+10.0*M_h2o_c --> 10.0*M_pro_L_c_tRNA_c+10.0*M_amp_c+10.0*M_ppi_c
		0.0	;	# 196 tRNA_c_charging_M_ser_L_c_gfp::10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA_c+10.0*M_h2o_c --> 10.0*M_ser_L_c_tRNA_c+10.0*M_amp_c+10.0*M_ppi_c
		0.0	;	# 197 tRNA_c_charging_M_thr_L_c_gfp::16.0*M_thr_L_c+16.0*M_atp_c+16.0*tRNA_c+16.0*M_h2o_c --> 16.0*M_thr_L_c_tRNA_c+16.0*M_amp_c+16.0*M_ppi_c
		0.0	;	# 198 tRNA_c_charging_M_trp_L_c_gfp::M_trp_L_c+M_atp_c+tRNA_c+M_h2o_c --> M_trp_L_c_tRNA_c+M_amp_c+M_ppi_c
		0.0	;	# 199 tRNA_c_charging_M_tyr_L_c_gfp::11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA_c+11.0*M_h2o_c --> 11.0*M_tyr_L_c_tRNA_c+11.0*M_amp_c+11.0*M_ppi_c
		0.0	;	# 200 tRNA_c_charging_M_val_L_c_gfp::18.0*M_val_L_c+18.0*M_atp_c+18.0*tRNA_c+18.0*M_h2o_c --> 18.0*M_val_L_c_tRNA_c+18.0*M_amp_c+18.0*M_ppi_c
		0.0	;	# 201 tRNA_c_c_exchange::[] --> tRNA_c
		0.0	;	# 202 tRNA_c_c_exchange_reverse::tRNA_c --> []
		0.0	;	# 203 PROTEIN_export_GFP::PROTEIN_gfp --> gfp_e
		0.0	;	# 204 M_h2s_c_exchange::[] --> M_h2s_c
		0.0	;	# 205 M_h2s_c_exchange_reverse::M_h2s_c --> []
		0.0	;	# 206 M_h2o_c_exchange::[] --> M_h2o_c
		0.0	;	# 207 M_h2o_c_exchange_reverse::M_h2o_c --> []
		0.0	;	# 208 M_pi_c_exchange::[] --> M_pi_c
		0.0	;	# 209 M_pi_c_exchange_reverse::M_pi_c --> []
		0.0	;	# 210 M_nh3_c_exchange::[] --> M_nh3_c
		0.0	;	# 211 M_nh3_c_exchange_reverse::M_nh3_c --> []
		0.0	;	# 212 M_hco3_c_exchange::[] --> M_hco3_c
		0.0	;	# 213 M_etoh_c_exchange::M_etoh_c --> []
		0.0	;	# 214 M_mglx_c_exchange::M_mglx_c --> []
		0.0	;	# 215 M_prop_c_exchange::M_prop_c --> []
		0.0	;	# 216 M_indole_c_exchange::M_indole_c --> []
		0.0	;	# 217 M_cadav_c_exchange::M_cadav_c --> []
		0.0	;	# 218 M_gaba_c_exchange::M_gaba_c --> []
		0.0	;	# 219 M_glycoA_c_exchange::M_glycoA_c --> []
		0.0	;	# 220 M_o2_c_exchange::M_o2_e --> M_o2_c
		0.0	;	# 221 M_co2_c_exchange::M_co2_c --> M_co2_e
		0.0	;	# 222 M_h_c_exchange::M_h_c --> M_h_e
		0.0	;	# 223 M_h_c_exchange_reverse::M_h_e --> M_h_c
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"R_gp_1::M_maltodextrin6_c+M_pi_c --> M_maltodextrin5_c+M_g1p_c"
		"R_gp_2::M_maltodextrin5_c+M_pi_c --> M_maltodextrin4_c+M_g1p_c"
		"R_gp_3::M_maltodextrin4_c+M_pi_c --> M_maltodextrin3_c+M_g1p_c"
		"R_gp_4::M_maltodextrin3_c+M_pi_c --> M_maltose_c+M_g1p_c"
		"R_gp::M_maltose_c+M_pi_c --> M_glc_D_c+M_g1p_c"
		"R_pgm::M_g1p_c --> M_g6p_c"
		"R_glk_atp::M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c"
		"R_pgi::M_g6p_c --> M_f6p_c"
		"R_pgi_reverse::M_f6p_c --> M_g6p_c"
		"R_pfk::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c"
		"R_fdp::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c"
		"R_fbaA::M_fdp_c --> M_dhap_c+M_g3p_c"
		"R_fbaA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c"
		"R_tpiA::M_dhap_c --> M_g3p_c"
		"R_tpiA_reverse::M_g3p_c --> M_dhap_c"
		"R_gpsA::M_dhap_c+M_nadph_c+M_h_c --> M_gly3p_c+M_nadp_c"
		"R_gpsA_reverse::M_gly3p_c+M_nadp_c --> M_dhap_c+M_nadph_c+M_h_c"
		"R_gapA::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c"
		"R_gapA_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c"
		"R_pgk::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c"
		"R_pgk_reverse::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c"
		"R_gpm::M_3pg_c --> M_2pg_c"
		"R_gpm_reverse::M_2pg_c --> M_3pg_c"
		"R_eno::M_2pg_c --> M_h2o_c+M_pep_c"
		"R_eno_reverse::M_h2o_c+M_pep_c --> M_2pg_c"
		"R_pyk::M_adp_c+M_pep_c --> M_atp_c+M_pyr_c"
		"R_pck::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c"
		"R_ppc::M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c"
		"R_pdh::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c"
		"R_pps::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c"
		"R_zwf::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c"
		"R_zwf_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c"
		"R_pgl::M_6pgl_c+M_h2o_c --> M_6pgc_c"
		"R_gnd::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c"
		"R_rpe::M_ru5p_D_c --> M_xu5p_D_c"
		"R_rpe_reverse::M_xu5p_D_c --> M_ru5p_D_c"
		"R_rpi::M_r5p_c --> M_ru5p_D_c"
		"R_rpi_reverse::M_ru5p_D_c --> M_r5p_c"
		"R_talAB::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c"
		"R_talAB_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c"
		"R_tkt1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c"
		"R_tkt1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c"
		"R_tkt2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c"
		"R_tkt2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c"
		"R_edd::M_6pgc_c --> M_2ddg6p_c+M_h2o_c"
		"R_eda::M_2ddg6p_c --> M_g3p_c+M_pyr_c"
		"R_gltA::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c"
		"R_acn::M_cit_c --> M_icit_c"
		"R_acn_reverse::M_icit_c --> M_cit_c"
		"R_icd::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c"
		"R_icd_reverse::M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c"
		"R_sucAB::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c"
		"R_sucCD::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c"
		"R_sdh::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c"
		"R_frd::M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c"
		"R_fum::M_fum_c+M_h2o_c --> M_mal_L_c"
		"R_fum_reverse::M_mal_L_c --> M_fum_c+M_h2o_c"
		"R_mdh::M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c"
		"R_mdh_reverse::M_oaa_c+M_h_c+M_nadh_c --> M_mal_L_c+M_nad_c"
		"R_cyd::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_hx_c"
		"R_cyo::4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_hx_c"
		"R_app::2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_hx_c"
		"R_atp::M_adp_c+M_pi_c+4.0*M_hx_c --> M_atp_c+4.0*M_h_c+M_h2o_c"
		"R_DNP_leak::M_hx_c --> M_h_c"
		"R_nuo::3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_hx_c"
		"R_pnt1::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c"
		"R_pnt2::M_nadh_c+M_nadp_c+2.0*M_hx_c --> 2.0*M_h_c+M_nad_c+M_nadph_c"
		"R_ndh1::M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c"
		"R_ndh2::M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c"
		"R_ppa::M_ppi_c+M_h2o_c --> 2.0*M_pi_c"
		"R_aceA::M_icit_c --> M_glx_c+M_succ_c"
		"R_aceB::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c"
		"R_maeA::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c"
		"R_maeB::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c"
		"R_pta::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c"
		"R_pta_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c"
		"R_ackA::M_actp_c+M_adp_c --> M_ac_c+M_atp_c"
		"R_ackA_reverse::M_ac_c+M_atp_c --> M_actp_c+M_adp_c"
		"R_acs::M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c"
		"R_adhE::M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c"
		"R_adhE_reverse::M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c"
		"R_ldh::M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c"
		"R_ldh_reverse::M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c"
		"R_pflAB::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c"
		"R_alaAC::M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c"
		"R_alaAC_reverse::M_ala_L_c+M_akg_c --> M_pyr_c+M_glu_L_c"
		"R_arg::M_accoa_c+2.0*M_glu_L_c+3.0*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c --> M_coa_c+2.0*M_adp_c+2.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c"
		"R_aspC::M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c"
		"R_asnB::M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c"
		"R_asnA::M_asp_L_c+M_atp_c+M_nh3_c --> M_asn_L_c+M_ppi_c+M_amp_c"
		"R_cysEMK::M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_ac_c"
		"R_gltBD::M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c"
		"R_gdhA::M_akg_c+M_nadph_c+M_nh3_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c"
		"R_gdhA_reverse::M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh3_c+M_h_c"
		"R_glnA::M_glu_L_c+M_atp_c+M_nh3_c --> M_gln_L_c+M_adp_c+M_pi_c"
		"R_glyA::M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c"
		"R_his::M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c"
		"R_ile::M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c"
		"R_leu::2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c"
		"R_lys::M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c"
		"R_met::M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c"
		"R_phe::M_chor_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c"
		"R_pro::M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c"
		"R_serABC::M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c"
		"R_thr::M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c"
		"R_trp::M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c"
		"R_tyr::M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c"
		"R_val::2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c"
		"R_arg_deg::M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c --> M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c"
		"R_asp_deg::M_asp_L_c --> M_fum_c+M_nh3_c"
		"R_asn_deg::M_asn_L_c+M_amp_c+M_ppi_c --> M_nh3_c+M_asp_L_c+M_atp_c"
		"R_gly_deg::M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c"
		"R_mglx_deg::M_mglx_c+M_nad_c+M_h2o_c --> M_pyr_c+M_nadh_c+M_h_c"
		"R_ser_deg::M_ser_L_c --> M_nh3_c+M_pyr_c"
		"R_pro_deg::M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c"
		"R_thr_deg1::M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c"
		"R_thr_deg2::M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c"
		"R_thr_deg3::M_thr_L_c+M_pi_c+M_adp_c --> M_nh3_c+M_for_c+M_atp_c+M_prop_c"
		"R_trp_deg::M_trp_L_c+M_h2o_c --> M_indole_c+M_nh3_c+M_pyr_c"
		"R_cys_deg::M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh3_c+M_pyr_c"
		"R_lys_deg::M_lys_L_c --> M_co2_c+M_cadav_c"
		"R_gln_deg::M_gln_L_c+M_h2o_c --> M_nh3_c+M_glu_L_c"
		"R_glu_deg::M_glu_L_c --> M_co2_c+M_gaba_c"
		"R_gaba_deg1::M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadh_c"
		"R_gaba_deg2::M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadph_c"
		"R_chor::M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c"
		"R_fol_e::M_gtp_c+4.0*M_h2o_c --> M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c"
		"R_fol_1::M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c"
		"R_fol_2a::M_4adochor_c --> M_4abz_c+M_pyr_c"
		"R_fol_2b::M_4abz_c+M_78mdp_c --> M_78dhf_c+M_h2o_c"
		"R_fol_3::M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c"
		"R_fol_4::M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c"
		"R_gly_fol::M_gly_L_c+M_thf_c+M_nad_c --> M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c"
		"R_gly_fol_reverse::M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c --> M_gly_L_c+M_thf_c+M_nad_c"
		"R_mthfd::M_mlthf_c+M_nadp_c --> M_methf_c+M_nadph_c"
		"R_mthfd_reverse::M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c"
		"R_mthfc::M_h2o_c+M_methf_c --> M_10fthf_c+M_h_c"
		"R_mthfc_reverse::M_10fthf_c+M_h_c --> M_h2o_c+M_methf_c"
		"R_mthfr2a::M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c"
		"R_mthfr2b::M_mlthf_c+M_h_c+M_nadph_c --> M_5mthf_c+M_nadp_c"
		"R_prpp_syn::M_r5p_c+M_atp_c --> M_prpp_c+M_amp_c"
		"R_or_syn_1::2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c+2.0*M_h_c"
		"R_or_syn_2::M_clasp_c+M_asp_L_c+M_q8_c --> M_or_c+M_q8h2_c+M_h2o_c+M_pi_c"
		"R_omp_syn::M_prpp_c+M_or_c --> M_omp_c+M_ppi_c"
		"R_ump_syn::M_omp_c --> M_ump_c+M_co2_c"
		"R_ctp_1::M_utp_c+M_atp_c+M_nh3_c --> M_ctp_c+M_adp_c+M_pi_c"
		"R_ctp_2::M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c"
		"R_A_syn_1::M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c"
		"R_A_syn_2::M_atp_c+M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c"
		"R_A_syn_3::M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c"
		"R_A_syn_4::M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c"
		"R_A_syn_5::M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c"
		"R_A_syn_6::M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c"
		"R_A_syn_7::M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c"
		"R_A_syn_8::M_saicar_c --> M_fum_c+M_aicar_c"
		"R_A_syn_9::M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c"
		"R_A_syn_10::M_faicar_c --> M_imp_c+M_h2o_c"
		"R_A_syn_12::M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c"
		"R_xmp_syn::M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c"
		"R_gmp_syn::M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c"
		"R_atp_amp::M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c"
		"R_utp_ump::M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c"
		"R_ctp_cmp::M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c"
		"R_gtp_gmp::M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c"
		"R_atp_adp::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c"
		"R_utp_adp::M_utp_c+M_h2o_c --> M_udp_c+M_pi_c"
		"R_ctp_adp::M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c"
		"R_gtp_adp::M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c"
		"R_udp_utp::M_udp_c+M_atp_c --> M_utp_c+M_adp_c"
		"R_cdp_ctp::M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c"
		"R_gdp_gtp::M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c"
		"R_atp_ump::M_atp_c+M_ump_c --> M_adp_c+M_udp_c"
		"R_atp_cmp::M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c"
		"R_atp_gmp::M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c"
		"R_adk_atp::M_amp_c+M_atp_c --> 2.0*M_adp_c"
		"transcriptional_initiation_gfp::GENE_gfp+RNAP --> OPEN_GENE_gfp"
		"transcription_gfp::OPEN_GENE_gfp+195.0*M_gtp_c+157.0*M_ctp_c+157.0*M_utp_c+208.0*M_atp_c+717.0*M_h2o_c --> mRNA_gfp+GENE_gfp+RNAP+717.0*M_ppi_c"
		"mRNA_degradation_gfp::mRNA_gfp --> 195.0*M_gmp_c+157.0*M_cmp_c+157.0*M_ump_c+208.0*M_amp_c"
		"translation_initiation_gfp::mRNA_gfp+RIBOSOME --> RIBOSOME_START_gfp"
		"translation_gfp::RIBOSOME_START_gfp+478.0*M_gtp_c+478.0*M_h2o_c+8.0*M_ala_L_c_tRNA_c+6.0*M_arg_L_c_tRNA_c+13.0*M_asn_L_c_tRNA_c+18.0*M_asp_L_c_tRNA_c+2.0*M_cys_L_c_tRNA_c+16.0*M_glu_L_c_tRNA_c+8.0*M_gln_L_c_tRNA_c+22.0*M_gly_L_c_tRNA_c+9.0*M_his_L_c_tRNA_c+12.0*M_ile_L_c_tRNA_c+21.0*M_leu_L_c_tRNA_c+20.0*M_lys_L_c_tRNA_c+6.0*M_met_L_c_tRNA_c+12.0*M_phe_L_c_tRNA_c+10.0*M_pro_L_c_tRNA_c+10.0*M_ser_L_c_tRNA_c+16.0*M_thr_L_c_tRNA_c+M_trp_L_c_tRNA_c+11.0*M_tyr_L_c_tRNA_c+18.0*M_val_L_c_tRNA_c --> RIBOSOME+mRNA_gfp+PROTEIN_gfp+478.0*M_gdp_c+478.0*M_pi_c+239.0*tRNA_c"
		"tRNA_c_charging_M_ala_L_c_gfp::8.0*M_ala_L_c+8.0*M_atp_c+8.0*tRNA_c+8.0*M_h2o_c --> 8.0*M_ala_L_c_tRNA_c+8.0*M_amp_c+8.0*M_ppi_c"
		"tRNA_c_charging_M_arg_L_c_gfp::6.0*M_arg_L_c+6.0*M_atp_c+6.0*tRNA_c+6.0*M_h2o_c --> 6.0*M_arg_L_c_tRNA_c+6.0*M_amp_c+6.0*M_ppi_c"
		"tRNA_c_charging_M_asn_L_c_gfp::13.0*M_asn_L_c+13.0*M_atp_c+13.0*tRNA_c+13.0*M_h2o_c --> 13.0*M_asn_L_c_tRNA_c+13.0*M_amp_c+13.0*M_ppi_c"
		"tRNA_c_charging_M_asp_L_c_gfp::18.0*M_asp_L_c+18.0*M_atp_c+18.0*tRNA_c+18.0*M_h2o_c --> 18.0*M_asp_L_c_tRNA_c+18.0*M_amp_c+18.0*M_ppi_c"
		"tRNA_c_charging_M_cys_L_c_gfp::2.0*M_cys_L_c+2.0*M_atp_c+2.0*tRNA_c+2.0*M_h2o_c --> 2.0*M_cys_L_c_tRNA_c+2.0*M_amp_c+2.0*M_ppi_c"
		"tRNA_c_charging_M_glu_L_c_gfp::16.0*M_glu_L_c+16.0*M_atp_c+16.0*tRNA_c+16.0*M_h2o_c --> 16.0*M_glu_L_c_tRNA_c+16.0*M_amp_c+16.0*M_ppi_c"
		"tRNA_c_charging_M_gln_L_c_gfp::8.0*M_gln_L_c+8.0*M_atp_c+8.0*tRNA_c+8.0*M_h2o_c --> 8.0*M_gln_L_c_tRNA_c+8.0*M_amp_c+8.0*M_ppi_c"
		"tRNA_c_charging_M_gly_L_c_gfp::22.0*M_gly_L_c+22.0*M_atp_c+22.0*tRNA_c+22.0*M_h2o_c --> 22.0*M_gly_L_c_tRNA_c+22.0*M_amp_c+22.0*M_ppi_c"
		"tRNA_c_charging_M_his_L_c_gfp::9.0*M_his_L_c+9.0*M_atp_c+9.0*tRNA_c+9.0*M_h2o_c --> 9.0*M_his_L_c_tRNA_c+9.0*M_amp_c+9.0*M_ppi_c"
		"tRNA_c_charging_M_ile_L_c_gfp::12.0*M_ile_L_c+12.0*M_atp_c+12.0*tRNA_c+12.0*M_h2o_c --> 12.0*M_ile_L_c_tRNA_c+12.0*M_amp_c+12.0*M_ppi_c"
		"tRNA_c_charging_M_leu_L_c_gfp::21.0*M_leu_L_c+21.0*M_atp_c+21.0*tRNA_c+21.0*M_h2o_c --> 21.0*M_leu_L_c_tRNA_c+21.0*M_amp_c+21.0*M_ppi_c"
		"tRNA_c_charging_M_lys_L_c_gfp::20.0*M_lys_L_c+20.0*M_atp_c+20.0*tRNA_c+20.0*M_h2o_c --> 20.0*M_lys_L_c_tRNA_c+20.0*M_amp_c+20.0*M_ppi_c"
		"tRNA_c_charging_M_met_L_c_gfp::6.0*M_met_L_c+6.0*M_atp_c+6.0*tRNA_c+6.0*M_h2o_c --> 6.0*M_met_L_c_tRNA_c+6.0*M_amp_c+6.0*M_ppi_c"
		"tRNA_c_charging_M_phe_L_c_gfp::12.0*M_phe_L_c+12.0*M_atp_c+12.0*tRNA_c+12.0*M_h2o_c --> 12.0*M_phe_L_c_tRNA_c+12.0*M_amp_c+12.0*M_ppi_c"
		"tRNA_c_charging_M_pro_L_c_gfp::10.0*M_pro_L_c+10.0*M_atp_c+10.0*tRNA_c+10.0*M_h2o_c --> 10.0*M_pro_L_c_tRNA_c+10.0*M_amp_c+10.0*M_ppi_c"
		"tRNA_c_charging_M_ser_L_c_gfp::10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA_c+10.0*M_h2o_c --> 10.0*M_ser_L_c_tRNA_c+10.0*M_amp_c+10.0*M_ppi_c"
		"tRNA_c_charging_M_thr_L_c_gfp::16.0*M_thr_L_c+16.0*M_atp_c+16.0*tRNA_c+16.0*M_h2o_c --> 16.0*M_thr_L_c_tRNA_c+16.0*M_amp_c+16.0*M_ppi_c"
		"tRNA_c_charging_M_trp_L_c_gfp::M_trp_L_c+M_atp_c+tRNA_c+M_h2o_c --> M_trp_L_c_tRNA_c+M_amp_c+M_ppi_c"
		"tRNA_c_charging_M_tyr_L_c_gfp::11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA_c+11.0*M_h2o_c --> 11.0*M_tyr_L_c_tRNA_c+11.0*M_amp_c+11.0*M_ppi_c"
		"tRNA_c_charging_M_val_L_c_gfp::18.0*M_val_L_c+18.0*M_atp_c+18.0*tRNA_c+18.0*M_h2o_c --> 18.0*M_val_L_c_tRNA_c+18.0*M_amp_c+18.0*M_ppi_c"
		"tRNA_c_c_exchange::[] --> tRNA_c"
		"tRNA_c_c_exchange_reverse::tRNA_c --> []"
		"PROTEIN_export_GFP::PROTEIN_gfp --> gfp_e"
		"M_h2s_c_exchange::[] --> M_h2s_c"
		"M_h2s_c_exchange_reverse::M_h2s_c --> []"
		"M_h2o_c_exchange::[] --> M_h2o_c"
		"M_h2o_c_exchange_reverse::M_h2o_c --> []"
		"M_pi_c_exchange::[] --> M_pi_c"
		"M_pi_c_exchange_reverse::M_pi_c --> []"
		"M_nh3_c_exchange::[] --> M_nh3_c"
		"M_nh3_c_exchange_reverse::M_nh3_c --> []"
		"M_hco3_c_exchange::[] --> M_hco3_c"
		"M_etoh_c_exchange::M_etoh_c --> []"
		"M_mglx_c_exchange::M_mglx_c --> []"
		"M_prop_c_exchange::M_prop_c --> []"
		"M_indole_c_exchange::M_indole_c --> []"
		"M_cadav_c_exchange::M_cadav_c --> []"
		"M_gaba_c_exchange::M_gaba_c --> []"
		"M_glycoA_c_exchange::M_glycoA_c --> []"
		"M_o2_c_exchange::M_o2_e --> M_o2_c"
		"M_co2_c_exchange::M_co2_c --> M_co2_e"
		"M_h_c_exchange::M_h_c --> M_h_e"
		"M_h_c_exchange_reverse::M_h_e --> M_h_c"
	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"GENE_gfp"
		"M_10fthf_c"
		"M_13dpg_c"
		"M_2ddg6p_c"
		"M_2pg_c"
		"M_3pg_c"
		"M_4abz_c"
		"M_4adochor_c"
		"M_5mthf_c"
		"M_5pbdra"
		"M_6pgc_c"
		"M_6pgl_c"
		"M_78dhf_c"
		"M_78mdp_c"
		"M_ac_c"
		"M_accoa_c"
		"M_actp_c"
		"M_adp_c"
		"M_aicar_c"
		"M_air_c"
		"M_akg_c"
		"M_ala_L_c"
		"M_ala_L_c_tRNA_c"
		"M_amp_c"
		"M_arg_L_c"
		"M_arg_L_c_tRNA_c"
		"M_asn_L_c"
		"M_asn_L_c_tRNA_c"
		"M_asp_L_c"
		"M_asp_L_c_tRNA_c"
		"M_atp_c"
		"M_cadav_c"
		"M_cair_c"
		"M_cdp_c"
		"M_chor_c"
		"M_cit_c"
		"M_clasp_c"
		"M_cmp_c"
		"M_co2_c"
		"M_coa_c"
		"M_ctp_c"
		"M_cys_L_c"
		"M_cys_L_c_tRNA_c"
		"M_dhap_c"
		"M_dhf_c"
		"M_e4p_c"
		"M_etoh_c"
		"M_f6p_c"
		"M_faicar_c"
		"M_fdp_c"
		"M_fgam_c"
		"M_fgar_c"
		"M_for_c"
		"M_fum_c"
		"M_g1p_c"
		"M_g3p_c"
		"M_g6p_c"
		"M_gaba_c"
		"M_gar_c"
		"M_gdp_c"
		"M_glc_D_c"
		"M_gln_L_c"
		"M_gln_L_c_tRNA_c"
		"M_glu_L_c"
		"M_glu_L_c_tRNA_c"
		"M_glx_c"
		"M_gly3p_c"
		"M_gly_L_c"
		"M_gly_L_c_tRNA_c"
		"M_glycoA_c"
		"M_gmp_c"
		"M_gtp_c"
		"M_h2o2_c"
		"M_h2o_c"
		"M_h2s_c"
		"M_h_c"
		"M_hco3_c"
		"M_his_L_c"
		"M_his_L_c_tRNA_c"
		"M_hx_c"
		"M_icit_c"
		"M_ile_L_c"
		"M_ile_L_c_tRNA_c"
		"M_imp_c"
		"M_indole_c"
		"M_lac_D_c"
		"M_leu_L_c"
		"M_leu_L_c_tRNA_c"
		"M_lys_L_c"
		"M_lys_L_c_tRNA_c"
		"M_mal_L_c"
		"M_maltodextrin3_c"
		"M_maltodextrin4_c"
		"M_maltodextrin5_c"
		"M_maltodextrin6_c"
		"M_maltose_c"
		"M_met_L_c"
		"M_met_L_c_tRNA_c"
		"M_methf_c"
		"M_mglx_c"
		"M_mlthf_c"
		"M_mql8_c"
		"M_mqn8_c"
		"M_nad_c"
		"M_nadh_c"
		"M_nadp_c"
		"M_nadph_c"
		"M_nh3_c"
		"M_o2_c"
		"M_oaa_c"
		"M_omp_c"
		"M_or_c"
		"M_pep_c"
		"M_phe_L_c"
		"M_phe_L_c_tRNA_c"
		"M_pi_c"
		"M_ppi_c"
		"M_pro_L_c"
		"M_pro_L_c_tRNA_c"
		"M_prop_c"
		"M_prpp_c"
		"M_pyr_c"
		"M_q8_c"
		"M_q8h2_c"
		"M_r5p_c"
		"M_ru5p_D_c"
		"M_s7p_c"
		"M_saicar_c"
		"M_ser_L_c"
		"M_ser_L_c_tRNA_c"
		"M_succ_c"
		"M_succoa_c"
		"M_thf_c"
		"M_thr_L_c"
		"M_thr_L_c_tRNA_c"
		"M_trp_L_c"
		"M_trp_L_c_tRNA_c"
		"M_tyr_L_c"
		"M_tyr_L_c_tRNA_c"
		"M_udp_c"
		"M_ump_c"
		"M_utp_c"
		"M_val_L_c"
		"M_val_L_c_tRNA_c"
		"M_xmp_c"
		"M_xu5p_D_c"
		"OPEN_GENE_gfp"
		"PROTEIN_gfp"
		"RIBOSOME"
		"RIBOSOME_START_gfp"
		"RNAP"
		"mRNA_gfp"
		"tRNA_c"
		"M_co2_e"
		"M_h_e"
		"M_o2_e"
		"gfp_e"
	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
