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

# ----------------------------------------------------------------------------------- #
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-02-12T18:49:52.564
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control(TXTL_dictionary::Dict{String,Any}, vargs...)

	# initialize the control_array -
	u = zeros(1)

	# alias the model parameters -
	Promoter = TXTL_dictionary["Promoter"]
	K1 = Promoter["K1"]
	K2 = Promoter["K2"]
	n = Promoter["coop_parameter"]
	k_binding = Promoter["k_binding"]
	I = Promoter["inducer_concentration"]

	f_function = (I^n)/(k_binding^n+I^n)
	u = (K1+K2*f_function)/(1+K1+K2*f_function)

	# return -
	return u[1]
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2021-02-12T18:49:52.621
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control(TXTL_dictionary::Dict{String,Any},vargs...)

	# initialize -
	w = ones(1)

	# User: override the default translation control impl here -
	# TODO: override ...

	# return -
	return w[1]
end
# ----------------------------------------------------------------------------------- #
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-02-12T18:49:52.658
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_allosteric_control(state::Array{Float64,1}, data_dictionary::Dict{String,Any}, vargs...)

	# initialize the v-vector -
	v = ones(250)

	# User: override the default allosteric control impl here -
	# TODO: override ...

	# return -
	return v
end
