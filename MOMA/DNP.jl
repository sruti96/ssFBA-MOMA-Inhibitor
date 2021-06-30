# ----------------------------------------------------------------------------------- #
# Copyright (c) 2020 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

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

# includes -
include("Include.jl")

"""
    solve_dynamic_problem()
"""
function solve_dynamic_problem(dt,tEND,case)

    #Setup Sampling on TXTL Parameters
    TXTL_parameters = zeros(4,1)
    TXTL_parameters[1] = 0.060 + rand(1)[1]*(0.075-0.060); #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
    TXTL_parameters[2] = 15.0 + rand(1)[1]*(10-5);  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
    TXTL_parameters[3] = 2.0 + rand(1)[1]*(2.3-2.0); #RIBOSOME_concentration #uM  <0.0023mM (ACS SynBio Garamella 2016)
    TXTL_parameters[4] = 1.0 + rand(1)[1]*(2-1.5); #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

	# TXTL_parameters[1] = .0675; #RNAP_concentration # uM # 60-75nM (ACS SynBio Garamella 2016)
    # TXTL_parameters[2] = 20; #max_transcription_rate # >5 NT/s (ACS SynBio Garamella 2016)
    # TXTL_parameters[3] = 2.15; #RIBOSOME_concentration #uM #0.0016mM with 72% MaxActive (Underwood, Swartz, Puglisi 2005 Biotech Bioeng) & <0.0023mM (ACS SynBio Garamella 2016)
    # TXTL_parameters[4] = 1.5; #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

    TXTL_dictionary = generate_TXTL_dictionary(TXTL_parameters)
    data_dictionary = generate_default_data_dictionary()

    # optional: do we need to change the problem setup to maximize the objective function?
    data_dictionary["is_minimum_flag"] = update_is_minimum_flag(data_dictionary)
    data_dictionary["objective_coefficient_array"] = update_objective_coefficient_array(data_dictionary, case)

    number_of_species, number_of_fluxes = size(data_dictionary["stoichiometric_matrix"])

    data_dictionary["Enzyme_sample"] = 1.25 + rand(1)[1]*(0.5)
	# data_dictionary["Enzyme_sample"] = 1.25

	sample = 1.5 + rand(1)[1]*(0.5)
	# sample = 1.5
    data_dictionary["AA_pos"] = sample .*ones(20,1)
    data_dictionary["AA_neg"] = sample .*ones(20,1)
    data_dictionary["TCA_pos"] = sample .*ones(12,1)
    data_dictionary["TCA_neg"] = sample .*ones(12,1)
    data_dictionary["GLY_pos"] = sample .*ones(13,1)
    data_dictionary["GLY_neg"] = sample .*ones(13,1)
    data_dictionary["ENERGY_pos"] = sample .*ones(12,1)
    data_dictionary["ENERGY_neg"] = sample .*ones(12,1)
    data_dictionary["REDUCING_pos"] = sample .*ones(4,1)
    data_dictionary["REDUCING_neg"] = sample .*ones(4,1)

    t_vec = collect(0:dt:tEND*dt)

    # generate flux data constraints
    # using Mike's files here
    TXTL_dictionary = BuildFluxData(TXTL_dictionary,case,t_vec)
    species_array  = BuildSpeciesArray(TXTL_dictionary,tEND,number_of_species)

    # initialize all arrays
    objective_time = Any[]
    flux_time = zeros(tEND,number_of_fluxes)
    uptake_time = zeros(tEND,number_of_species)
	status_array = Any[]

    # solve the lp problem -
    for time_index = 1:tEND

        # update bounds
        data_dictionary = update_bounds_array(data_dictionary, TXTL_dictionary, species_array, time_index, case)

        # create a problem object, solve the problem and save solution
        problem_object = build_fba_problem_object(data_dictionary)
        problem_solution = solve_simulation_problem(problem_object)

        # update species, objective, flux, uptake and exit array
        species_array[time_index+1,:] = species_array[time_index,:] + problem_solution.uptake_array.*dt;
        push!(objective_time, problem_solution.objective_value)
        flux_time[time_index,:] = problem_solution.flux
        uptake_time[time_index,:] = problem_solution.uptake_array
		push!(status_array, problem_solution.status_flag)

    end

     return objective_time, flux_time, uptake_time, species_array, status_array
end


"""
    solve_moma_problem()
"""
function solve_moma_problem(dt,tEND,case)

	# solve control problem to get wild type flux distribution

	# initialize all arrays
    objective_time = Any[]
    flux_time = Any[]
    uptake_time = Any[]
	species_time = Any[]
	status_array = Any[]
	exit = false

	while exit==false
    	objective_time, flux_time, uptake_time, species_time, status_array = solve_dynamic_problem(dt,tEND,"control");
		if isempty(findall(status_array .!= 5))
			exit = true
		end
	end

    #Setup TXTL Parameters
	TXTL_parameters = zeros(4,1)
    TXTL_parameters[1] = 0.060 + rand(1)[1]*(0.075-0.060); #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
    TXTL_parameters[2] = 15.0 + rand(1)[1]*(10-5);  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
    TXTL_parameters[3] = 2.0 + rand(1)[1]*(2.3-2.0); #RIBOSOME_concentration #uM  <0.0023mM (ACS SynBio Garamella 2016)
    TXTL_parameters[4] = 1.0 + rand(1)[1]*(2-1.5); #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

    TXTL_dictionary = generate_TXTL_dictionary(TXTL_parameters)
    data_dictionary = generate_default_data_dictionary()
    number_of_species, number_of_fluxes = size(data_dictionary["stoichiometric_matrix"])

    data_dictionary["Enzyme_sample"] = 1.25 + rand(1)[1]*(0.5)
	sample = 1.0 + rand(1)[1]*(1)
    data_dictionary["AA_pos"] = sample .*ones(20,1)
    data_dictionary["AA_neg"] = sample .*ones(20,1)
    data_dictionary["TCA_pos"] = sample .*ones(12,1)
    data_dictionary["TCA_neg"] = sample .*ones(12,1)
    data_dictionary["GLY_pos"] = sample .*ones(13,1)
    data_dictionary["GLY_neg"] = sample .*ones(13,1)
    data_dictionary["ENERGY_pos"] = sample .*ones(12,1)
    data_dictionary["ENERGY_neg"] = sample .*[1,1,1,1,1,1,1,1,1,1,1,4]
    data_dictionary["REDUCING_pos"] = sample .*ones(4,1)
    data_dictionary["REDUCING_neg"] = sample .*ones(4,1)

    t_vec = collect(0:dt:tEND*dt)

    # Generate flux data constraints
    # using Mike's files here
    TXTL_dictionary = BuildFluxData(TXTL_dictionary,case,t_vec)
    species_array  = BuildSpeciesArray(TXTL_dictionary,tEND,number_of_species)

    # initialize arrays
    moma_objective_time = Any[]
    moma_flux_time = zeros(tEND,number_of_fluxes)
    moma_uptake_time = zeros(tEND,number_of_species)
    status_array = Any[]

    # solve the lp problem -
    for time_index = 1:tEND

        # update MOMA bounds
        data_dictionary = update_bounds_array(data_dictionary, TXTL_dictionary, species_array, time_index, case)

        # grab wild type solution from control fba solution for the specific time_index
        wild_type_flux_array = flux_time[time_index, :]

        # create a MOMA problem object (use the factory method in utilities ...)
        # solve the problem and save solution -
        MOMA_problem_object = build_moma_problem_object(data_dictionary, wild_type_flux_array)
        problem_solution = solve_simulation_problem(MOMA_problem_object)

        # update species, flux, uptake and exit arrays
        species_array[time_index+1,:] = species_array[time_index,:] + problem_solution.uptake_array.*dt;
        push!(moma_objective_time, problem_solution.objective_value)
        moma_flux_time[time_index,:] = problem_solution.flux
        moma_uptake_time[time_index,:] = problem_solution.uptake_array
        push!(status_array, problem_solution.status_flag)

    end

    return moma_objective_time, moma_flux_time, moma_uptake_time, species_array, flux_time, species_time, status_array
end

function run_ensemble(dt,tEND,case)

	i = 1
	while i <= no_samples
		moma_objective_time, moma_flux_time, moma_uptake_time, species_array, flux_time, species_time, status_array = solve_moma_problem(dt,tEND,case);

		#Save results if no errors in solution
		if isempty(findall(status_array .!= 1))
			writedlm("Ensemble/$(case)/Flux_$i.txt",moma_flux_time)
			writedlm("Ensemble/$(case)/Species_$i.txt",species_array)

			i += 1
		end
	end

end

# which case? dnp/tta
case = "dnp"
dt = 0.1 # hr
tEND = convert(Int64,16/dt)
no_samples = 100
t_sim = collect(0:dt:tEND*dt)

# execute solver -
@info("Start: solving the MOMA problem ...")
run_ensemble(dt,tEND,case)
@info("Stop: calculation completed.")

include("plot_ensemble.jl")
