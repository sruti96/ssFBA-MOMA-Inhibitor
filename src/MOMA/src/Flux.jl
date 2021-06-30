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

function calculate_moma_flux_distribution(stoichiometric_matrix::Array{Float64,2}, flux_bounds_array::Array{Float64,2}, species_bounds_array::Array{Float64,2}, wild_type_flux_array::Array{Float64,1})

    # print(default_bounds_array)
    number_of_species, number_of_fluxes = size(stoichiometric_matrix)

    # Read in v initial
    # This is the wild type flux array
    vi = wild_type_flux_array

    # Set up variables
    vd = Variable(number_of_fluxes)

    # Setup problem
    problem = minimize(sumsquares(vd-vi))

    species_upper_bounds = species_bounds_array[:, 2]
    species_lower_bounds = species_bounds_array[:, 1]

    problem.constraints += stoichiometric_matrix*vd <= species_upper_bounds
    problem.constraints += stoichiometric_matrix*vd >= species_lower_bounds

    flux_upper_bounds = flux_bounds_array[:, 2]
    flux_lower_bounds = flux_bounds_array[:, 1]

    problem.constraints += vd <= flux_upper_bounds
    problem.constraints += vd >= flux_lower_bounds

    # Call optimizer -
    solve!(problem, Gurobi.Optimizer())

    status_flag = problem.status

    if string(status_flag) == "OPTIMAL"
        status_flag = 1
    else
        status_flag = 0
    end

    objective_value = problem.optval
    flux = vd.value[:,1]
    uptake_array = stoichiometric_matrix*flux

    # Formulate the return tuple -
    return (objective_value, flux, uptake_array, status_flag)

end

# function calculate_fba_flux_distribution(stoichiometric_matrix::Array{Float64,2}, flux_bounds_array::Array{Float64,2}, species_bounds_array::Array{Float64,2}, objective_coefficient_array::Array{Float64,1}; min_flag::Bool = true)
#
#     # print(default_bounds_array)
#     number_of_species, number_of_fluxes = size(stoichiometric_matrix)
#
#     model = Model(GLPK.Optimizer)
#     @variable(model, v[1:number_of_fluxes])
#     @objective(model, Max, dot(objective_coefficient_array,v))
#
#
#     flux_upper_bounds = flux_bounds_array[:, 2]
#     flux_lower_bounds = flux_bounds_array[:, 1]
#
#     for i = 1:number_of_fluxes
#         @constraint(model, flux_lower_bounds[i] <= v[i] <= flux_upper_bounds[i])
#     end
#
#     species_upper_bounds = species_bounds_array[:, 2]
#     species_lower_bounds = species_bounds_array[:, 1]
#
#     b = stoichiometric_matrix*v
#
#     for i = 1:number_of_species
#         @constraint(model, species_lower_bounds[i] <= b[i] <= species_upper_bounds[i])
#     end
#
#     optimize!(model)
#
#     status_flag = termination_status(model)
#     objective_value = JuMP.objective_value(model)
#     flux = JuMP.value.(v)
#     uptake_array = stoichiometric_matrix*flux
#
#     # Formulate the return tuple -
#     return (objective_value, flux, uptake_array, status_flag, status_flag)
#
# end

# function calculate_fba_flux_distribution(stoichiometric_matrix::Array{Float64,2}, flux_bounds_array::Array{Float64,2}, species_bounds_array::Array{Float64,2}, objective_coefficient_array::Array{Float64,1}; min_flag::Bool = true)
#
#     # print(default_bounds_array)
#     number_of_species, number_of_fluxes = size(stoichiometric_matrix)
#
#     # Set up variables
#     vd = Variable(number_of_fluxes)
#
#     problem = maximize(dot(objective_coefficient_array,vd))
#
#     species_upper_bounds = species_bounds_array[:, 2]
#     species_lower_bounds = species_bounds_array[:, 1]
#
#     problem.constraints += stoichiometric_matrix*vd <= species_upper_bounds
#     problem.constraints += stoichiometric_matrix*vd >= species_lower_bounds
#
#     flux_upper_bounds = flux_bounds_array[:, 2]
#     flux_lower_bounds = flux_bounds_array[:, 1]
#
#     problem.constraints += vd <= flux_upper_bounds
#     problem.constraints += vd >= flux_lower_bounds
#
#     # Call optimizer -
#     solve!(problem, ECOS.Optimizer())
#
#     status_flag = problem.status
#     objective_value = problem.optval
#     flux = vd.value[:,1]
#     uptake_array = stoichiometric_matrix*flux
#
#     # Formulate the return tuple -
#     return (objective_value, flux, uptake_array, status_flag, status_flag)
#
# end

"""
    calculate_optimal_flux_distribution(S,[Lv,Uv],[Lx,Ux],c; min_flag=true)

    Computes the optimal metabolic flux distribution given the constraints.

    Inputs:
    `S` - stoichiometric_matrix (M x R)
    `[Lv,Uv]` - (R x 2) array of the flux lower (L) and upper (U) bounds
    `[Lx,Ux]` - (M x 2) array of species lower (L) and upper (U) bounds (if at steady-state, L = U = 0)
    `c` - R x 1 vector holding indexes for objective vector

    Outputs:
    `objective_value` - value of the objective function at the optimum
    `calculated_flux_array` - R x 1 flux array at the optimum
    `dual_value_array` - R x 1 dual values
    `uptake_array` - M x 1 array of S*v
    `exit_flag` = 0 if optimal
    `status_flag` = 5 if optimal
"""

function calculate_fba_flux_distribution(stoichiometric_matrix::Array{Float64,2}, default_bounds_array::Array{Float64,2}, species_bounds_array::Array{Float64,2}, objective_coefficient_array::Array{Float64,1}; min_flag::Bool = true)

    min_flag = false

    # set some constants -
    TIME_RESTART_LIM = 60

    # Get the stoichiometric_matrix from data_dictionary -
    (number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

    # # Setup the GLPK problem -
    lp_problem = GLPK.glp_create_prob();
    GLPK.glp_set_prob_name(lp_problem, "sample");
    GLPK.glp_set_obj_name(lp_problem, "objective")

    # Set solver parameters
    solver_parameters = GLPK.glp_smcp();
    GLPK.glp_init_smcp(solver_parameters);
    solver_parameters.msg_lev = GLPK.GLP_MSG_OFF;

    # Are we doing min -or- max?
    if min_flag == true
    	GLPK.glp_set_obj_dir(lp_problem, GLPK.GLP_MIN);
    else
    	GLPK.glp_set_obj_dir(lp_problem, GLPK.GLP_MAX);
    end

    # Set the number of constraints and fluxes -
    GLPK.glp_add_rows(lp_problem, number_of_species);
    GLPK.glp_add_cols(lp_problem, number_of_fluxes);

    # Setup flux bounds, and objective function -
    (number_of_fluxes,number_of_bounds) = size(default_bounds_array)
    for flux_index = 1:number_of_fluxes

    	flux_lower_bound = default_bounds_array[flux_index,1]
    	flux_upper_bound = default_bounds_array[flux_index,2]

    	# Check bounds type ... default is DB -
    	if (flux_upper_bound == flux_lower_bound)
    		flux_constraint_type = GLPK.GLP_FX
    	else
    		flux_constraint_type = GLPK.GLP_DB
    	end

    	# flux symbol? (later use name - for now, fake it)
    	flux_symbol = "R_"*string(flux_index)

    	# Set the bounds in GLPK -
    	GLPK.glp_set_col_name(lp_problem, flux_index, flux_symbol);
    	GLPK.glp_set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);
    end

    # Setup objective function -
    for (flux_index,obj_coeff) in enumerate(objective_coefficient_array)

    	# Set the objective function value in GLPK -
    	GLPK.glp_set_obj_coef(lp_problem, flux_index, obj_coeff);
    end

    # Setup problem constraints for the metabolites -
    for species_index = 1:number_of_species

    	species_lower_bound = species_bounds_array[species_index,1]
    	species_upper_bound = species_bounds_array[species_index,2]

    	# default
    	species_constraint_type = GLPK.GLP_FX
    	if (species_lower_bound != species_upper_bound)
    		species_constraint_type = GLPK.GLP_DB
    	end

    	# set the symbol -
    	species_symbol = "x_"*string(species_index)

    	# Set the species bounds in GLPK -
    	GLPK.glp_set_row_name(lp_problem, species_index, species_symbol);
    	GLPK.glp_set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);

    end

    # Setup the stoichiometric array -
    counter = 1;
    row_index_array = zeros(Cint,number_of_species*number_of_fluxes);
    col_index_array = zeros(Cint,number_of_species*number_of_fluxes);
    species_index_vector = collect(1:number_of_species);
    flux_index_vector = collect(1:number_of_fluxes);
    flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
    for species_index in species_index_vector
    	for flux_index in flux_index_vector
    		row_index_array[counter] = species_index;
    		col_index_array[counter] = flux_index;
    		flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
    		counter+=1;
    	end
    end

    GLPK.glp_load_matrix(lp_problem, number_of_species*number_of_fluxes, GLPK.offset(row_index_array), GLPK.offset(col_index_array), GLPK.offset(flat_stoichiometric_array));

    # Call the solver -
    exit_flag = GLPK.glp_simplex(lp_problem, solver_parameters);

    # Get the objective function value -
    objective_value = GLPK.glp_get_obj_val(lp_problem);

    # Get the calculated flux values from GLPK -
    calculated_flux_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	calculated_flux_array[flux_index] = GLPK.glp_get_col_prim(lp_problem, flux_index);
    end

    # Get the dual values -
    dual_value_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	dual_value_array[flux_index] = GLPK.glp_get_col_dual(lp_problem, flux_index);
    end

    # is this solution optimal?
    status_flag = GLPK.glp_get_status(lp_problem)

    # Calculate the uptake array -
    uptake_array = stoichiometric_matrix*calculated_flux_array;

    # Formulate the return tuple -
    return (objective_value, calculated_flux_array, uptake_array, exit_flag, status_flag);
end
