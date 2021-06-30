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
function solve_fba_problem(problem::VLFBAProblem, vargs...)::VLFBASolution

    # get stuff to solve the static fba problem -
    stoichiometric_matrix = problem.stoichiometric_matrix
    flux_bounds_array = problem.flux_bounds_array
    species_bounds_array = problem.species_bounds_array
    obj_coeff_vector = problem.obj_coeff_vector
    is_minimum_flag = problem.is_minimum_flag

    # solve the static fba problem -
    (objective_value, calculated_flux_array, uptake_array, exit_flag, status_flag) = calculate_fba_flux_distribution(stoichiometric_matrix,flux_bounds_array,
        species_bounds_array, obj_coeff_vector; min_flag=is_minimum_flag)

    # create a solution object and return to the user -
    return VLFBASolution(objective_value, calculated_flux_array, uptake_array, exit_flag, status_flag)
end

function solve_moma_problem(problem::VLMOMAProblem, vargs...)::VLMOMASolution

    # get stuff to solve the MOMA problem -
    stoichiometric_matrix = problem.stoichiometric_matrix
    flux_bounds_array = problem.flux_bounds_array
    species_bounds_array = problem.species_bounds_array
    wild_type_flux_array = problem.wild_type_flux_array

    # solve the MOMA problem -
    (objective_value, flux, uptake_array, status_flag) = calculate_moma_flux_distribution(stoichiometric_matrix, flux_bounds_array, species_bounds_array, wild_type_flux_array)

    # create a solution object and return to the user -
    return VLMOMASolution(objective_value, flux, uptake_array, status_flag)

end

function solve_simulation_problem(problem::VLAbstractProblem, vargs...)::VLAbstractSolution

    # switch on the type of problem -
    if typeof(problem) == VLFBAProblem
        return solve_fba_problem(problem, vargs...)
    elseif (typeof(problem) == VLMOMAProblem)
        return solve_moma_problem(problem, vargs...)
    end
end
