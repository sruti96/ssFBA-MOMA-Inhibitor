function build_biophysical_dictionary(path_to_biophysical_constants_file::String)::Dict{String,Any}

    # check the path - is it legit?
    check_file_existence(path_to_biophysical_constants_file)

    # TODO: we need to check if the biophysical_constants_file is a TOML file
    # TODO: is this a TOML file?

    # Load the model dictionary -
    model_dictionary = TOML.parsefile(path_to_biophysical_constants_file)

    # initialize -
    biophysical_dictionary = Dict{String,Any}()

    # convert -
    list_of_constant_dictionaries = model_dictionary["biophysical_constants"]
    for (key_string, local_dictionary) in list_of_constant_dictionaries

        # get the value -
        value = parse(Float64, local_dictionary["value"])

        # cache -
        biophysical_dictionary[key_string] = value
    end

    # return -
    return biophysical_dictionary
end

function build_ec_data_dictionary(path_to_ec_data_file::String)::Dict{String,Any}

    # check the path - is it legit?
    check_file_existence(path_to_ec_data_file)

    # TODO: we need to check if the path_to_ec_data_file is a TOML file
    # TODO: is this a TOML file?

    # load the ec model dictionary -
    ec_data_dictionary = TOML.parsefile(path_to_ec_data_file)

    # return -
    return ec_data_dictionary
end

function build_control_constants_dictionary(path_to_ec_data_file::String)::Dict{String,Any}

    # check the path - is it legit?
    check_file_existence(path_to_ec_data_file)

    # TODO: we need to check if the path_to_ec_data_file is a TOML file
    # TODO: is this a TOML file?

    # load the ec model dictionary -
    control_parameter_data_dictionary = TOML.parsefile(path_to_ec_data_file)

    # return -
    return control_parameter_data_dictionary
end

function build_fba_problem_object(data_dictionary::Dict{String,Any}; is_minimum_flag::Bool=true)::VLFBAProblem

    # takes stuff from the dd and bundles it up as a static problem object -

    # create a problem, and populate its fields -
    problem = VLFBAProblem()
    problem.stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
    problem.flux_bounds_array = data_dictionary["flux_bounds_array"]
    problem.species_bounds_array = data_dictionary["species_bounds_array"]
    problem.obj_coeff_vector = data_dictionary["objective_coefficient_array"]
    problem.is_minimum_flag = data_dictionary["is_minimum_flag"]
    return problem
end

function build_moma_problem_object(data_dictionary::Dict{String,Any}, wild_type_flux_array::Array{Float64,1})::VLMOMAProblem

    # takes stuff from the dd and bundles it up as a static problem object -
    problem = VLMOMAProblem()
    problem.stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
    problem.flux_bounds_array = data_dictionary["flux_bounds_array"]
    problem.species_bounds_array = data_dictionary["species_bounds_array"]
    problem.wild_type_flux_array = wild_type_flux_array
    return problem
end
