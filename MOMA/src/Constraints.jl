function update_is_minimum_flag(data_dictionary::Dict{String,Any}, vargs...)::Bool

    """
        update_is_minimum_flag(data_dictionary, vargs...)
        Updates where the problem is setup as either a min (is_minimum_flag = true) or max (is_minimum_flag = false)
    """

    return false
end


function update_bounds_array(data_dictionary::Dict{String,Any}, TXTL_dictionary, species_array, time_index, case)::Dict{String,Any}

    # For both MOMA and FBA problems
    # using Mike's files here
    if case == "control"
        data_dictionary = Bounds_control(data_dictionary,TXTL_dictionary,species_array,time_index);
    elseif case == "dnp"
        data_dictionary = Bounds_DNP(data_dictionary,TXTL_dictionary,species_array,time_index);
    elseif case == "tta"
        data_dictionary = Bounds_TTA(data_dictionary,TXTL_dictionary,species_array,time_index);
    end

    return data_dictionary
end


function update_objective_coefficient_array(data_dictionary::Dict{String,Any}, case, vargs...)::Array{Float64,1}

    # For FBA problem

    data_dictionary["objective_coefficient_array"][1] = 1;
    data_dictionary["objective_coefficient_array"][176:180] .= 1;

    # only for DNP case -
    if case == "dnp"
      data_dictionary["objective_coefficient_array"][64] = 1;
    end

    return data_dictionary["objective_coefficient_array"]
end
