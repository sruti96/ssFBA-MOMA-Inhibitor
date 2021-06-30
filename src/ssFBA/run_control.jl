include("Include.jl")
include("BuildFluxData.jl")
include("BuildSpeciesArray.jl")
include("Bounds_control.jl")
include("Bounds_DNP.jl")
include("Bounds_TTA.jl")
# load the data dictionary -
function run(dt,tEND,case)

    TXTL_parameters = zeros(4,1)
    TXTL_parameters[1] = .0675; #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
    TXTL_parameters[2] = 20;  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
    TXTL_parameters[3] = 2.15; #RIBOSOME_concentration #uM #0.0016mM with 72% MaxActive (Underwood, Swartz, Puglisi 2005 Biotech Bioeng) & <0.0023mM (ACS SynBio Garamella 2016)
    TXTL_parameters[4] = 1.5; #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

    data_dictionary = DataDictionary(0,0,0)
    TXTL_dictionary = TXTLDictionary(TXTL_parameters)

    data_dictionary["objective_coefficient_array"][1] = 1;
    data_dictionary["objective_coefficient_array"][176:180] .= 1;

    if case == "dnp"
      data_dictionary["objective_coefficient_array"][64] = 1;
    end

    number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
    number_of_species = length(data_dictionary["species_bounds_array"][:,1])

    #=============Create Ensemble Parameters===================#
    data_dictionary["Enzyme_sample"] =1.25
    sample = 1.5
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


    #=============Generate flux data constraints===================#
    t_vec = collect(0:dt:tEND*dt)
    TXTL_dictionary = BuildFluxData(TXTL_dictionary,case,t_vec)
    Species_array  = BuildSpeciesArray(TXTL_dictionary,tEND,number_of_species)

    #Initialize Arrays
    flux_time = zeros(tEND,number_of_fluxes)
    uptake_time = zeros(tEND,number_of_species)
    dual_time = zeros(tEND,number_of_fluxes)
    exit_array = Any[]
    # solve the lp problem -
    for time_index = 1:tEND
        #Setup Bounds
        if case == "control"
            # data_dictionary = Bounds_control(data_dictionary,TXTL_dictionary,Species_array,time_index);
            data_dictionary["flux_bounds_array"] = update_flux_bounds_array(data_dictionary, TXTL, Species_vector, idx)
            data_dictionary["species_bounds_array"] = update_species_bounds_array(data_dictionary, TXTL, Species_vector, idx)
        elseif case == "dnp"
            data_dictionary = Bounds_DNP(data_dictionary,TXTL_dictionary,Species_array,time_index);
        elseif case == "tta"
            data_dictionary = Bounds_TTA(data_dictionary,TXTL_dictionary,Species_array,time_index);
        end
        (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

        #update Species, Flux, and Exit flag
        Species_array[time_index+1,:] = Species_array[time_index,:] + uptake_array.*dt;
        flux_time[time_index,:] = flux_array
        uptake_time[time_index,:] = uptake_array
        dual_time[time_index,:] = dual_array
        push!(exit_array, exit_flag)
    end
    return Species_array, flux_time, uptake_time, exit_array, dual_time
end

#options: "control", "dnp", "tta"
case = "dnp"
dt = 0.1 #time step (hr)
tEND = convert(Int64,16/dt)
t_sim = collect(0:dt:tEND*dt)
Species_array, flux_time, uptake_time, exit_array, dual_time = run(dt,tEND,case)


using PyPlot
include("plot_sim.jl")
