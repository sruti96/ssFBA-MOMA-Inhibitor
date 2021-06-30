include("Include.jl")
include("BuildFluxData.jl")
include("BuildSpeciesArray.jl")
include("Bounds_control.jl")
include("Bounds_DNP.jl")
include("Bounds_TTA.jl")

function run_ensemble(dt,tEND,case,no_samples)
    run_i = 1
    while run_i <= no_samples

        #Setup Sampling on TXTL Parameters
        TXTL_parameters = zeros(4,1)
        TXTL_parameters[1] = 0.060 + rand(1)[1]*(0.075-0.060); #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
        TXTL_parameters[2] = 15.0 + rand(1)[1]*(10);  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
        TXTL_parameters[3] = 2.0 + rand(1)[1]*(2.3-2.0); #RIBOSOME_concentration #uM  <0.0023mM (ACS SynBio Garamella 2016)
        TXTL_parameters[4] = 1.0 + rand(1)[1]*(2-1); #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)

        #======================Set up objective function=======================#
        data_dictionary = DataDictionary(0,0,0)
        TXTL_dictionary = TXTLDictionary(TXTL_parameters)
        data_dictionary["objective_coefficient_array"][1] = 1;
        data_dictionary["objective_coefficient_array"][176:180] .= 1;

        if case == "dnp"
            data_dictionary["objective_coefficient_array"][64] = 1;
        end
        number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
        number_of_species = length(data_dictionary["species_bounds_array"][:,1])

        #==================Create Ensemble Parameters=========================#
        data_dictionary["Enzyme_sample"] = 1.0 + rand(1)[1]*(0.5)
        LowerBound = 1.0
        UpperBound = 2.0
        Diff = UpperBound - LowerBound
        data_dictionary["AA_pos"] = LowerBound .+ rand(20).*(Diff)
        data_dictionary["AA_neg"] = LowerBound .+ rand(20).*(Diff)
        data_dictionary["TCA_pos"] = LowerBound .+ rand(12).*(Diff)
        data_dictionary["TCA_neg"] = LowerBound .+ rand(12).*(Diff)
        data_dictionary["GLY_pos"] = LowerBound .+ rand(13).*(Diff)
        data_dictionary["GLY_neg"] = LowerBound .+ rand(13).*(Diff)
        data_dictionary["ENERGY_pos"] = LowerBound .+ rand(12).*(Diff)
        data_dictionary["ENERGY_neg"] = LowerBound .+ rand(12).*(Diff)
        data_dictionary["REDUCING_pos"] = LowerBound .+ rand(4).*(Diff)
        data_dictionary["REDUCING_neg"] = LowerBound .+ rand(4).*(Diff)

        #=============Generate flux data constraints===================#
        t_vec = collect(0:dt:tEND*dt)
        TXTL_dictionary = BuildFluxData(TXTL_dictionary,case,t_vec)
        Species_array  = BuildSpeciesArray(TXTL_dictionary,tEND,number_of_species)

        flux_time = zeros(tEND,number_of_fluxes)
        uptake_time = zeros(tEND,number_of_species)
        exit_array = Any[]
        #solve the lp problem -
        for time_index = 1:tEND
            #Setup Bounds
            if case == "control"
                data_dictionary = Bounds_control(data_dictionary,TXTL_dictionary,Species_array,time_index);
            elseif case == "dnp"
                data_dictionary = Bounds_DNP(data_dictionary,TXTL_dictionary,Species_array,time_index);
            elseif case == "tta"
                data_dictionary = Bounds_TTA(data_dictionary,TXTL_dictionary,Species_array,time_index);
            end

            #Run solver
            (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

            #Update Species, Flux, and Exit flag
            Species_array[time_index+1,:] = Species_array[time_index,:] + uptake_array.*dt;
            flux_time[time_index,:] = flux_array
            uptake_time[time_index,:] = uptake_array
            push!(exit_array, exit_flag)
        end

        #Save results if no errors in solution
        if isempty(findall(exit_array .!= 5))
            writedlm("Ensemble/$(case)/Flux_$run_i.txt",flux_time)
            writedlm("Ensemble/$(case)/Species_$run_i.txt",Species_array)
            run_i += 1
        end
    end
end


#options: "control", "dnp", "tta"
case = "tta"
dt = 0.1 #time step (hr)
tEND = convert(Int64,16/dt)
no_samples = 100

run_ensemble(dt,tEND,case,no_samples)


using PyPlot
include("plot_ensemble.jl")
