using DataFrames
using CSV
using DelimitedFiles
using SmoothingSplines
using PyCall
#np = pyimport("numpy")
@pyimport numpy as np

function BuildFluxData(TXTL,case,t_array)

    data_t = [0.0;2.0;4.0;8.0;16.0;]; #Time points of data
    data_met = CSV.read("config/SpeciesDict/$(case)_metabolites.dat")[1:5,:]
    data_aa = CSV.read("config/SpeciesDict/$(case)_aa.dat")[1:5,:]
    deletecols!(data_met,:FAD)

    no_mets = length(data_met[1,:])
    no_aa = length(data_aa[1,:])

    #===============Interpolate data for Smoothing function==========#
    t_vec = collect(0.0:0.01:16.0)
    #np.interp(time points where you want interpolation, given time vector, given data)
    sim_met = zeros(length(t_vec),no_mets)
    for i = 1:no_mets
        sim_met[:,i] = np.interp(t_vec,data_t,data_met[:,i])
    end

    sim_aa = zeros(length(t_vec),no_aa)
    for i = 1:no_aa
        sim_aa[:,i] = np.interp(t_vec,data_t,data_aa[:,i])
    end

    #===============Smooth concentration values==========#
    smooth_met = zeros(length(t_vec),no_mets)
    for i = 1:no_mets
        spl_met = fit(SmoothingSpline,t_vec,sim_met[:,i],100.0)
        y_met = SmoothingSplines.predict(spl_met)
        diff = y_met[1] - sim_met[1,i]
        smooth_met[:,i] = y_met .- diff
    end

    smooth_aa = zeros(length(t_vec),no_aa)
    for i = 1:no_aa
        spl_aa = fit(SmoothingSpline,t_vec,sim_aa[:,i],100.0)
        y_aa = SmoothingSplines.predict(spl_aa)
        diff = y_aa[1] - sim_aa[1,i]
        smooth_aa[:,i] = y_aa .- diff
    end

    #===============Interpolate Smooth Curves for Flux==========#
    #np.interp(time points where you want interpolation, given time vector, given data)
    sim_met_smooth = zeros(length(t_array),no_mets)
    for i = 1:no_mets
        sim_met_smooth[:,i] = np.interp(t_array,t_vec,smooth_met[:,i])
    end

    sim_aa_smooth = zeros(length(t_array),no_aa)
    for i = 1:no_aa
        sim_aa_smooth[:,i] = np.interp(t_array,t_vec,smooth_aa[:,i])
    end

    no_data_points = length(sim_met_smooth[:,1])
    #=============Extract Flux Values=========================#
    flux_met = zeros(no_data_points-1,no_mets)
    flux_aa = zeros(no_data_points-1,no_aa)
    for flux_idx = 1:(no_data_points-1)
        flux_met[flux_idx,:] .= (sim_met_smooth[flux_idx+1,:] - sim_met_smooth[flux_idx,:])/(t_array[flux_idx+1]-t_array[flux_idx])
        flux_aa[flux_idx,:] .= (sim_aa_smooth[flux_idx+1,:] - sim_aa_smooth[flux_idx,:])/(t_array[flux_idx+1]-t_array[flux_idx])
    end
    #================Find Species Index=======================#
    network_idx = CSV.read("config/Reactions.txt",header=[Symbol("idx"),Symbol("rxn_name"),Symbol("substrate"),Symbol("arrow"),Symbol("product")])

    AA_name = readdlm("config/SpeciesDict/AA_name_internal.txt")
    Met_name  = readdlm("config/SpeciesDict/Met_name_internal.txt")
    Species_idx = readdlm("config/SpeciesDict/Species_index.dat")

    Met_Idx = zeros(no_mets,1)
    MetIdx = convert(Array{Int64,2},Met_Idx)
    for i = 1:no_mets
        MetIdx[i,1] = findall(Met_name[i].==Species_idx[:,2])[1]
    end

    AA_Idx = zeros(no_aa,1)
    AAIdx = convert(Array{Int64,2},AA_Idx)
    for i = 1:no_aa
        AAIdx[i,1] = findall(AA_name[i].==Species_idx[:,2])[1]
    end

    #=====================Separate Metabolite Index into subsections==============================#
    #Energy
    energy = readdlm("config/SpeciesDict/energy_idx.dat")
    no_energy = length(energy)
    Energy_Idx = convert(Array{Int64,2},zeros(no_energy,1))
    for i = 1:no_energy
        Energy_Idx[i] = getindex(findall(energy[i].==Met_name))[1]
    end
    e_idx = sort(Energy_Idx, dims = 1)
    EnergyIdx = MetIdx[e_idx]
    flux_energy = flux_met[:,e_idx]

    #Glycolysis
    glycolysis = readdlm("config/SpeciesDict/gly_index.dat")
    no_gly = length(glycolysis)
    Gly_Idx = convert(Array{Int64,2},zeros(no_gly,1))
    for i = 1:no_gly
        Gly_Idx[i] = getindex(findall(glycolysis[i].==Met_name))[1]
    end
    g_idx = sort(Gly_Idx, dims = 1)
    GlyIdx = MetIdx[g_idx]
    flux_gly = flux_met[:,g_idx]

    #Reducing
    reducing = readdlm("config/SpeciesDict/reducing_index.dat")
    no_reducing = length(reducing)
    Reducing_Idx = convert(Array{Int64,2},zeros(no_reducing,1))
    for i = 1:no_reducing
        Reducing_Idx[i] = getindex(findall(reducing[i].==Met_name))[1]
    end
    r_idx = sort(Reducing_Idx, dims = 1)
    ReducingIdx = MetIdx[r_idx]
    flux_reducing = flux_met[:,r_idx]

    #TCA
    tca = readdlm("config/SpeciesDict/tca_index.dat")
    no_tca = length(tca)
    TCA_Idx = convert(Array{Int64,2},zeros(no_tca,1))
    for i = 1:no_tca
        TCA_Idx[i] = getindex(findall(tca[i].==Met_name))[1]
    end
    t_idx = sort(TCA_Idx, dims = 1)
    TCAIdx = MetIdx[t_idx]
    flux_tca = flux_met[:,t_idx]

    #===========================Load in Enzyme Activity Assays =======================================#
    rxn_enzyme_data = convert(Array{Int64,2},readdlm("config/EnzymeDict/rxn_index_enzyme_data.dat"))
    vmax_2h = readdlm("config/EnzymeDict/$(case)_enzyme_vmax.dat")[:,1]
    vmax_8h = readdlm("config/EnzymeDict/$(case)_enzyme_vmax.dat")[:,2]

    #==================Store In Dictionary========================#
    ConcData_dictionary  = Dict{Symbol,Any}()
    ConcData_dictionary[:Species_index_metabolite] = MetIdx
    ConcData_dictionary[:Species_index_aa] = AAIdx
    ConcData_dictionary[:Conc_met] = sim_met_smooth
    ConcData_dictionary[:Conc_aa] = sim_aa_smooth
    ConcData_dictionary[:Flux_met] = flux_met
    ConcData_dictionary[:Flux_aa] = flux_aa
    ConcData_dictionary[:Flux_energy] = flux_energy
    ConcData_dictionary[:Flux_gly] = flux_gly
    ConcData_dictionary[:Flux_reducing] = flux_reducing
    ConcData_dictionary[:Flux_tca] = flux_tca
    ConcData_dictionary[:index_energy] = EnergyIdx
    ConcData_dictionary[:index_gly] = GlyIdx
    ConcData_dictionary[:index_reducing] = ReducingIdx
    ConcData_dictionary[:index_tca] = TCAIdx
    ConcData_dictionary[:index_enzyme_activity] = rxn_enzyme_data
    ConcData_dictionary[:enzyme_vmax_2h] = vmax_2h
    ConcData_dictionary[:enzyme_vmax_8h] = vmax_8h

    TXTL[:ConcDictionary] = ConcData_dictionary
    return TXTL
end
