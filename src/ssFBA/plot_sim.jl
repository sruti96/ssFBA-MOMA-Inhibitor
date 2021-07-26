using PyPlot
using CSV
using DataFrames
Species_idx = readdlm("config/SpeciesDict/Species_index.dat")[:,2]
t = [0;2;4;8;15.9]

#--------------------mRNA & GFP----------------------------#
gfp = CSV.read("FinalDATA/Data/$(case)_mRNA_protein.dat")
mean_gfp = gfp[1:5,:]
error_gfp = gfp[6:10,:]
gfp_sim_idx = [152;148] #Species index for mRNA and GFP

row = 1
col = 2

figure(1,figsize=(8,4))
for idx = 1:2
    subplot(row,col,idx)
    ylabel(names(gfp)[idx])
    errorbar(t, mean_gfp[:,idx], yerr=error_gfp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="black")
    if idx == 1
        plot(t_sim,Species_array[:,gfp_sim_idx[idx]]*1e6,color="black")
        axis([0,16.0,0,1000])
        yticks([0,200,400,600,800,1000])
    else
        plot(t_sim,Species_array[:,gfp_sim_idx[idx]]*1e3,color="black")
        axis([0,16.0,0,30])
        xticks([0,4,8,12,16])
    end

    # if names(gfp)[idx] == ("mRNA (nM)")
    #
    # end
end
tight_layout()

# --------------------AminoAcids----------------------------#
row = 5
col = 4

aa = CSV.read("FinalDATA/Data/$(case)_aa.dat")
mean_aa = aa[1:5,:]
error_aa = aa[6:10,:]
aa_sim_idx = convert(Array{Int64},readdlm("config/SpeciesDict/aa_index.dat"))
figure(2,figsize=(10,8))
for idx = 1:20
    subplot(row,col,idx)
    ylabel(names(aa)[idx])
    errorbar(t, mean_aa[:,idx], yerr=error_aa[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="black")
    plot(t_sim,Species_array[:,aa_sim_idx[idx]],color="black")
    axis([0,16.0,0,3])
    xticks([0,4,8,12,16])
    yticks([0,1,2,3])
    if names(aa)[idx] == ("ALA (mM)")
        axis([0,16.0,0,10])
        yticks([0,5,10])
    elseif names(aa)[idx] == ("GLU (mM)")
        axis([0,16.0,0,150])
        yticks([0,50,100,150])
    elseif names(aa)[idx] == ("VAL (mM)")
        axis([0,16.0,0,10])
        yticks([0,5,10])
    elseif names(aa)[idx] == ("SER (mM)")
        axis([0,16.0,0,4])
        yticks([0,2,4])
    elseif names(aa)[idx] == ("TRP (mM)")
        axis([0,16.0,0,4])
        yticks([0,2,4])
    elseif names(aa)[idx] == ("CYS (mM)")
        axis([0,16.0,0,1])
        yticks([0,0.5,1])
    elseif names(aa)[idx] == ("GLN (mM)")
        axis([0,16.0,0,1])
        yticks([0,0.5,1])
    end
end
tight_layout()

#--------------------Energy Species---------------------------#
row = 3
col = 4
Energy_index = readdlm("config/SpeciesDict/energy_idx.dat")
energy = CSV.read("FinalDATA/Data/$(case)_energy.dat")
mean_energy = energy[1:5,:]
error_energy = energy[6:10,:]
energy_sim_idx = Array{Int64}[]
for i = 1:length(Energy_index)
    push!(energy_sim_idx,findall(Energy_index[i] .== Species_idx))
end
figure(3,figsize=(10,5))
for idx = 1:12
    subplot(row,col,idx)
    ylabel(names(energy)[idx])
    errorbar(t, mean_energy[:,idx], yerr=error_energy[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="black")
    plot(t_sim,Species_array[:,energy_sim_idx[idx]],color="black")
    axis([0,16.0,0,2.0])
    xticks([0,4,8,12,16])
    if names(energy)[idx] == ("ADP (mM)")
        axis([0,16.0,0,4])
        yticks([0,2,4])
    end
end
tight_layout()

#--------------------Glycolysis---------------------------#
row = 5
col = 3

gly = CSV.read("FinalDATA/Data/$(case)_gly.dat")
mean_gly = gly[1:5,:]
error_gly = gly[6:10,:]

gly_index = readdlm("config/SpeciesDict/gly_index.dat")
gly_sim_idx = Array{Int64}[]
for i = 1:length(gly_index)
    push!(gly_sim_idx,findall(gly_index[i] .== Species_idx))
end

figure(4, figsize=(8,6))
for idx = 1:13
    subplot(row,col,idx)
    ylabel(names(gly)[idx])
    errorbar(t, mean_gly[:,idx], yerr=error_gly[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="black")
    plot(t_sim,Species_array[:,gly_sim_idx[idx]],color="black")
    axis([0,16.0,0,0.1])
    xticks([0,4,8,12,16])
    #yticks([0,0.05,0.1])
    if names(gly)[idx] == ("3PG (mM)")
        axis([0,16.0,0,8])
        yticks([0,4,8])
    elseif names(gly)[idx] == ("F16P (mM)")
        axis([0,16.0,0,6])
        yticks([0,2,4,6])
    elseif names(gly)[idx] == ("GAP (mM)")
        axis([0,16.0,0,7])
        yticks([0,3,6])
    elseif names(gly)[idx] == ("Gly3P (mM)")
        axis([0,16.0,0,2])
        yticks([0,1,2])
    elseif names(gly)[idx] == ("E4P (mM)")
        axis([0,16.0,0,2])
        yticks([0,1,2])
    elseif names(gly)[idx] == ("RL5P (mM)")
        axis([0,16.0,0,0.2])
    elseif names(gly)[idx] == ("R5P (mM)")
        axis([0,16.0,0,0.2])
    elseif names(gly)[idx] == ("G6P (mM)")
        axis([0,16.0,0,0.1])
    elseif names(gly)[idx] == ("Maltose (mM)")
        axis([0,16.0,0,0.5])
        #yticks([0,0.5,0.1])
    end
end
tight_layout()

#--------------------Reducing Power---------------------------#
row = 2
col = 2

reducing = CSV.read("FinalDATA/Data/$(case)_reducing.dat")
deletecols!(reducing,("FAD (mM)"))
mean_reducing = reducing[1:5,:]
error_reducing = reducing[6:10,:]

reducing_index = readdlm("config/SpeciesDict/reducing_index.dat")
reducing_sim_idx = Array{Int64}[]
for i = 1:length(reducing_index)
    push!(reducing_sim_idx,findall(reducing_index[i] .== Species_idx))
end

figure(5,figsize=(6,3))
for idx = 1:4
    subplot(row,col,idx)
    ylabel(names(reducing)[idx])
    errorbar(t, mean_reducing[:,idx], yerr=error_reducing[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="black")
    plot(t_sim,Species_array[:,reducing_sim_idx[idx]],color="black")
    axis([0,16.0,0,0.2])
    xticks([0,4,8,12,16])
    if names(reducing)[idx] == ("NADH (mM)")
        axis([0,16.0,0,2])
    elseif names(reducing)[idx] == ("NAD (mM)")
        axis([0,16.0,0,0.5])
    end
end
tight_layout()

#-------------------TCA---------------------------#
row = 4
col = 3

tca = CSV.read("FinalDATA/Data/$(case)_tca.dat")
mean_tca = tca[1:5,:]
error_tca = tca[6:10,:]

tca_index = readdlm("config/SpeciesDict/tca_index.dat")
tca_sim_idx = Array{Int64}[]
for i = 1:length(tca_index)
    push!(tca_sim_idx,findall(tca_index[i] .== Species_idx))
end

figure(6,figsize=(8,6))
for idx = 1:12
    subplot(row,col,idx)
    ylabel(names(tca)[idx])
    errorbar(t, mean_tca[:,idx], yerr=error_tca[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="black")
    plot(t_sim,Species_array[:,tca_sim_idx[idx]],color="black")
    axis([0,16.0,0,1])
    xticks([0,4,8,12,16])
    yticks([0,0.5,1])
    if names(tca)[idx] == ("ACE (mM)")
        axis([0,16.0,0,60])
        yticks([0,20,40,60])
    elseif names(tca)[idx] == ("LAC (mM)")
        axis([0,16.0,0,10])
        yticks([0,5,10])
    elseif names(tca)[idx] == ("PYR (mM)")
        axis([0,16.0,0,15])
        yticks([0,5,10,15])
    elseif names(tca)[idx] == ("SUCC (mM)")
        axis([0,16.0,0,20])
        yticks([0,5,10,15,20])
    elseif names(tca)[idx] == ("CIT (mM)")
        axis([0,16.0,0,1])
        yticks([0,0.5,1])
    elseif names(tca)[idx] == ("aKG (mM)")
        axis([0,16.0,0,2])
        yticks([0,1,2])
    elseif names(tca)[idx] == ("OAA (mM)")
        axis([0,16.0,0,0.2])
        yticks([0,0.1,0.2])
    elseif names(tca)[idx] == ("FUM (mM)")
        axis([0,16.0,0,3])
        yticks([0,1.5,3])
    elseif names(tca)[idx] == ("MAL (mM)")
        axis([0,16.0,0,6])
        yticks([0,3,6])
    elseif names(tca)[idx] == ("PEP (mM)")
        axis([0,16.0,0,6])
        yticks([0,3,6])
    end
end
tight_layout()
