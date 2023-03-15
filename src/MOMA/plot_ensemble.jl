using PyPlot
using CSV
using DataFrames
using Statistics
using DelimitedFiles

#===================Find mean and CI for Ensemble==================#
t_sim = collect(0.0:0.1:16)

case = "tta"
S1 = readdlm("Ensemble/$(case)/Species_1.txt")

no_time, no_species = size(S1)
no_samples = 100
Species_ensemble = zeros(no_time,no_species,no_samples)

for i = 1:no_samples
    Species_ensemble[:,:,i] = readdlm("Ensemble/$(case)/Species_$i.txt")
end

Species_mean = mean(Species_ensemble, dims = 3)
Species_err = std(Species_ensemble, dims = 3)
Species_pos = Species_mean + 1.96*Species_err
Species_neg = Species_mean - 1.96*Species_err

Species_idx = readdlm("config/SpeciesDict/Species_index.dat")[:,2]
t = [0;2;4;8;16]

# if case == "tta"
#     markercolor = "#FC766AFF"
#     shade = "#FC766AFF"
#     lcolor = "#FC766AFF"
#
# elseif case == "dnp"
#     markercolor = "#783937FF"
#     shade = "#783937FF"
#     lcolor = "#783937FF"
# end
# markercolor = "black"
# shade = "powderblue"
# lcolor = "black"

if case == "tta"
    markercolor = "dimgrey"
    shade = "lightgrey"
    lcolor = "dimgrey"
elseif case == "dnp"
    markercolor = "darkred"
    shade = "lightpink"
    lcolor = "darkred"
end

#0 = no
#1 = yes
plotline = 1
plotsave =  0
filedir = "Figures"

#--------------------mRNA & GFP----------------------------#
gfp = DataFrame(CSV.File("FinalDATA/Data/$(case)_mRNA_protein.dat"))
mean_gfp = gfp[1:5,:]
error_gfp = gfp[6:10,:]
gfp_sim_idx = [152;148] #Species index for mRNA and GFP

row = 1
col = 2

figure(1,figsize=(6.5,2.5))
for idx = 1:2
    subplot(row,col,idx)
    ylabel(names(gfp)[idx])
    errorbar(t, mean_gfp[:,idx], yerr=error_gfp[:,idx],fmt="o",markersize = 3,capsize=2,elinewidth=1,color=markercolor, label="Experimental Data")
    if idx == 1
        if plotline == 1
             plot(t_sim,Species_mean[:,gfp_sim_idx[idx]]*1e6,color=lcolor, label="TTA")
             axis([0,16.0,0,1000])
             xticks([0,4,8,12,16])
             yticks([0,300,600,900])
        end
        fill_between(t_sim,Species_pos[:,gfp_sim_idx[idx]]*1e6,Species_neg[:,gfp_sim_idx[idx]]*1e6,color=shade,alpha=0.5,linewidth=0, label="95% CI of Ensemble")
        legend(loc=1, frameon=false, fontsize=8)
    else
        if plotline == 1
            plot(t_sim,Species_mean[:,gfp_sim_idx[idx]]*1e3,color=lcolor)
            axis([0,16.0,0,30])
            xticks([0,4,8,12,16])
            yticks([0,10,20])
        end
        fill_between(t_sim,Species_pos[:,gfp_sim_idx[idx]]*1e3,Species_neg[:,gfp_sim_idx[idx]]*1e3,color=shade,alpha=0.5,linewidth=0)
    end
    xlabel("Time (h)")
end
tight_layout()

if plotsave == 1
    savefig("$(filedir)/$(case)_protein.pdf")
    close()
end

# # -------------------AminoAcids----------------------------#
# row = 5
# col = 4
#
# aa = CSV.read("FinalDATA/Data/$(case)_aa.dat")
# mean_aa = aa[1:5,:]
# error_aa = aa[6:10,:]
# aa_sim_idx = convert(Array{Int64},readdlm("config/SpeciesDict/aa_index.dat"))
# figure(2,figsize=(11,7))
# for idx = 1:20
#     subplot(row,col,idx)
#     ylabel(names(aa)[idx])
#     errorbar(t, mean_aa[:,idx], yerr=error_aa[:,idx],fmt="o",markersize=3,capsize=2,elinewidth=1,color=markercolor)
#     fill_between(t_sim,Species_pos[:,aa_sim_idx[idx]],Species_neg[:,aa_sim_idx[idx]],color=shade,linewidth=0,alpha = 0.2)
#     if plotline == 1
#         plot(t_sim,Species_mean[:,aa_sim_idx[idx]],color=lcolor)
#     end
#     axis([0,16.0,0,3])
#     xticks([0,4,8,12,16])
#     yticks([0,1,2,3])
#     if names(aa)[idx] == "ALA (mM)"
#         axis([0,16.0,0,10])
#         yticks([0,5,10])
#     elseif names(aa)[idx] == "GLU (mM)"
#         axis([0,16.0,0,150])
#         yticks([0,50,100,150])
#     elseif names(aa)[idx] == "VAL (mM)"
#         axis([0,16.0,0,10])
#         yticks([0,5,10])
#     elseif names(aa)[idx] == "SER (mM)"
#         axis([0,16.0,0,4])
#         yticks([0,2,4])
#     elseif names(aa)[idx] == "TRP (mM)"
#         axis([0,16.0,0,4])
#         yticks([0,2,4])
#     elseif names(aa)[idx] == "CYS (mM)"
#         axis([0,16.0,0,1])
#         yticks([0,0.5,1])
#     elseif names(aa)[idx] == "GLN (mM)"
#         axis([0,16.0,0,1])
#         yticks([0,0.5,1])
#     end
#     if idx > 16
#         xlabel("Time (h)")
#     end
# end
# tight_layout()
# if plotsave == 1
#     savefig("$(filedir)/$(case)_aa.pdf")
#     close()
# end
#
# #--------------------Energy Species---------------------------#
# row = 3
# col = 4
# Energy_index = readdlm("config/SpeciesDict/energy_idx.dat")
# energy = CSV.read("FinalDATA/Data/$(case)_energy.dat")
# mean_energy = energy[1:5,:]
# error_energy = energy[6:10,:]
# energy_sim_idx = Array{Int64}[]
# for i = 1:length(Energy_index)
#     push!(energy_sim_idx,findall(Energy_index[i] .== Species_idx))
# end
# figure(3,figsize=(11,5))
# for idx = 1:12
#     subplot(row,col,idx)
#     ylabel(names(energy)[idx])
#     errorbar(t, mean_energy[:,idx], yerr=error_energy[:,idx],fmt="o",markersize=3,capsize=2,elinewidth=1,color=markercolor)
#     fill_between(t_sim,Species_pos[:,energy_sim_idx[idx][1]],Species_neg[:,energy_sim_idx[idx][1]],color=shade,linewidth=0,alpha = 0.2)
#     if plotline == 1
#         plot(t_sim,Species_mean[:,energy_sim_idx[idx][1]],color=lcolor)
#     end
#     axis([0,16.0,0,2.0])
#     xticks([0,4,8,12,16])
#     if names(energy)[idx] == "ADP (mM)"
#         axis([0,16.0,0,4])
#         yticks([0,2,4])
#     end
#     if idx >= 9
#         xlabel("Time (h)")
#     end
# end
# tight_layout()
# if plotsave == 1
#     savefig("$(filedir)/$(case)_energy.pdf")
#     close()
# end
#
# #--------------------Glycolysis---------------------------#
# row = 4
# col = 3
#
# gly = CSV.read("FinalDATA/Data/$(case)_gly.dat")
# deletecols!(gly,"E4P (mM)")
# mean_gly = gly[1:5,:]
# error_gly = gly[6:10,:]
#
# gly_index = readdlm("config/SpeciesDict/gly_index_plot.dat")
# gly_sim_idx = Array{Int64}[]
# for i = 1:length(gly_index)
#     push!(gly_sim_idx,findall(gly_index[i] .== Species_idx))
# end
#
# figure(4, figsize=(9,6))
# for idx = 1:12
#     subplot(row,col,idx)
#     ylabel(names(gly)[idx])
#     errorbar(t, mean_gly[:,idx], yerr=error_gly[:,idx],fmt="o",markersize=3,capsize=2,elinewidth=1,color=markercolor)
#     fill_between(t_sim,Species_pos[:,gly_sim_idx[idx][1]],Species_neg[:,gly_sim_idx[idx][1]],color=shade,linewidth=0,alpha = 0.2)
#     if plotline == 1
#         plot(t_sim,Species_mean[:,gly_sim_idx[idx][1]],color=lcolor)
#     end
#     axis([0,16.0,0,0.1])
#     xticks([0,4,8,12,16])
#     #yticks([0,0.05,0.1])
#     if names(gly)[idx] == "3PG (mM)"
#         axis([0,16.0,0,8])
#         yticks([0,4,8])
#     elseif names(gly)[idx] == "F16P (mM)"
#         axis([0,16.0,0,6])
#         yticks([0,2,4,6])
#     elseif names(gly)[idx] == "GAP (mM)"
#         axis([0,16.0,0,7])
#         yticks([0,3,6])
#     elseif names(gly)[idx] == "Gly3P (mM)"
#         axis([0,16.0,0,2])
#         yticks([0,1,2])
#     elseif names(gly)[idx] == "E4P (mM)"
#         axis([0,16.0,0,2])
#         yticks([0,1,2])
#     elseif names(gly)[idx] == "RL5P (mM)"
#         axis([0,16.0,0,0.2])
#     elseif names(gly)[idx] == "R5P (mM)"
#         axis([0,16.0,0,0.2])
#     elseif names(gly)[idx] == "G6P (mM)"
#         axis([0,16.0,0,0.1])
#     elseif names(gly)[idx] == "Maltose (mM)"
#         axis([0,16.0,0,0.5])
#         #yticks([0,0.5,0.1])
#     end
#     if idx > 9
#         xlabel("Time (h)")
#     end
# end
# tight_layout()
# if plotsave == 1
#     savefig("$(filedir)/$(case)_gly.pdf")
#     close()
# end
#
# #--------------------Reducing Power---------------------------#
# row = 2
# col = 2
#
# reducing = CSV.read("FinalDATA/Data/$(case)_reducing.dat")
# deletecols!(reducing,"FAD (mM)")
# mean_reducing = reducing[1:5,:]
# error_reducing = reducing[6:10,:]
#
# reducing_index = readdlm("config/SpeciesDict/reducing_index.dat")
# reducing_sim_idx = Array{Int64}[]
# for i = 1:length(reducing_index)
#     push!(reducing_sim_idx,findall(reducing_index[i] .== Species_idx))
# end
#
# figure(5,figsize=(6,4))
# for idx = 1:4
#     subplot(row,col,idx)
#     ylabel(names(reducing)[idx])
#     errorbar(t, mean_reducing[:,idx], yerr=error_reducing[:,idx],fmt="o",markersize=3,capsize=2,elinewidth=1,color=markercolor)
#     fill_between(t_sim,Species_pos[:,reducing_sim_idx[idx][1]],Species_neg[:,reducing_sim_idx[idx][1]],color=shade,linewidth=0,alpha = 0.2)
#     if plotline == 1
#         plot(t_sim,Species_mean[:,reducing_sim_idx[idx][1]],color=lcolor)
#     end
#     axis([0,16.0,0,0.2])
#     xticks([0,4,8,12,16])
#     if names(reducing)[idx] == "NADH (mM)"
#         axis([0,16.0,0,2])
#     elseif names(reducing)[idx] == "NAD (mM)"
#         axis([0,16.0,0,0.5])
#     end
#     if idx > 2
#         xlabel("Time (h)")
#     end
# end
#
# tight_layout()
# if plotsave == 1
#     savefig("$(filedir)/$(case)_nad.pdf")
#     close()
# end
#
# #-------------------TCA---------------------------#
# row = 4
# col = 3
#
# tca = CSV.read("FinalDATA/Data/$(case)_tca.dat")
# mean_tca = tca[1:5,:]
# error_tca = tca[6:10,:]
#
# tca_index = readdlm("config/SpeciesDict/tca_index.dat")
# tca_sim_idx = Array{Int64}[]
# for i = 1:length(tca_index)
#     push!(tca_sim_idx,findall(tca_index[i] .== Species_idx))
# end
#
# figure(6,figsize=(9,6))
# for idx = 1:12
#     subplot(row,col,idx)
#     ylabel(names(tca)[idx])
#     errorbar(t, mean_tca[:,idx], yerr=error_tca[:,idx],fmt="o",markersize=3,capsize=2,elinewidth=1,color=markercolor)
#     fill_between(t_sim,Species_pos[:,tca_sim_idx[idx][1]],Species_neg[:,tca_sim_idx[idx][1]],color=shade,linewidth=0,alpha = 0.2)
#     if plotline == 1
#         plot(t_sim,Species_mean[:,tca_sim_idx[idx][1]],color=lcolor)
#     end
#     axis([0,16.0,0,1])
#     xticks([0,4,8,12,16])
#     yticks([0,0.5,1])
#     if names(tca)[idx] == "ACE (mM)"
#         axis([0,16.0,0,80])
#         yticks([0,20,40,60,80])
#     elseif names(tca)[idx] == "LAC (mM)"
#         axis([0,16.0,0,12])
#         yticks([0,4,8,12])
#     elseif names(tca)[idx] == "PYR (mM)"
#         axis([0,16.0,0,15])
#         yticks([0,5,10,15])
#     elseif names(tca)[idx] == "SUCC (mM)"
#         axis([0,16.0,0,20])
#         yticks([0,5,10,15,20])
#     elseif names(tca)[idx] == "CIT (mM)"
#         axis([0,16.0,0,1])
#         yticks([0,0.5,1])
#     elseif names(tca)[idx] == "aKG (mM)"
#         axis([0,16.0,0,2])
#         yticks([0,1,2])
#     elseif names(tca)[idx] == "OAA (mM)"
#         axis([0,16.0,0,0.2])
#         yticks([0,0.1,0.2])
#     elseif names(tca)[idx] == "FUM (mM)"
#         axis([0,16.0,0,3])
#         yticks([0,1.5,3])
#     elseif names(tca)[idx] == "MAL (mM)"
#         axis([0,16.0,0,6])
#         yticks([0,3,6])
#     elseif names(tca)[idx] == "PEP (mM)"
#         axis([0,16.0,0,6])
#         yticks([0,3,6])
#     end
#     if idx > 9
#         xlabel("Time (h)")
#     end
# end
# tight_layout()
# if plotsave == 1
#     savefig("$(filedir)/$(case)_tca.pdf")
#     close()
# end
