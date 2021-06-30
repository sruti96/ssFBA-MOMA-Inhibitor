
using PyPlot
using CSV
using DataFrames
using DelimitedFiles
t = [2;4;]
row = 5
col = 4

control = CSV.read("Data/control_enzyme.dat")
c_data =readdlm("Data/control_enzyme.dat")[2:end,1:19]
mean_control = c_data[1:2,:] * 0.001*60
error_control = c_data[3:4,:] .* 0.001*60

idx_plot = [1;2;16]
width = 0.5
row = 3
col = 1
a = collect(1:2:2*2)
fig,ax = plt[:subplots](col,row,figsize=(6,2.2))
for idx = 1:3
    ax[idx][:set_title](names(control)[idx_plot[idx]],fontsize=10)
    ax[idx][:set_ylabel]("Activity (mM/hr)")
    ax[idx][:set_xlabel]("Time (h)")
    b_control = ax[idx][:bar](a.-1,mean_control[:,idx_plot[idx]],yerr=error_control[:,idx_plot[idx]],align="center",color="lightgrey",ecolor="black",capsize=4)
    if idx == 2
        ax[idx][:set_yticks]([0,200,400,600,800,1000])
        ax[idx][:set_ylim]([0,1000])
    elseif idx == 3
        ax[idx][:set_yticks]([0,4,8,12])
        ax[idx][:set_ylim]([0,12])

    end
    ax[idx][:set_xticks]([])
    plt[:margins](0.015)
    plt[:tick_params](axis="x",which="both",bottom="off",top="off")

    tight_layout()
end
savefig("enzyme_1.pdf")


no_plot = [1;2;16]
all = collect(1:1:19)
idx_plot = Any[]
for i = 1:19
    idx = findall(all[i].==no_plot)
    if isempty(idx)
        push!(idx_plot,i)
    end
end


width = 0.6
row = 4
col = 4
a = collect(1:2:2*2)
fig,ax = plt[:subplots](col,row,figsize=(8,15))
for idx = 1:16
    ax[idx][:set_title](names(control)[idx_plot[idx]],fontsize=10)
    ax[idx][:set_ylabel]("Activity (mM/hr)")
    ax[idx][:set_xlabel]("Time (h)")
    b_control = ax[idx][:bar](a.-1, mean_control[:,idx_plot[idx]],yerr=error_control[:,idx_plot[idx]],align="center",color="lightgrey",ecolor="black",capsize=4)
    ax[idx][:set_xticks]([])
    plt[:margins](0.015)
    plt[:tick_params](axis="x",which="both",bottom="on",top="off")
    tight_layout()
end
savefig("enzyme_2.pdf")
