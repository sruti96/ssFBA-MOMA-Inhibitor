using PyPlot
using CSV
using DataFrames
t = [0;2;4;8;16.0]
row = 5
col = 3

control = CSV.read("Data/control_gly.dat")
mean_control = control[1:5,:]
error_control = control[6:10,:]

DNP = CSV.read("Data/DNP_gly.dat")
mean_dnp = DNP[1:5,:]
error_dnp = DNP[6:10,:]

TTA = CSV.read("Data/TTA_gly.dat")
mean_tta = TTA[1:5,:]
error_tta = TTA[6:10,:]
figure(figsize=(8,8))
for idx = 1:13
    subplot(row,col,idx)
    ylabel(names(control)[idx])
    errorbar(t, mean_tta[:,idx], yerr=error_tta[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:green")
    plot(t,mean_tta[:,idx],color="tab:green")
    errorbar(t, mean_dnp[:,idx], yerr=error_dnp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:red")
    plot(t,mean_dnp[:,idx],color="tab:red")
    errorbar(t, mean_control[:,idx], yerr=error_control[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="dimgrey")
    plot(t,mean_control[:,idx],color="dimgrey")
    axis([0,16.2,0,0.05])
    xticks([0,4,8,12,16])
    #yticks([0,0.05,0.1])
    if names(control)[idx] == ("3PG (mM)")
        axis([0,16.2,0,8])
        yticks([0,4,8])
    elseif names(control)[idx] == ("F16P (mM)")
        axis([0,16.2,0,4])
        yticks([0,2,4])
    elseif names(control)[idx] == ("GAP (mM)")
        axis([0,16.2,0,7])
        yticks([0,3,6])
    elseif names(control)[idx] == ("Gly3P (mM)")
        axis([0,16.2,0,2])
        yticks([0,1,2])
    elseif names(control)[idx] == ("E4P (mM)")
        axis([0,16.2,0,2])
        yticks([0,1,2])
    elseif names(control)[idx] == ("RL5P (mM)")
        axis([0,16.2,0,0.2])
    elseif names(control)[idx] == ("R5P (mM)")
        axis([0,16.2,0,0.2])
    elseif names(control)[idx] == ("G6P (mM)")
        axis([0,16.2,0,0.1])
    elseif names(control)[idx] == ("Maltose (mM)")
        axis([0,16.2,0,0.5])
        #yticks([0,0.5,0.1])
    end
    xlabel("Time (h)")
end
tight_layout()
savefig("plot_gly.pdf")
