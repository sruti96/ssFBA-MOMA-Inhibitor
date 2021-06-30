using PyPlot
using CSV
using DataFrames
t = [0;2;4;8;16.0]
row = 5
col = 4

control = CSV.read("Data/control_aa.dat")
mean_control = control[1:5,:]
error_control = control[6:10,:]

DNP = CSV.read("Data/DNP_aa.dat")
mean_dnp = DNP[1:5,:]
error_dnp = DNP[6:10,:]

TTA = CSV.read("Data/TTA_aa.dat")
mean_tta = TTA[1:5,:]
error_tta = TTA[6:10,:]
figure(figsize=(10,8))
for idx = 1:20
    subplot(row,col,idx)
    ylabel(names(control)[idx])
    errorbar(t, mean_tta[:,idx], yerr=error_tta[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:green")
    plot(t,mean_tta[:,idx],color="tab:green")
    errorbar(t, mean_dnp[:,idx], yerr=error_dnp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:red")
    plot(t,mean_dnp[:,idx],color="tab:red")
    errorbar(t, mean_control[:,idx], yerr=error_control[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="dimgrey")
    plot(t,mean_control[:,idx],color="dimgrey")
    axis([0,16.2,0,3])
    xticks([0,4,8,12,16])
    yticks([0,1,2,3])
    if names(control)[idx] == "ALA (mM)"
        axis([0,16.2,0,10])
        yticks([0,5,10])
    elseif names(control)[idx] == "GLU (mM)"
        axis([0,16.2,0,150])
        yticks([0,50,100,150])
    elseif names(control)[idx] == "VAL (mM)"
        axis([0,16.2,0,10])
        yticks([0,5,10])
    elseif names(control)[idx] == "SER (mM)"
        axis([0,16.2,0,4])
        yticks([0,2,4])
    elseif names(control)[idx] == "TRP (mM)"
        axis([0,16.2,0,4])
        yticks([0,2,4])
    elseif names(control)[idx] == "CYS (mM)"
        axis([0,16.2,0,1])
        yticks([0,0.5,1])
    elseif names(control)[idx] == "GLN (mM)"
        axis([0,16.2,0,1])
        yticks([0,0.5,1])
    end
    xlabel("Time (h)")

end
tight_layout()
savefig("plot_aa.pdf")
