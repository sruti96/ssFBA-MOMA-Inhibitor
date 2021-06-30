using PyPlot
using CSV
using DataFrames
t = [0;2;4;8;16.0]
row = 3
col = 4

control = CSV.read("Data/control_energy.dat")
mean_control = control[1:5,:]
error_control = control[6:10,:]

DNP = CSV.read("Data/dnp_energy.dat")
mean_dnp = DNP[1:5,:]
error_dnp = DNP[6:10,:]

TTA = CSV.read("Data/tta_energy.dat")
mean_tta = TTA[1:5,:]
error_tta = TTA[6:10,:]
figure(1,figsize=(10,5))
for idx = 1:12
    subplot(row,col,idx)
    ylabel(names(control)[idx])
    errorbar(t, mean_tta[:,idx], yerr=error_tta[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:green")
    plot(t,mean_tta[:,idx],color="tab:green")
    errorbar(t, mean_dnp[:,idx], yerr=error_dnp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:red")
    plot(t,mean_dnp[:,idx],color="tab:red")
    errorbar(t, mean_control[:,idx], yerr=error_control[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="dimgrey")
    plot(t,mean_control[:,idx],color="dimgrey")
    axis([0,16.2,0,1.0])
    xticks([0,4,8,12,16])

    if names(control)[idx] == "ATP (mM)"
        axis([0,16.2,0,2])
        yticks([0,1,2])
    elseif names(control)[idx] == "ADP (mM)"
        axis([0,16.2,0,2.5])
        yticks([0,1,2])
    elseif names(control)[idx] == "UTP (mM)"
        axis([0,16.2,0,2])
        yticks([0,1,2])
    elseif names(control)[idx] == "AMP (mM)"
        axis([0,16.2,0,2])
        yticks([0,1,2])
    end
    xlabel("Time (h)")

end
tight_layout()
savefig("plot_energy.pdf")
