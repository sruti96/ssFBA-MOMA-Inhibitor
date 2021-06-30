using PyPlot
using CSV
using DataFrames
t = [0;2;4;8;16.0]
row = 2
col = 3

control = CSV.read("Data/control_reducing.dat")
mean_control = control[1:5,:]
error_control = control[6:10,:]

DNP = CSV.read("Data/DNP_reducing.dat")
mean_dnp = DNP[1:5,:]
error_dnp = DNP[6:10,:]

TTA = CSV.read("Data/TTA_reducing.dat")
mean_tta = TTA[1:5,:]
error_tta = TTA[6:10,:]
figure(figsize=(8,4))
for idx = 1:5
    subplot(row,col,idx)
    ylabel(names(control)[idx])
    errorbar(t, mean_tta[:,idx], yerr=error_tta[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:green")
    plot(t,mean_tta[:,idx],color="tab:green")
    errorbar(t, mean_dnp[:,idx], yerr=error_dnp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:red")
    plot(t,mean_dnp[:,idx],color="tab:red")
    errorbar(t, mean_control[:,idx], yerr=error_control[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="dimgrey")
    plot(t,mean_control[:,idx],color="dimgrey")
    axis([0,16.2,0,0.2])
    xticks([0,4,8,12,16])
    if names(control)[idx] == ("NADH (mM)")
        axis([0,16.2,0,2])

    elseif names(control)[idx] == ("FAD (mM)")
        axis([0,16.2,0,0.1])
    elseif names(control)[idx] == ("NAD (mM)")
        axis([0,16.2,0,0.5])
    end
    xlabel("Time (h)")

end
tight_layout()
savefig("plot_reducing.pdf")
