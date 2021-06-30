using PyPlot
using CSV
using DataFrames
t = [0;2;4;8;16.0]
row = 4
col = 3

control = CSV.read("Data/control_tca.dat")
mean_control = control[1:5,:]
error_control = control[6:10,:]

DNP = CSV.read("Data/DNP_tca.dat")
mean_dnp = DNP[1:5,:]
error_dnp = DNP[6:10,:]

TTA = CSV.read("Data/TTA_tca.dat")
mean_tta = TTA[1:5,:]
error_tta = TTA[6:10,:]
figure(figsize=(8,6))
for idx = 1:12
    subplot(row,col,idx)
    ylabel(names(control)[idx])
    errorbar(t, mean_tta[:,idx], yerr=error_tta[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:green")
    plot(t,mean_tta[:,idx],color="tab:green")
    errorbar(t, mean_dnp[:,idx], yerr=error_dnp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:red")
    plot(t,mean_dnp[:,idx],color="tab:red")
    errorbar(t, mean_control[:,idx], yerr=error_control[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="dimgrey")
    plot(t,mean_control[:,idx],color="dimgrey")
    axis([0,16.2,0,1])
    xticks([0,4,8,12,16])
    yticks([0,0.5,1])
    if names(control)[idx] == ("ACE (mM)")
        axis([0,16.2,0,60])
        yticks([0,20,40,60])
    elseif names(control)[idx] == ("LAC (mM)")
        axis([0,16.2,0,10])
        yticks([0,5,10])
    elseif names(control)[idx] == ("PYR (mM)")
        axis([0,16.2,0,15])
        yticks([0,5,10,15])
    elseif names(control)[idx] == ("SUCC (mM)")
        axis([0,16.2,0,20])
        yticks([0,5,10,15,20])
    elseif names(control)[idx] == ("CIT (mM)")
        axis([0,16.2,0,0.1])
        yticks([0,0.05,0.1])
    elseif names(control)[idx] == ("aKG (mM)")
        axis([0,16.2,0,2])
        yticks([0,1,2])
    elseif names(control)[idx] == ("OAA (mM)")
        axis([0,16.2,0,0.2])
        yticks([0,0.1,0.2])
    elseif names(control)[idx] == ("FUM (mM)")
        axis([0,16.2,0,3])
        yticks([0,1.5,3])
    elseif names(control)[idx] == ("MAL (mM)")
        axis([0,16.2,0,6])
        yticks([0,3,6])
    elseif names(control)[idx] == ("PEP (mM)")
        axis([0,16.2,0,6])
        yticks([0,3,6])
    end
    xlabel("Time (h)")
end
tight_layout()
savefig("plot_tca.pdf")
