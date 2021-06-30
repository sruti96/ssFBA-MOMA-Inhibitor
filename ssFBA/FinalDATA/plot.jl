using PyPlot
using CSV
using DataFrames
t = [0;2;4;8;16.0]
row = 1
col = 2

control = CSV.read("Data/control_mRNA_protein.dat")
mean_control = control[1:5,:]
error_control = control[6:10,:]

DNP = CSV.read("Data/DNP_mRNA_protein.dat")
mean_dnp = DNP[1:5,:]
error_dnp = DNP[6:10,:]

TTA = CSV.read("Data/TTA_mRNA_protein.dat")
mean_tta = TTA[1:5,:]
error_tta = TTA[6:10,:]
figure(1,figsize=(6.5,2.5))
for idx = 1:2
    subplot(row,col,idx)
    ylabel(names(control)[idx])
    # errorbar(t, mean_tta[:,idx], yerr=error_tta[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:green")
    # plot(t,mean_tta[:,idx],color="tab:green")
    # errorbar(t, mean_dnp[:,idx], yerr=error_dnp[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="tab:red")
    # plot(t,mean_dnp[:,idx],color="tab:red")
    errorbar(t, mean_control[:,idx], yerr=error_control[:,idx],fmt="o",markersize = 4,capsize=2,elinewidth=1,color="dimgrey")
    plot(t,mean_control[:,idx],color="dimgrey")
    axis([0,16.2,0,30])
    xticks([0,4,8,12,16])
    ylabel("GFP Protein (μM)")

    if names(control)[idx] == "mRNA (nM)"
        axis([0,16.2,0,1000])
        yticks([0,200,400,600,800,1000])
        ylabel("GFP mRNA (nM)")
    end
    xlabel("Time (h)")
end
tight_layout()
savefig("protein_mrna.pdf")
