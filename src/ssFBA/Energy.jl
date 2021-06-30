using DelimitedFiles
using Statistics
network = readdlm("config/Network.dat")
prod_network = deepcopy(network)
cons_network = deepcopy(network)
for i = 1:223
    idx = findall(prod_network[:,i].<0.0)
    idx_c = findall(cons_network[:,i].>0.0)
    if ~isempty(idx)
        prod_network[idx,i] .= 0
    end
    if ~isempty(idx_c)
        cons_network[idx_c,i] .= 0
    end
end
no_species,no_flux = size(network)

case = "control"
f1 = readdlm("Ensemble/$(case)/Flux_1.txt")
no_time, no_flux = size(f1)
no_samples = 100

flux_ensemble = zeros(no_time,no_flux,no_samples)
for i = 1:no_samples
    flux_ensemble[:,:,i] = readdlm("Ensemble/$(case)/Flux_$i.txt")
end

#ATP phosphate bond loss
gly = [7; 10; 21]
gly_idx = [1;1;1]
ana = [27; 30; 78; 79]
ana_idx = [1; 2; 1; 2]
aa = [87; 89; 90; 95; 97; 100; 101; 103; 105; 106]
aa_idx = [4; 2; 2; 1; 5; 1; 1; 1; 2; 2]
chor = [126;  131]
chor_idx = [1;1]
purine = [141; 142; 146; 147; 149; 151; 152; 153; 154; 160]
purine_idx = [2; 2; 1; 1; 1; 1; 1; 1; 1; 2]
deg = [161; 165; 169; 170; 171; 172; 173; 174; 175]
deg_idx = [2; 1; 1; 1; 1; 1; 1; 1; 1]

#GTP bond loss
g_chor = [127]
g_chor_idx = [3]
g_purine = [158]
g_purine_idx = [1]
g_deg = [164; 168]
g_deg_idx = [2; 1]

#CTP_bond loss
c_deg = [163; 167]
c_deg_idx = [2; 1]

#UTP loss
utp_deg = [162; 166]
utp_deg_idx = [2;1]

#ATP gain
a_gly = [20; 26;]
a_gly_idx = [1; 1]
a_tca = [53]
a_tca_idx = [1]
a_ox = [63]
a_ox_idx = [1]
a_overflow = [77]
a_overflow_idx = [1]
a_aa =[111; 118]
a_aa_idx =[2; 1]

n_gly = zeros(no_samples,1)
n_tca = zeros(no_samples,1)
n_ox = zeros(no_samples,1)
n_overflow = zeros(no_samples,1)
n_aa = zeros(no_samples,1)

u_gly = zeros(no_samples,1)
u_ana = zeros(no_samples,1)
u_aa = zeros(no_samples,1)
u_chor = zeros(no_samples,1)
u_purine = zeros(no_samples,1)
u_deg = zeros(no_samples,1)

u_gtp_chor  = zeros(no_samples,1)
u_gtp_purine = zeros(no_samples,1)
u_gtp_deg = zeros(no_samples,1)
u_ctp_deg  = zeros(no_samples,1)
u_utp_deg  = zeros(no_samples,1)

dt = 0.1
for i = 1:no_samples

    n_gly[i] = sum((flux_ensemble[:,a_gly,i]*a_gly_idx .*dt))
    n_tca[i] = sum(flux_ensemble[:,a_tca,i]*a_tca_idx .*dt)
    n_ox[i] = sum(flux_ensemble[:,a_ox,i]*a_ox_idx .*dt)
    n_overflow[i] = sum(flux_ensemble[:,a_overflow,i]*a_overflow_idx .*dt)
    n_aa[i] = sum(flux_ensemble[:,a_aa,i]*a_aa_idx .*dt)

    u_gly[i] = sum(flux_ensemble[:,gly,i]*gly_idx .*dt)
    u_ana[i] = sum(flux_ensemble[:,ana,i]*ana_idx .*dt)
    u_aa[i] = sum(flux_ensemble[:,aa,i]*aa_idx .*dt)
    u_chor[i] = sum(flux_ensemble[:,chor,i]*chor_idx .*dt)
    u_purine[i] = sum(flux_ensemble[:,purine,i]*purine_idx .*dt)
    u_deg[i] = sum(flux_ensemble[:,deg,i]*deg_idx .*dt)

    u_gtp_chor[i] = sum(flux_ensemble[:,g_chor,i]*g_chor_idx .*dt)
    u_gtp_purine[i] = sum(flux_ensemble[:,g_purine,i]*g_purine_idx .*dt)

    u_gtp_deg[i] = sum(flux_ensemble[:,g_deg,i]*g_deg_idx .*dt)
    u_ctp_deg[i] = sum(flux_ensemble[:,c_deg,i]*c_deg_idx .*dt)
    u_utp_deg[i] = sum(flux_ensemble[:,utp_deg,i]*utp_deg_idx .*dt)


end

#====================ATP in TXTL =======================#
tl_amt = Any[]
for j= 181:200
    aa_idx = findall(network[:,j].>0.0)[1]
    push!(tl_amt, network[aa_idx,j])
end

TX_cost = sum(2*[195; 157; 157; 208])
TL_cost = sum(2*tl_amt) + 478

ProdEff = Any[]

for i = 1:no_samples
    push!(ProdEff, TX_cost*sum(flux_ensemble[:,177,i]*dt)+TL_cost*sum(flux_ensemble[:,180,i].*dt))
end

atp_prod = n_gly + n_tca + n_ox + n_overflow + n_aa #ATP Produced

protein_eff = ProdEff./atp_prod

gly_eff = u_gly./atp_prod

ana_eff = u_ana./atp_prod

aa_eff = u_aa./atp_prod

chor_eff = (u_chor+u_gtp_chor)./atp_prod

purine_eff = (u_purine+u_gtp_purine)./atp_prod

deg_eff = (u_deg+u_gtp_deg+u_ctp_deg+u_utp_deg)./atp_prod


protein_eff_mean = mean(protein_eff)*100
protein_eff_std = std(protein_eff)*100

gly_eff_mean = mean(gly_eff)*100
gly_eff_std = std(gly_eff)*100

ana_eff_mean = mean(ana_eff)*100
ana_eff_std = std(ana_eff)*100

aa_eff_mean = mean(aa_eff)*100
aa_eff_std = std(aa_eff)*100

chor_eff_mean = mean(chor_eff)*100
chor_eff_std = std(chor_eff)*100


purine_eff_mean = mean(purine_eff)*100
purine_eff_std = std(purine_eff)*100

deg_eff_mean = mean(deg_eff)*100
deg_eff_std = std(deg_eff)*100


total_eff = protein_eff_mean + gly_eff_mean + ana_eff_mean +chor_eff_mean+ purine_eff_mean+deg_eff_mean
println("EnergyEff = ",protein_eff_mean, " +/- ", protein_eff_std)
println("Gly = ",gly_eff_mean, " +/- ", gly_eff_std)
println("Anaplerosis = ",ana_eff_mean, " +/- ", ana_eff_std)
println("AminoAcid = ",aa_eff_mean, " +/- ", aa_eff_std)
println("Chorismate = ",chor_eff_mean, " +/- ", chor_eff_std)
println("Purine = ",purine_eff_mean, " +/- ", purine_eff_std)
println("Degradation = ",deg_eff_mean, " +/- ", deg_eff_std)
println("TotalEff = ",total_eff)



using DataFrames
df = DataFrame()
df.A = ["TXTL", "Glycolysis", "Anaplerosis", "AminoAcids", " Degradation", "Chorismate", "Purine"]
df.Mean = [protein_eff_mean, gly_eff_mean, ana_eff_mean, aa_eff_mean,deg_eff_mean, chor_eff_mean, purine_eff_mean ]
df.Error = [protein_eff_std, gly_eff_std, ana_eff_std, aa_eff_std,deg_eff_std, chor_eff_std, purine_eff_std ]
CSV.write("Results/EnergyEff_res_$case.csv",df)
