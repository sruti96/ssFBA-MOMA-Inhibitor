using DelimitedFiles
using Statistics

case = "tta"
CarbonNum = readdlm("config/SpeciesDict/CarbonNum.txt")

S1 = readdlm("Ensemble/$(case)/Species_1.txt")
no_time, no_species = size(S1)
no_samples = 100
Species_ensemble = zeros(no_time,no_species,no_samples)
for i = 1:no_samples
    Species_ensemble[:,:,i] = readdlm("Ensemble/$(case)/Species_$i.txt")
end
net = zeros(no_species,no_samples)
for j = 1:no_samples
    net[:,j] = Species_ensemble[end,:,j] - Species_ensemble[1,:,j]
end
sub_idx = findall(net.<0.0) #Substrates (Depleted)
prod_idx = findall(net.>0.0) #Products (Accumulated)
substrate = zeros(no_species,no_samples)
product = zeros(no_species,no_samples)
exclude_energy = [18; 24;31;34;38;41;60;71;72;140;141;142] #Exclude Enery Species from Carbon Calculation
CarbonNum[exclude_energy] .= 0.0
substrate[sub_idx] = abs.(net[sub_idx])
product[prod_idx] = net[prod_idx]

#Calculate Total Carbon consumed and produced
substrate_carbon = substrate.*CarbonNum
product_carbon = product.*CarbonNum

#Divide into sections
organic = [15; 16; 17; 21; 36; 47; 53; 54; 66; 81; 86; 91; 110; 122; 131; 132; 40]
gly = [3; 5; 6; 44;; 48; 50; 55; 56; 57; 61; 92; 93; 94; 95; 96; 113]
ppp = [11; 12; 46; 125; 126; 127; 146]
energy = [18; 24; 31; 34; 38; 41; 60; 71; 72; 140; 141; 142]

aa = Any[]
aa_index = convert(Array{Int64},readdlm("config/SpeciesDict/aa_index.dat")) #Amino Acids
aa_trna = aa_index .+ 1 #Amino Acids charged tRNA
append!(aa,aa_index)
append!(aa,aa_trna)

all = Any[]
append!(all,aa)
append!(all,organic)
append!(all,gly)
append!(all,energy)
append!(all,ppp)

remain = collect(1:1:157)
other = Any[]
for i = 1:157
    idx = findall(remain[i].==all)
    if isempty(idx)
        push!(other,i)
    end
end

CarbonYield = Any[]
AAYield = Any[]
CO2Yield = Any[]
OrgYield = Any[]
GlyYield = Any[]
pppYield = Any[]
energyYield = Any[]
otherYield = Any[]
totalYield = Any[]
for i = 1:no_samples
    push!(CarbonYield, product_carbon[148,i]/sum(substrate_carbon[:,i]))
    push!(AAYield, sum(product_carbon[aa_index,i])/sum(substrate_carbon[:,i]))
    push!(CO2Yield, product_carbon[39,i]/sum(substrate_carbon[:,i]))
    push!(OrgYield, sum(product_carbon[organic,i])/sum(substrate_carbon[:,i]))
    push!(GlyYield, sum(product_carbon[gly,i])/sum(substrate_carbon[:,i]))
    push!(pppYield, sum(product_carbon[ppp,i])/sum(substrate_carbon[:,i]))
    push!(energyYield, sum(product_carbon[energy,i])/sum(substrate_carbon[:,i]))
    push!(otherYield, sum(product_carbon[other,i])/sum(substrate_carbon[:,i]))
    push!(totalYield, sum(product_carbon[:,i])/sum(substrate_carbon[:,i]))
end

yield_mean = mean(CarbonYield)*100
yield_error = std(CarbonYield)*100


aa_yield_mean = mean(AAYield)*100
aa_yield_error = std(AAYield)*100
co2_yield_mean = mean(CO2Yield)*100
co2_yield_error = std(CO2Yield)*100
org_yield_mean = mean(OrgYield)*100
org_yield_error = std(OrgYield)*100
gly_yield_mean = mean(GlyYield)*100
gly_yield_error = std(GlyYield)*100
ppp_yield_mean = mean(pppYield)*100
ppp_yield_error = std(pppYield)*100
energy_yield_mean = mean(energyYield)*100
energy_yield_error = std(energyYield)*100
other_yield_mean = mean(otherYield)*100
other_yield_error = std(otherYield)*100
total_yield_mean = mean(totalYield)*100
total_yield_error = std(totalYield)*100
println("CarbonYield = ",yield_mean, " +/- ", yield_error)
println("CO2Yield = ",co2_yield_mean, " +/- ", co2_yield_error)
println("OrgYield = ",org_yield_mean, " +/- ", org_yield_error)
println("GlyYield = ",gly_yield_mean, " +/- ", gly_yield_error)
println("PPPYield = ",ppp_yield_mean, " +/- ", ppp_yield_error)
println("AAYield = ",aa_yield_mean, " +/- ", aa_yield_error)
println("EnergyYield = ",energy_yield_mean, " +/- ", energy_yield_error)
println("OtherYield = ",other_yield_mean, " +/- ", other_yield_error)
println("TotalYield = ",total_yield_mean, " +/- ", total_yield_error)

using DataFrames
using CSV
df = DataFrame()
df.A = ["Protein","CO2","Organic", "Glycolysis", "PPP", "AA", "Energy", "Other"]
df.Mean = [yield_mean, co2_yield_mean, org_yield_mean, gly_yield_mean, ppp_yield_mean, aa_yield_mean, energy_yield_mean, other_yield_mean]
df.Error =  [yield_error, co2_yield_error, org_yield_error, gly_yield_error, ppp_yield_error, aa_yield_error, energy_yield_error, other_yield_error]
CSV.write("Results/CarbonYield_res_$case.csv",df)
