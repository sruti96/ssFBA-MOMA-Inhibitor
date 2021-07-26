include("Include.jl")

function calculate_flux_difference(fba_flux, moma_flux)
    sum = 0
    for i = 1:length(fba_flux)
        sum += (fba_flux[i] - moma_flux[i])^2
    end
    l2_norm = (sum)^0.5
    return l2_norm
end

function calculate_normalized_flux(case,type)

    f1 = readdlm("/Users/srutidammalapati/Desktop/ssFBA-MOMA-Inhibitor/src/$(type)/Ensemble/$(case)/Flux_1.txt")
    no_time, no_flux = size(f1)
    no_samples = 100
    flux_ensemble = zeros(no_time,no_flux,no_samples)
    for i = 1:no_samples
        flux_ensemble[:,:,i] = readdlm("/Users/srutidammalapati/Desktop/ssFBA-MOMA-Inhibitor/src/$(type)/Ensemble/$(case)/Flux_$i.txt")
    end

    norm_flux_ensemble = zeros(no_time,no_flux,no_samples)
    for i = 1:no_samples
        norm_flux_ensemble[:,:,i] = flux_ensemble[:,:,i]./flux_ensemble[1,1,i] .* 100
    end

    norm_flux_mean = mean(norm_flux_ensemble, dims= 3)
    norm_flux_error = std(norm_flux_ensemble, dims= 3)

    return norm_flux_mean
end

case = "dnp"

fba_f_mean = calculate_normalized_flux(case, "ssFBA")
moma_f_mean = calculate_normalized_flux(case, "MOMA")
diff_at_i = zeros(160)

for i = 1:160

    fba_at_i = fba_f_mean[i,:]
    moma_at_i = moma_f_mean[i,:]
    diff_at_i[i] = calculate_flux_difference(fba_at_i, moma_at_i)
end
diff_at_i = diff_at_i./maximum(diff_at_i)
plot(diff_at_i)
