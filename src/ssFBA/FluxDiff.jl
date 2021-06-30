#Find relative difference of fluxes between treatment and control
using PyPlot
idx = [2; 3; 4; 11; 12; 13; 16; 18; 21; 22; 24; 28; 31; 33; 35; 36; 39; 40; 41; 42]#;32;30]
scale = 751.1
f2_control = readdlm("Results/Flux/control_flux_2h.txt")[idx] ./scale
f2_control_err = readdlm("Results/Flux/control_flux_2h_error.txt")[idx] ./scale .*100

f8_control = readdlm("Results/Flux/control_flux_8h.txt")[idx] ./scale
f8_control_err = readdlm("Results/Flux/control_flux_8h_error.txt")[idx] ./scale .*100

f10_control = readdlm("Results/Flux/control_flux_10h.txt")[idx] ./scale
f10_control_err = readdlm("Results/Flux/control_flux_10h_error.txt")[idx] ./scale .*100

f2_dnp = readdlm("Results/Flux/dnp_flux_2h.txt")[idx] ./scale
f2_dnp_err = readdlm("Results/Flux/dnp_flux_2h_error.txt")[idx] ./scale .*100

f8_dnp = readdlm("Results/Flux/dnp_flux_8h.txt")[idx] ./scale
f8_dnp_err = readdlm("Results/Flux/dnp_flux_8h_error.txt")[idx] ./scale .*100

f10_dnp = readdlm("Results/Flux/dnp_flux_10h.txt")[idx] ./scale
f10_dnp_err = readdlm("Results/Flux/dnp_flux_10h_error.txt")[idx] ./scale .*100

f2_tta = readdlm("Results/Flux/tta_flux_2h.txt")[idx] ./scale
f2_tta_err = readdlm("Results/Flux/tta_flux_2h_error.txt")[idx] ./scale .*100

f8_tta = readdlm("Results/Flux/tta_flux_8h.txt")[idx] ./scale
f8_tta_err = readdlm("Results/Flux/tta_flux_8h_error.txt")[idx] ./scale .*100

f10_tta = readdlm("Results/Flux/tta_flux_10h.txt")[idx] ./scale
f10_tta_err = readdlm("Results/Flux/tta_flux_10h_error.txt")[idx] ./scale .*100

#f2_dnp_err = f2_dnp_err ./ f2_control .*100
f2_dnp_norm = (f2_dnp .- f2_control) .* 100
f8_dnp_norm = (f8_dnp .- f8_control) .*100
f10_dnp_norm = (f10_dnp .- f10_control) .* 100
f2_tta_norm = (f2_tta .- f2_control) .* 100
f8_tta_norm = (f8_tta .- f8_control) .* 100
f10_tta_norm = (f10_tta .- f10_control) .* 100 #./ f10_control .* 100
a = collect(1:1:length(f2_dnp_norm))

#Choose which flux diff to plot
exp = vec(f8_dnp_norm)
exp_err = vec(f8_dnp_err)

fig,ax = plt[:subplots](figsize=(6,3))
f_pos = findall(exp.>0.0)
f_neg = findall(exp.<0.0)
dnp = ax[:bar](a[f_pos],exp[f_pos],yerr =exp_err[f_pos],align="center",color="#5a8fe6",ecolor="black",capsize=2.5)
dnp = ax[:bar](a[f_neg],exp[f_neg],yerr = exp_err[f_neg],align="center",color="#d64f4f",ecolor="black",capsize=2.5)

ax[:set_ylabel]("Flux difference from control (A.U.)")
ax[:set_xticks]([])
plt[:margins](0.015)
plt[:ylim]([-75,25])
ax[:set_yticks]([-75,-50,-25,0,25])

#plt[:ylim]([-105,500])
#ax[:set_yticks]([-100,-50,0,50])
plt[:tick_params](axis="x",which="both",bottom="off",top="off")
tight_layout()

#=====================Normalize Flux to scalar for flux distribution figure ==============#
f2_control = readdlm("Results/Flux/control_flux_2h.txt") ./scale .* 100
f10_control = readdlm("Results/Flux/control_flux_10h.txt") ./scale .* 100
f8_control = readdlm("Results/Flux/control_flux_8h.txt") ./scale .* 100

f2_dnp = readdlm("Results/Flux/dnp_flux_2h.txt") ./scale .* 100
f10_dnp = readdlm("Results/Flux/dnp_flux_10h.txt") ./scale .* 100
f8_dnp = readdlm("Results/Flux/dnp_flux_8h.txt") ./scale .* 100


f2_tta = readdlm("Results/Flux/tta_flux_2h.txt") ./scale .* 100
f10_tta = readdlm("Results/Flux/tta_flux_10h.txt") ./scale .* 100
f8_tta = readdlm("Results/Flux/tta_flux_8h.txt") ./scale .* 100


writedlm("Results/NormalizedFlux_control_2h_fig.txt",f2_control)
writedlm("Results/NormalizedFlux_control_10h_fig.txt",f10_control)
writedlm("Results/NormalizedFlux_control_8h_fig.txt",f8_control)

writedlm("Results/NormalizedFlux_dnp_2h_fig.txt",f2_dnp)
writedlm("Results/NormalizedFlux_dnp_10h_fig.txt",f10_dnp)
writedlm("Results/NormalizedFlux_dnp_8h_fig.txt",f8_dnp)

writedlm("Results/NormalizedFlux_tta_2h_fig.txt",f2_tta)
writedlm("Results/NormalizedFlux_tta_10h_fig.txt",f10_tta)
writedlm("Results/NormalizedFlux_tta_8h_fig.txt",f8_tta)
