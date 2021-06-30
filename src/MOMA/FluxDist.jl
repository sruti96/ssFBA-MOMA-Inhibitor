using Statistics

case = "tta"
f1 = readdlm("Ensemble/$(case)/Flux_1.txt")
no_time, no_flux = size(f1)
no_samples = 100
flux_ensemble = zeros(no_time,no_flux,no_samples)
for i = 1:no_samples
    flux_ensemble[:,:,i] = readdlm("Ensemble/$(case)/Flux_$i.txt")
end

norm_flux_ensemble = zeros(no_time,no_flux,no_samples)
for i = 1:no_samples
    norm_flux_ensemble[:,:,i] = flux_ensemble[:,:,i]./1 .* 100
end

norm_flux_mean = mean(norm_flux_ensemble, dims= 3)
norm_flux_error = std(norm_flux_ensemble, dims= 3)

f_mean = mean(flux_ensemble, dims =3)
f_error = std(flux_ensemble, dims =3)

flux_ox = norm_flux_mean[:,63]
flux_ox_error = norm_flux_error[:,63]

writedlm("Results/Flux/$(case)_flux_ox.txt",flux_ox)
writedlm("Results/Flux/$(case)_flux_ox_err.txt",flux_ox_error)

flux_2 = norm_flux_mean[20,:]
flux_2_error = norm_flux_error[20,:]
flux_8 = norm_flux_mean[80,:]
flux_8_error = norm_flux_error[80,:]
flux_10 = norm_flux_mean[160,:]
flux_10_error = norm_flux_error[160,:]


reversible_fluxes = [
    9; #pgi f6p -> G6P
    13; #fbA DHAP + G3P -> F16P
    15; #tpiA g3p -> dhap
    19; #gapA 13dpg -> g3p
    21; #pgk 3pg -> 13dpg
    23; #gpm 2pg -> 3pg
    25; #eno pep -> 2pg
    32; #zwf 6pgl -> g6p
    36; #rpe xu5p -> ru5p
    38; #rpi ru5p -> r5p
    40; #talAB e4p + f6p -> g3p + s7p
    42; #tkt1 g3p + s7P ->r5p + xu5p
    44; #tkt2 f6p + g3p -> e4p + xu5p
    49; #acn icit -> cit
    51; #icd akg -> icit
    57; #fum mal -> fum
    59; #mdh oaa -> mal
    76; #pta actp -> accoa
    78; #ackA ac -> actp
    83; #lac -> pyr
    ];

flux_2_single = deepcopy(flux_2)
flux_8_single = deepcopy(flux_8)
flux_10_single = deepcopy(flux_10)

flux_2_e = deepcopy(flux_2_error)
flux_8_e = deepcopy(flux_8_error)
flux_10_e = deepcopy(flux_10_error)


for i = 1:length(reversible_fluxes)
    flux_2_single[reversible_fluxes[i]-1] = flux_2_single[reversible_fluxes[i]-1] - flux_2_single[reversible_fluxes[i]]
    flux_8_single[reversible_fluxes[i]-1] = flux_8_single[reversible_fluxes[i]-1] - flux_8_single[reversible_fluxes[i]]
    flux_10_single[reversible_fluxes[i]-1] = flux_10_single[reversible_fluxes[i]-1] - flux_10_single[reversible_fluxes[i]]

    flux_2_e[reversible_fluxes[i]-1] = flux_2_single[reversible_fluxes[i]-1] + flux_2_single[reversible_fluxes[i]]
    flux_8_e[reversible_fluxes[i]-1] = flux_8_single[reversible_fluxes[i]-1] + flux_8_single[reversible_fluxes[i]]
    flux_10_e[reversible_fluxes[i]-1] = flux_10_single[reversible_fluxes[i]-1] + flux_10_single[reversible_fluxes[i]]

end
STM = readdlm("config/Network.dat");

r5p_idx = [97;106;141]
e4p_idx = [126]
pg3_idx = [104]
pep_idx = [126]
pyr_idx = [85;98;99;100;108; 86;106;113;114;119;120;129]
aca_idx = [87;91;99;112;116]
oaa_idx = [88]

glu_idx = [95;103;123;131; 89;106;115;122;128;142;147;148;151;160]
akg_idx = [85;87;88;94;98;99;100;102;104;107;108; 86;92;93;109;124;125]

akg_s_idx = STM[21,akg_idx]
glu_s_idx = STM[64,glu_idx]
r5p_s_idx = STM[125,r5p_idx]
e4p_s_idx = STM[46,e4p_idx]
pg3_s_idx = STM[6,pg3_idx]
pep_s_idx = STM[113,pep_idx]
pyr_s_idx = STM[122,pyr_idx]
aca_s_idx = STM[16,aca_idx]
oaa_s_idx = STM[110,oaa_idx]


glu_con = sum(glu_s_idx.*flux_2_single[glu_idx])
akg_con = sum(akg_s_idx.*flux_2_single[akg_idx])
r5p_con = sum(r5p_s_idx.*flux_2_single[r5p_idx])
e4p_con = sum(e4p_s_idx.*flux_2_single[e4p_idx])
pg3_con = sum(pg3_s_idx.*flux_2_single[pg3_idx])
pep_con = sum(pep_s_idx.*flux_2_single[pep_idx])
pyr_con = sum(pyr_s_idx.*flux_2_single[pyr_idx])
aca_con = sum(aca_s_idx.*flux_2_single[aca_idx])
oaa_con = sum(oaa_s_idx.*flux_2_single[oaa_idx])

glu_con_8 = sum(glu_s_idx.*flux_8_single[glu_idx])
akg_con_8 = sum(akg_s_idx.*flux_8_single[akg_idx])
r5p_con_8 = sum(r5p_s_idx.*flux_8_single[r5p_idx])
e4p_con_8 = sum(e4p_s_idx.*flux_8_single[e4p_idx])
pg3_con_8 = sum(pg3_s_idx.*flux_8_single[pg3_idx])
pep_con_8 = sum(pep_s_idx.*flux_8_single[pep_idx])
pyr_con_8 = sum(pyr_s_idx.*flux_8_single[pyr_idx])
aca_con_8 = sum(aca_s_idx.*flux_8_single[aca_idx])
oaa_con_8 = sum(oaa_s_idx.*flux_8_single[oaa_idx])

glu_con_10 = sum(glu_s_idx.*flux_10_single[glu_idx])
akg_con_10 = sum(akg_s_idx.*flux_10_single[akg_idx])
r5p_con_10 = sum(r5p_s_idx.*flux_10_single[r5p_idx])
e4p_con_10 = sum(e4p_s_idx.*flux_10_single[e4p_idx])
pg3_con_10 = sum(pg3_s_idx.*flux_10_single[pg3_idx])
pep_con_10 = sum(pep_s_idx.*flux_10_single[pep_idx])
pyr_con_10 = sum(pyr_s_idx.*flux_10_single[pyr_idx])
aca_con_10 = sum(aca_s_idx.*flux_10_single[aca_idx])
oaa_con_10 = sum(oaa_s_idx.*flux_10_single[oaa_idx])

aa_up = [glu_con; akg_con; r5p_con; e4p_con; pg3_con; pep_con; pyr_con; aca_con; oaa_con]
aa_up_8 = [glu_con_8; akg_con_8; r5p_con_8; e4p_con_8; pg3_con_8; pep_con_8; pyr_con_8; aca_con_8; oaa_con_8]
aa_up_10 = [glu_con_10; akg_con_10; r5p_con_10; e4p_con_10; pg3_con_10; pep_con_10; pyr_con_10; aca_con_10; oaa_con_10]


flux_2_single[77] = flux_2_single[77] - flux_2_single[79] # ACA -> AC
flux_2_single[29] = flux_2_single[29] + flux_2_single[84] # PYR -> ACA
flux_2_single[73] = flux_2_single[73] + flux_2_single[74] #MAL -> PYR

flux_8_single[77] = flux_8_single[77] - flux_8_single[79] # ACA -> AC
flux_8_single[29] = flux_8_single[29] + flux_8_single[84] # PYR -> ACA
flux_8_single[73] = flux_8_single[73] + flux_8_single[74] #MAL -> PYR

flux_10_single[77] = flux_10_single[77] - flux_10_single[79] # ACA -> AC
flux_10_single[29] = flux_10_single[29] + flux_10_single[84] # PYR -> ACA
flux_10_single[73] = flux_10_single[73] + flux_10_single[74] #MAL -> PYR

idx = [
    1; #gp  M -> Maltose 1
    6; #pgm G1p -> G6P 2
    7; #hk Glc -> G6P 3
    8; #pgi G6P -> F6P 4
    10; #pfk F6P -> F16P 5
    11; #fdp F16P -> F6P 6
    12; #fbaA F16P -> DHAP + G3P 7
    14; #tpiA DHAP -> G3P 8
    18; #gapA G3P -> 13DPG 9
    20; #pgk 13DPG -> 3PG 10
    22; #gpm 3PG -> 2PG 11
    24; #eno 2PG -> 3PG 12
    26; #pyk PEP -> PYR
    27; #pck OAA -> PEP
    28; #ppc PEP -> OAA 15
    29; #pdh PYR -> ACA
    30; #pps PYR -> PEP
    31; #zwf G6P -> 6PGL
    33; #pgl 6PGL -> 6PGC
    34; #gnd 6PGC -> ru5p 20
    35; #rpe ru5P -> xu5P
    37; #rpi r5p -> ru5P
    39; #talAB g3p + s7P -> e4p + f6p
    41; #tkt1 r5p + xu5p -> g3p + s7p
    43; #tkt2 e4p + xu5p -> f6p + g3p 25
    45; #edd 6pgc -> 2kddg6p
    46; #eda 2ddg6p -> g3p + pyr
    47; #gltA ACA + OAA -> cit
    48; #acn cit -> icit
    50; #icd icit -> akg 30
    52; #sucAB akg -> succoa
    53; #sucCD succoa -> succ
    54; #sdh succ -> fum
    55; #frd fum -> succ
    56; #fum fum -> mal 35
    58; #mdh mal -> oaa
    71; #aceA icit -> glx + succ
    72; #aceB ACA + glx -> mal
    73; #maeAB mal -> pyr
    77; #ackA actp -> ac 40
    82; #ldh pyr -> lac
    63; #atp
    ]

f2 = flux_2_single[idx]
f8 = flux_8_single[idx]
f10 = flux_10_single[idx]

append!(f2,aa_up)
append!(f8,aa_up_8)
append!(f10,aa_up_10)

f2_error = flux_2_error[idx]
f8_error = flux_8_error[idx]
f10_error = flux_10_error[idx]

writedlm("Results/Flux/$(case)_flux_2h.txt",f2)
writedlm("Results/Flux/$(case)_flux_8h.txt",f8)
writedlm("Results/Flux/$(case)_flux_10h.txt",f10)

writedlm("Results/Flux/$(case)_flux_2h_error.txt",f2_error)
writedlm("Results/Flux/$(case)_flux_8h_error.txt",f8_error)
writedlm("Results/Flux/$(case)_flux_10h_error.txt",f10_error)
