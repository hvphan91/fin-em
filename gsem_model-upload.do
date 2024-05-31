/// prep data

// cd...
// capture log close 
// log using [log_file], replace text


import excel fin-em-clean, firstrow case(lower)


for var co gdp pop wpop: gen lX = ln(X+0.01)
for var fin_fi fin_re fire gdpg m3 pcre cr2p sval scap manva ren gini indem trd urb: replace X = X/100 

gen p9950 = p99/p50
gen p9050 = p90/p50
gen lman = ln(manva*gdp + 0.01)

encode country, gen(cid)
xtset cid year


drop if inc == 2
drop if year < 1995 | year > 2020


/// fixed effect models
quietly xtreg lco l.fire l.lco lgdp p9950 lpop trd urb i.year, fe vce(cluster cid)
est store fire
quietly xtreg lco l.pfp l.lco lgdp p9950 lpop trd urb i.year, fe vce(cluster cid)
est store pfp
quietly xtreg lco l.pcre l.lco lgdp p9950 lpop trd urb i.year, fe vce(cluster cid)
est store pcre
quietly xtreg lco gdpg ren l.fire l.lco lgdp p9950 lpop trd urb i.year, fe vce(cluster cid)
est store fe_fire
quietly xtreg lco gdpg ren l.pfp l.lco lgdp p9950 lpop trd urb i.year, fe vce(cluster cid)
est store fe_pfp
quietly xtreg lco gdpg ren l.pcre l.lco lgdp p9950 lpop trd urb i.year, fe vce(cluster cid)
est store fe_pcre

esttab fire pfp pcre fe_fire fe_pfp fe_pcre, obslast l se noomit ar2 drop(*.year) nonum o(l.fire l.pfp l.pcre gdpg ren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst


quietly xtreg lco l.fire l.lco lgdp i.year, fe vce(cluster cid)
est store fire_s
quietly xtreg lco l.pfp l.lco lgdp i.year, fe vce(cluster cid)
est store pfp_s
quietly xtreg lco l.pcre l.lco lgdp i.year, fe vce(cluster cid)
est store pcre_s
quietly xtreg lco gdpg ren l.fire l.lco lgdp i.year, fe vce(cluster cid)
est store fe_fire_s
quietly xtreg lco gdpg ren l.pfp l.lco lgdp i.year, fe vce(cluster cid)
est store fe_pfp_s
quietly xtreg lco gdpg ren l.pcre l.lco lgdp i.year, fe vce(cluster cid)
est store fe_pcre_s

esttab fire_s pfp_s pcre_s fe_fire_s fe_pfp_s fe_pcre_s, obslast l se noomit ar2 drop(*.year) nonum o(l.fire l.pfp l.pcre gdpg ren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst




/// main model
foreach var in fire pfp pcre{
	quietly gsem (gdpg <- l.`var' l.lco lgdp i.year i.cid) (ren <- l.`var' l.lco lgdp i.year i.cid)(lco <- gdpg ren l.`var' l.lco lgdp i.year i.cid), vce(cluster cid)
	est store sem_`var'_l1
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var']
	nlcom [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[l.`var'] + [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
} 

esttab sem_fire_l1 sem_pfp_l1 sem_pcre_l1, obslast l se noomit ar2 drop(*.year *.cid L.lco) mti("SEM1" "SEM2" "SEM3") nonum o(l.fire l.pfp l.pcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst



/// robustness 1: 2010-2020
preserve
drop if year < 2010

foreach var in fire pfp pcre{
	quiet gsem (gdpg <- l.`var' l.lco lgdp i.year i.cid) (ren <- l.`var' l.lco lgdp i.year i.cid)(lco <- gdpg ren l.`var' l.lco lgdp i.year i.cid), vce(cluster cid)
	est store sem_`var'_l1
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var']
	nlcom [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[l.`var'] + [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
} 

esttab sem_fire_l1 sem_pfp_l1 sem_pcre_l1, obslast l se noomit ar2 drop(*.year *.cid L.lco) mti("SEM1" "SEM2" "SEM3") nonum o(l.fire l.pfp l.pcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst

restore


/// robustness 2: 
preserve
drop if year > 2019
gen period = 5 * floor(year/5)
collapse lco fire pfp pcre ren gdpg fi fm p9950 lman lgdp lpop trd urb lwpop, by (cid period)

xtset cid period

foreach var in fire pfp pcre{
	quiet gsem (gdpg <- L5.`var' L5.lco lgdp i.period i.cid) (ren <- L5.`var' L5.lco lgdp i.period i.cid)(lco <- gdpg ren L5.`var' L5.lco lgdp i.period i.cid), vce(cluster cid)
	est store sem_`var'_l5	
	nlcom [lco]_b[gdpg]*[gdpg]_b[L5.`var']
	nlcom [lco]_b[ren]*[ren]_b[L5.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[L5.`var'] + [lco]_b[ren]*[ren]_b[L5.`var']
	nlcom [lco]_b[L5.`var'] + [lco]_b[gdpg]*[gdpg]_b[L5.`var'] + [lco]_b[ren]*[ren]_b[L5.`var']
} 

esttab sem_fire_l5 sem_pfp_l5 sem_pcre_l5, obslast l se noomit ar2 drop(*.period *.cid L5.lco) mti("SEM1" "SEM2" "SEM3") nonum o(L5.fire L5.pfp L5.pcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst

restore


/// robustness 3: 
foreach var in fi fm{
	quiet gsem (gdpg <- l.`var' l.lco lgdp i.year i.cid) (ren <- l.`var' l.lco lgdp i.year i.cid)(lco <- gdpg ren l.`var' l.lco lgdp i.year i.cid), vce(cluster cid)
	est store sem_`var'_l1
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var']
	nlcom [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[l.`var'] + [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
} 

esttab sem_fi_l1 sem_fm_l1, obslast l se noomit ar2 drop(*.year *.cid L.lco) mti("SEM1" "SEM2") nonum o(l.fi l.fm) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst






*******************************************************************************
/// appendix : fixed effects models with all possible variables
xtreg lco l.fire l.lco p9950 lgdp lman lpop trd urb lwpop i.year, fe vce(cluster cid)
est store fe_fire
xtreg lco l.pfp l.lco p9950 lgdp lman lpop trd urb lwpop i.year, fe vce(cluster cid)
est store fe_pfp
xtreg lco l.pcre l.lco p9950 lgdp lman lpop trd urb lwpop i.year, fe vce(cluster cid)
est store fe_pcre
xtreg lco gdpg ren l.fire l.lco p9950 lgdp lman lpop trd urb lwpop i.year, fe vce(cluster cid)
est store fe_fire_m
xtreg lco gdpg ren l.pfp l.lco p9950 lgdp lman lpop trd urb lwpop i.year, fe vce(cluster cid)
est store fe_pfp_m
xtreg lco gdpg ren l.pcre l.lco p9950 lgdp lman lpop trd urb lwpop i.year, fe vce(cluster cid)
est store fe_pcre_m
esttab fe_fire fe_pfp fe_pcre fe_fire_m fe_pfp_m fe_pcre_m, obslast l se noomit ar2 drop(*.year) nonum o(l.fire l.pfp l.pcre gdpg ren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst

/// main model, all covariates
foreach var in fire pfp pcre{
	quietly gsem (gdpg <- l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid) (ren <- l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid)(lco <- gdpg ren l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid), vce(cluster cid)
	est store sem_`var'_l1
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var']
	nlcom [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[l.`var'] + [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
} 

esttab sem_fire_l1 sem_pfp_l1 sem_pcre_l1, obslast l se noomit ar2 drop(*.year *.cid L.lco) mti("SEM1" "SEM2" "SEM3") nonum o(l.fire l.pfp l.pcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst


///r1
preserve
drop if year < 2010

foreach var in fire pfp pcre{
	quiet gsem (gdpg <- l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid) (ren <- l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid)(lco <- gdpg ren l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid), vce(cluster cid)
	est store sem_`var'_l1
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var']
	nlcom [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[l.`var'] + [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
} 

esttab sem_fire_l1 sem_pfp_l1 sem_pcre_l1, obslast l se noomit ar2 drop(*.year *.cid L.lco) mti("SEM1" "SEM2" "SEM3") nonum o(l.fire l.pfp l.pcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst

restore


///r2
preserve
drop if year > 2019
gen period = 5 * floor(year/5)
collapse lco fire pfp pcre ren gdpg fi fm p9950 lman lgdp lpop trd urb lwpop, by (cid period)

xtset cid period

foreach var in fire pfp pcre{
	quiet gsem (gdpg <- L5.`var' L5.lco lgdp p9950 lman lpop trd urb lwpop i.period i.cid) (ren <- L5.`var' L5.lco lgdp p9950 lman lpop trd urb lwpop i.period i.cid)(lco <- gdpg ren L5.`var' L5.lco lgdp p9950 lman lpop trd urb lwpop i.period i.cid), vce(cluster cid)
	est store sem_`var'_l5	
	nlcom [lco]_b[gdpg]*[gdpg]_b[L5.`var']
	nlcom [lco]_b[ren]*[ren]_b[L5.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[L5.`var'] + [lco]_b[ren]*[ren]_b[L5.`var']
	nlcom [lco]_b[L5.`var'] + [lco]_b[gdpg]*[gdpg]_b[L5.`var'] + [lco]_b[ren]*[ren]_b[L5.`var']
} 

esttab sem_fire_l5 sem_pfp_l5 sem_pcre_l5, obslast l se noomit ar2 drop(*.period *.cid L5.lco) mti("SEM1" "SEM2" "SEM3") nonum o(L5.fire L5.pfp L5.pcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst

restore

///r3
foreach var in fi fm{
	quiet gsem (gdpg <- l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid) (ren <- l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid)(lco <- gdpg ren l.`var' l.lco lgdp p9950 lman lpop trd urb lwpop i.year i.cid), vce(cluster cid)
	est store sem_`var'_l1
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var']
	nlcom [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
	nlcom [lco]_b[l.`var'] + [lco]_b[gdpg]*[gdpg]_b[l.`var'] + [lco]_b[ren]*[ren]_b[l.`var']
} 

esttab sem_fi_l1 sem_fm_l1, obslast l se noomit ar2 drop(*.year *.cid L.lco) mti("SEM1" "SEM2") nonum o(l.fi l.fm) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst
