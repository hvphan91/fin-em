/// prep data

// cd...
// capture log close 
// log using [log_file], replace text

import excel fin-em-clean, firstrow case(lower)

drop if inc == 2
drop if year < 1995 | year > 2020

gen p9950 = p99/p50*100
gen p9050 = p90/p50*100
gen man = manva*gdp
gen gdpc = gdp/pop
gen coc = co/pop
gen wp = wpop/pop*100

for var p99 p90 p50 pfp fi fm: replace X = X*100
for var co coc fire pfp pcre ren indem manva p9950 gdpc pop wp urb trd fi fm: gen lX = ln(X+0.01)

encode country, gen(cid)
xtset cid year

// trends
preserve
collapse fire pcre co gdpc ren, by(year)
graph twoway scatter fire year, msymbol(O) connect(l) yaxis(1) || scatter pcre year, msymbol(S) connect(l) yaxis(2)
graph twoway scatter gdpc year, msymbol(O) connect(l) yaxis(1) || scatter ren year, msymbol(S) connect(l) yaxis(2) 
graph twoway scatter co year, msymbol(O) connect(l) yaxis(1) 
restore


// fixed effects with and without mediator
foreach var in lfire lpcre{
	quiet xtreg lco `var' lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtnm_`var'_l
	quiet xtreg lco `var' lgdpc lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtm_`var'_l
	quiet xtreg lco c.`var'##c.`var' lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtnm_`var'_nl
	quiet xtreg lco c.`var'##c.`var' lgdpc lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtm_`var'_nl
}

esttab xtnm_lfire_l xtnm_lfire_nl xtnm_lpcre_l xtnm_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)
esttab xtm_lfire_l  xtm_lfire_nl  xtm_lpcre_l  xtm_lpcre_nl,  obslast l b(3) se(3) noomit ar2 drop(*.year) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) 


// graphs pcre squared models
quiet reg lco c.lpcre##c.lpcre lp9950 lurb ltrd lpop l.lco i.year i.cid, vce(cluster cid)
quiet margins, at(lpcre = (1.9(0.05)5.7)) saving(m1)
quiet reg lco c.lpcre##c.lpcre lgdpc lren lp9950 lurb ltrd lpop l.lco i.year i.cid, vce(cluster cid)
quiet margins, at(lpcre = (1.9(0.05)5.7)) saving(m2)

combomarginsplot m1 m2, recast(line) recastci(rarea)


// main mediation models
foreach var in lfire lpcre{
	quiet gsem (lgdpc <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- c.`var'##c.`var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_nl
	
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var'] + [lco]_b[c.`var'#c.`var']
	
	quiet gsem (lgdpc <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- `var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[`var']
}

esttab gsem_lfire_l gsem_lfire_nl gsem_lpcre_l gsem_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)


// r1: 5-year avg
preserve
drop if year > 2019
gen period = 5 * floor(year/5)
collapse lco lfire lpfp lpcre lren lgdpc lp9950 lpop ltrd lurb, by (cid period)

xtset cid period

foreach var in lfire lpcre{
	quietly gsem (lgdpc <- L5.`var' lp9950 ltrd lurb lpop L5.lco i.period i.cid) (lren <- L5.`var' lp9950 ltrd lurb lpop L5.lco i.period i.cid) (lco <- L5.`var' lgdpc lren lp9950 ltrd lurb lpop L5.lco i.period i.cid), vce(cluster cid)
	est store sem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[L5.`var']
	
	quietly gsem (lgdpc <- cL5.`var'##cL5.`var' lp9950 ltrd lurb lpop L5.lco i.period i.cid) (lren <- cL5.`var'##cL5.`var' lp9950 ltrd lurb lpop L5.lco i.period i.cid) (lco <- cL5.`var'##cL5.`var' lgdpc lren lp9950 ltrd lurb lpop L5.lco i.period i.cid), vce(cluster cid)
	est store sem_`var'_nl
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[L5.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[cL5.`var'#cL5.`var']
	
} 
esttab sem_lfire_l sem_lfire_nl sem_lpcre_l sem_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.period *.cid) nonum o(L5.lfire cL5.lfire#cL5.lfire L5.lpcre cL5.lpcre#cL5.lpcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) noconst
restore


// r2: alternative fin var.
foreach var in lpfp lfi lfm {
	quiet gsem (lgdpc <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- c.`var'##c.`var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_nl
	
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var'] + [lco]_b[c.`var'#c.`var']
	
	quiet gsem (lgdpc <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- `var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[`var']
}

esttab gsem_lpfp_l gsem_lpfp_nl gsem_lfi_l gsem_lfi_nl gsem_lfm_l gsem_lfm_nl, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lpfp c.lpfp#c.lpfp lfi c.lfi#c.lfi lfm c.lfm#c.lfm lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)


// T.B2
gsem (lgdpc <- c.lpcre##i.year lp9950 ltrd lurb lpop l.lco i.cid) (lren <- c.lpcre##i.year lp9950 ltrd lurb lpop l.lco i.cid) (lco <- c.lpcre##i.year lgdpc lren lp9950 ltrd lurb lpop l.lco i.cid), vce(cluster cid)
est store gsem_lpcre

forvalues num = 1996/2019 {
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.lpcre#`num'.year]
}

esttab gsem_lpcre, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lfire lpcre c.lpcre#c.lpcre lgdpc lren)


// T.B4
foreach var in lfire lpcre{
	quiet gsem (lgdpc <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lgdpc l.lren l.lco i.year i.cid) (lren <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lgdpc l.lren l.lco i.year i.cid) (lco <- c.`var'##c.`var' lgdpc lren lp9950 ltrd lurb lpop l.lgdpc l.lren l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_nl
	
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var'] + [lco]_b[c.`var'#c.`var']
	
	quiet gsem (lgdpc <- `var' lp9950 ltrd lurb lpop l.lgdpc l.lren l.lco i.year i.cid) (lren <- `var' lp9950 ltrd lurb lpop l.lgdpc l.lren l.lco i.year i.cid) (lco <- `var' lgdpc lren lp9950 ltrd lurb lpop l.lgdpc l.lren l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[`var']
}

esttab gsem_lfire_l gsem_lpcre_l gsem_lfire_nl gsem_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)


// T.C1-2
preserve
drop if year > 2019

foreach var in lfire lpcre{
	quiet xtreg lco `var' lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtnm_`var'_l
	quiet xtreg lco `var' lgdpc lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtm_`var'_l
	quiet xtreg lco c.`var'##c.`var' lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtnm_`var'_nl
	quiet xtreg lco c.`var'##c.`var' lgdpc lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store xtm_`var'_nl
}

esttab xtnm_lfire_l xtnm_lfire_nl xtnm_lpcre_l xtnm_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)
esttab xtm_lfire_l  xtm_lfire_nl  xtm_lpcre_l  xtm_lpcre_nl,  obslast l b(3) se(3) noomit ar2 drop(*.year) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001) 


foreach var in lfire lpcre{
	quiet gsem (lgdpc <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- c.`var'##c.`var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_nl
	
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var'] + [lco]_b[lren]*[lren]_b[c.`var'#c.`var'] + [lco]_b[c.`var'#c.`var']
	
	quiet gsem (lgdpc <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- `var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lren]*[lren]_b[`var'] + [lco]_b[lgdpc]*[lgdpc]_b[`var'] + [lco]_b[`var']
}

esttab gsem_lfire_l gsem_lfire_nl gsem_lpcre_l gsem_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)
restore


// T.D1-2
foreach var in lfire lpcre{
	quiet xtreg lco `var' c.lgdpc##c.lgdpc c.lren##c.lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store `var'_l
	quiet xtreg lco c.`var'##c.`var' c.lgdpc##c.lgdpc c.lren##c.lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store `var'_nl
	quiet xtreg lco `var' lgdpc c.lren##c.lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store `var'_ls
	quiet xtreg lco c.`var'##c.`var' lgdpc c.lren##c.lren lp9950 lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store `var'_nls
}
esttab lfire_l lfire_nl lpcre_l lpcre_nl lfire_ls lfire_nls lpcre_ls lpcre_nls, obslast l b(3) se(3) noomit ar2 drop(*.year) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc c.lgdpc#c.lgdpc lren c.lren#c.lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)


foreach var in lfire lpcre{
	quiet gsem (lgdpc <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- c.`var'##c.`var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- c.`var'##c.`var' lgdpc c.lren##c.lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_nl
	
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var']

	quiet gsem (lgdpc <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lren <- `var' lp9950 ltrd lurb lpop l.lco i.year i.cid) (lco <- `var' lgdpc c.lren##c.lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
}

esttab gsem_lfire_l gsem_lfire_nl gsem_lpcre_l gsem_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)


// T.E1
foreach var in lfire lpcre{
	quiet xtreg lp9950 `var' lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store m_`var'_1
	quiet xtreg lp9950 `var' lgdpc lren lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store m_`var'_2
	quiet xtreg lp9950 c.`var'##c.`var' lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store m_`var'_3
	quiet xtreg lp9950 c.`var'##c.`var' lgdpc lren lurb ltrd lpop l.lco i.year, fe vce(cluster cid)
	est store m_`var'_4
}

esttab m_lfire_1 m_lfire_2 m_lfire_3 m_lfire_4 m_lpcre_1 m_lpcre_2 m_lpcre_3 m_lpcre_4, obslast l b(3) se(3) noomit ar2 drop(*.year) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)

foreach var in lfire lpcre{
	quiet gsem (lgdpc <- c.`var'##c.`var' ltrd lurb lpop l.lco i.year i.cid) (lp9950 <- c.`var'##c.`var' ltrd lurb lpop l.lco i.year i.cid) (lren <- c.`var'##c.`var' ltrd lurb lpop l.lco i.year i.cid) (lco <- c.`var'##c.`var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_nl
	
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[c.`var'#c.`var']
	nlcom [lco]_b[lp9950]*[lp9950]_b[`var']
	nlcom [lco]_b[lp9950]*[lp9950]_b[c.`var'#c.`var']

	quiet gsem (lgdpc <- `var' ltrd lurb lpop l.lco i.year i.cid) (lp9950 <- `var' ltrd lurb lpop l.lco i.year i.cid) (lren <- `var' ltrd lurb lpop l.lco i.year i.cid) (lco <- `var' lgdpc lren lp9950 ltrd lurb lpop l.lco i.year i.cid), vce(cluster cid)
	est store gsem_`var'_l
	nlcom [lco]_b[lgdpc]*[lgdpc]_b[`var']
	nlcom [lco]_b[lp9950]*[lp9950]_b[`var']
}

esttab gsem_lfire_l gsem_lfire_nl gsem_lpcre_l gsem_lpcre_nl, obslast l b(3) se(3) noomit ar2 drop(*.year *.cid) nonum o(lfire c.lfire#c.lfire lpcre c.lpcre#c.lpcre lgdpc lren lp9950) star(+ 0.1 * 0.05 ** 0.01 *** 0.001)