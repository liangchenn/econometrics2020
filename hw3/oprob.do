//question 1


gen N = airlineAA +  airlineDL + airlineUA + airlineAL + airlineLCC  + airlineWN 
gen lnN = log(N)
replace lnN = 0 if lnN == .

eststo clear
eststo: oprobit N marketdistance marketsize lnN
