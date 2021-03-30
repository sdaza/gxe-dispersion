/* this do file reflects the data cleaning and anlaysis as of March 22, 2021 */
clear all
set maxvar 32767


global data /Users/yiyuehuangfu/Desktop/HRS PGS/data
global log /Users/yiyuehuangfu/Desktop/HRS PGS/log
/* please modify the file path for the global "data" to reflect the the position of replication
 package on your devide */ 

/* raw data used to construct the analytic sample:
  - RAND longitudinal HRS, 
  - Polygenic Score file for European and African ancestry, 
  - Global Aging Gateway Harmonized HRS */ 

/************1. prepare data for the analysis ************************************/
use "$data/randhrs1992_2016v2",clear
   rename hhid HHID 
   rename pn PN 
 
 merge 1:1 HHID PN using "$data/PGENSCORE3E_R"  // merge with PGS of European ancestry 
   drop _merge 
  
 merge 1:1 HHID PN using "$data/PGENSCORE3A_R", update  // merge with PGS of African ancestry 
   drop _merge
   
 merge 1:1 HHID PN using "$data/TRK2016TR_R"  // merge with the track file 
   drop _merge

   tempfile randhrs
   save `randhrs' 
   
 use "$data/H_HRS_b",clear /*get father's occupation, financial condition while 
  growing up, and health condition during childhood from Gateway harmonized file */ 
  keep hhid pn radadoccup rafinanch r*healthch  h*rural r*momeft r*momatt r*momgrela ralhchild 
   rename hhid HHID
   rename pn PN 
   merge 1:1 HHID PN using `randhrs'
    drop _merge 
	
   
 /********* 2. generate/recode  variables *********************/
   g bmi_giant15 = EA_PGS3_BMI_GIANT15     
    replace bmi_giant15 = AA_PGS3_BMI_GIANT15  if bmi_giant15==.
	
	sum bmi_giant15  // 15,190 respondents have BMI PGS 
	tab RACE, sum(bmi_giant15) // 12,090 white, 3,100 black 
	
  g t2b_pgs = EA_PGS3_T2D_DIAGRAM12 
   replace t2b_pgs = AA_PGS3_T2D_DIAGRAM12 if t2b_pgs==.
      
	     label var bmi_giant15 "BMI PGS" 
 /*construct childhood condition indicator with  mother's education, father's education
  childhood health, father's occupation when 16, and financial condition while growing up */ 
  
 egen healthch = rowmax(r*healthch)
     label define healthch 1"excellent" 2"very good" 3"above average" 4"fair" 5"poor"
	 label values healthch healthch
	 
 recode healthch (1=5) (2=4) (3=3) (4=2) (5=1), gen(healthch_recode)
recode rafinanch (1=3) (2=2) (3=1), gen(rafinanch_recode)
recode radadoccup (1=2) (2 3=1), gen(fatherocc2)

 recode rafeduc (0/11=1) (12=2) (13/max=3) 
 recode rameduc (0/11=1) (12=2) (13/max=3) 
 recode raeduc (1 =1) (2 3=2) (4 5=3)
  label define edu 1"less than HS" 2"HS" 3">HS" 
  label values raeduc edu 
     
g ec =rafeduc+rameduc+healthch_recode+rafinanch+ fatherocc2 
 label var ec "early life condition" 
 xtile ec_q = ec, n(3)
 tab ec_q, sum(ec)
   label var ec_q "early life condition tertiles" 
   
 label define ec 1"worst" 2"median" 3"best" 
   label values ec_q ec 
  
   
   /* construct a differnt childhood condition indicator with childhood health, financial 
  condition while growing up, mother's effort, mother's attention, good relationship with mother, 
  parental warmth   */ 
 recode r*momeft (1=4) (2=3) (3=2) (4=1)
  egen momeft = rowmax(r9momeft r10momeft r11momeft)
  
 recode r*momatt (1=4) (2=3) (3=2) (4=1)
  egen momatt = rowmax(r*momatt)
  
 egen momgrela = rowmax(r*momgrela)
 
 recode ralhchild (0=3) (1=2) (2=1) (3=0), gen(ralhchild_recode)
 
   mvpatterns rafinanch_recode momeft momatt momgrela ralhchild healthch_recode rafeduc rameduc 

   g ecalt = rafinanch_recode+momeft+momatt+momgrela+healthch_recode
     xtile ecalt_q = ecalt,n(3)
	  tab ecalt_q, sum(ecalt) 
   
  /* cohort dummies  */
g cohort = 1 if BIRTHYR >=1910 & BIRTHYR <= 1924 
replace cohort = 2 if BIRTHYR >= 1925 & BIRTHYR <= 1934
 replace cohort = 3 if BIRTHYR >=1935 & BIRTHYR <=1944
 replace cohort = 4 if BIRTHYR >1945 & BIRTHYR <1960
 
 label define cohort 1"before 1924" 2 "1925-1934" 3"1935-1944" 4"1945-1960",modify 
   label values cohort cohort 
 
   
   
   /* average BMI, both measured and self-reported */ 
  egen srbmi = rowmean(r1bmi-r13bmi)
  egen pmbmi = rowmean(r8pmbmi-r13pmbmi) 
  
   label var srbmi "average self-reported BMI" 
   label var pmbmi "average measured BMI" 

   /*** in order to run the multi level modeling, the data need to be shaped in long
 format */ 
 keep  r*bmi t2b_pgs bmi_giant15  rabplace cohort  r*diab BIRTHYR ///
   RACE GENDER raeduc PC*_* HHID hhidpn  srbmi pmbmi r*pmhghts r*heart h*rural ec* ///
   r*cancr r*shlt r*cenreg r*cendiv *educ rafinanch* healthch* mom*  ralhchild*
   
  drop if bmi_giant15 ==.
  
  rename r?bmi rbmi?
  rename r??bmi rbmi??
  
  rename r?pmbmi rmbmi?
  rename r??pmbmi rmbmi??
  
  rename r?pmhghts rpmhghts? 
  rename r??pmhghts rpmhghts??
  
  rename r?heart rheart?
  rename r??heart rheart??
  rename r?cancr rcancr?
  rename r??cancr rcancr??
  rename r?shlt rshlt?
  rename r??shlt rshlt??
  
  rename r?diab rdiab?
  rename r??diab rdiab??
  
  rename r?cenreg rcenreg? 
  rename r??cenreg rcenreg??
  rename r*cendiv rcendiv* 
  rename h*rural hrural* 
  
  reshape long rmbmi rbmi rpmhghts rheart rcancr rshlt rdiab rcenreg rcendiv hrural, i(hhidpn) j(wave)
   label var rmbmi "measured BMI" 
   label var rbmi "self-reported BMI" 
   
 g age = 1992+ (wave-1)*2 -BIRTHYR 
  g age2 = age*age
  sum age if bmi_giant15 < . & BIRTHYR <1956,d 
   g age_mc = age -66
   g age_mc2 = age_mc * age_mc
    
  drop if rbmi ==. | rbmi==.s | rbmi==.n | rbmi==.m | rbmi==.i | rbmi==.d | rbmi==.x

  bys hhidpn: g n=_n
    label var n "nth interview" 
  bys hhidpn: g N=_N
    label var N "total number of intervew an individual has" 
  
  tab N if n==1 
  
  /* disease */ 
  recode rheart (0 4 5=0) (1 3 6=1)
    label var rheart "report heart disease this wave" 
  recode rcancr (0 4 5 =0) (1 3 6=1)
     label var rcancr "report cancer this wave" 
  recode rdiab (0 4 5 =0) (1 3 =1) 
    label var rdiab "report diabetes this wave" 
	 
	 bys hhidpn: egen heart = max(rheart) 
bys hhidpn: egen cancer = max(rcancr) 
bys hhidpn: egen srh = max(rshlt) 
bys hhidpn: egen diab =max(rdiab)
  label var heart "ever had heart disease"
  label var cancer "ever had cancer" 
  label var srh "worst self-rated health" 
  label var diab "ever had diab" 
  
  /* obesity status */ 
  mark obe if rbmi >=30 
   bys hhidpn: egen obe_any = max(obe)
  mark overweight if rbmi >=25 
   bys hhidpn: egen ow_any = max(overweight) 
   
   //save "$data/working/sample03222021",replace
  
  /*************** descriptives **************************************************/
  /* I ran some descriptives to look at the sample */ 
 sum N if RACE ==1 & BIRTHYR<=1960 & ec_q <.  & n==1
 sum rbmi if RACE ==1 & BIRTHYR<=1960 & ec_q <.
 sum rmbmi if RACE ==1 & BIRTHYR<=1960 & ec_q <.
 tab ec_q if RACE ==1 & BIRTHYR<=1960 & ec_q <. & n==1 
 tab cohort if RACE ==1 & BIRTHYR<=1960 & ec_q <. & n==1 
 tab GENDER if RACE ==1 & BIRTHYR<=1960 & ec_q <.
 tab raeduc if RACE ==1 & BIRTHYR<=1960 & ec_q <.
 
 sum N if RACE ==2 & BIRTHYR<=1960 & ec_q <.  & n==1
 sum rbmi if RACE ==2 & BIRTHYR<=1960 & ec_q <.
 sum rmbmi if RACE ==2 & BIRTHYR<=1960 & ec_q <.
 tab ec_q if RACE ==2 & BIRTHYR<=1960 & ec_q <. & n==1 
 tab cohort if RACE ==2 & BIRTHYR<=1960 & ec_q <. & n==1 
 tab GENDER if RACE ==2 & BIRTHYR<=1960 & ec_q <.
 tab raeduc if RACE ==2 & BIRTHYR<=1960 & ec_q <.
 
   drop if rmbmi ==. 
   
   drop N n 
   bys hhidpn: g n =_n 
   bys hhidpn: g N=_n
   
   sum N if RACE ==1 & BIRTHYR<=1960 & ec_q <.  
   sum N if RACE ==2 & BIRTHYR<=1960 & ec_q <.  
  

  /***************3. replicate Walter (2016) *************************************/
  //log using "$results/models_replicate.log", replace

 // Walter model 4 , white
  mixed rbmi ib1.cohort i.GENDER  ib1.cohort#i.GENDER ///
   bmi_giant15 ib1.cohort#c.bmi_giant15 ///
   age_mc ib1.cohort#c.age_mc i.GENDER#c.age_mc  c.bmi_giant15#c.age_mc ///
   age_mc2 ib1.cohort#c.age_mc2 PC* ///
   ||hhidpn: age_mc age_mc2 if wave <=12 & RACE ==1 & rbmi <46 & BIRTHYR <=1958 , mle 
  outreg2 using "$results/gcm_replicatewalter.xls", append dec(3) alpha(0.001, 0.01, 0.05)
  
 // Walter model 4 , black
  mixed rbmi ib1.cohort i.GENDER  ib1.cohort#i.GENDER ///
   bmi_giant15 ib1.cohort#c.bmi_giant15 ///
   age_mc ib1.cohort#c.age_mc i.GENDER#c.age_mc  c.bmi_giant15#c.age_mc ///
   age_mc2 ib1.cohort#c.age_mc2 PC* ///
   ||hhidpn: age_mc age_mc2 if wave <=12 & RACE ==2  & BIRTHYR <=1958 , mle 
  outreg2 using "$results/gcm_replicatewalter.xls", append dec(3) alpha(0.001, 0.01, 0.05)
  
 // log close 
  
  
 
  /****** 4. models with self report BMI  **********************************************************/
  // log using "$results/models_sr.log", replace
  /* average sr bmi */ 
 reg srbmi  bmi_giant15  ib1.cohort   age_mc age_mc2 i.ec_q  ///
  ib1.cohort#c.bmi_giant15   ///
  i.raeduc i.GENDER PC*  ///
   if n==1 & RACE ==1  ,vce(cluster HHID)
   test [4.cohort#bmi_giant15] = [3.cohort#bmi_giant15]
     est sto OLS 
   outreg2 using "$results/ols_sr_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
   
 
 reg srbmi  bmi_giant15  ib1.cohort     age_mc age_mc2 i.ec_q  ///
  ib1.cohort#c.bmi_giant15 ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==2 & BIRTHYR <=1960  & ec_q <. ,vce(cluster HHID)
      outreg2 using "$results/ols_sr_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	  
   /* average sr bmi, EC X PGS */ 
 reg srbmi  bmi_giant15  ib1.cohort  age_mc age_mc2  ec  ///
  c.ec#c.bmi_giant15  ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==1 & BIRTHYR <=1960  ,vce(cluster HHID)
      outreg2 using "$results/ols_sr_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	
	  
 reg srbmi  bmi_giant15  ib1.cohort    age_mc age_mc2 i.ec2_q  ///
  ib1.ec2_q#c.bmi_giant15 ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==2 & BIRTHYR <=1960  ,vce(cluster HHID)
      outreg2 using "$results/ols_sr_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)

		
  /****************************************************/
 // Walter model 4 +controls , white 
  mixed rbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.cohort#c.bmi_giant15    ib1.cohort#c.age_mc ib1.cohort#c.age_mc2 ///
 c.bmi_giant15#c.age_mc  ///
     i.raeduc  i.GENDER   i.GENDER#c.age_mc i.cohort#i.GENDER   PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==1 & BIRTHYR <=1960 & ec_q <. , mle 
   test [4.cohort#bmi_giant15] = [3.cohort#bmi_giant15]
   
  est sto GC_SR
  
  coefplot OLS GC, keep(*bmi_giant15)  drop(bmi_giant15) label ///
     title("Coefficients on cohort X PGS")
   
  margins , at( bmi_giant15 = (-3 -2 -1 0 1 2 3) cohort = (1 2 3 4) )
   marginsplot 
   
   
   margins , at( bmi_giant15 = (-3 -2 -1 0 1 2 3) rabplace = (1  3  5 6 ) )
      marginsplot 
   
  twoway qfit rbmi age if e(sample) & cohort ==1 || ///    
         qfit rbmi age if e(sample) &   cohort ==2 || ///
		 qfit rbmi age if e(sample) &  cohort ==3 || ///
		 qfit rbmi age if e(sample) &  cohort ==4 , legend(order(1 "before 1924" 2 "1925-1934" 3 "1935-1944" 4 "after 1945") )
   	
  //outreg2 using "$results/gcm_sr_01282021.xls", append dec(3)  alpha(0.001, 0.01, 0.05)
   


 
// Walter model 4 + controls, black
  mixed rbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.cohort#c.bmi_giant15    ib1.cohort#c.age_mc ib1.cohort#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc  PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==2 & BIRTHYR <=1960   , mle 
   outreg2 using "$results/gcm_sr_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	

	
  /****************************************************/
  // Walter model 4 + control + EC X PGS, white 
mixed rbmi  bmi_giant15 ib1.cohort age_mc  age_mc2 i.ec2_q  ///
   ib1.ec_q#c.bmi_giant15    ib1.ec_q#c.age_mc ib1.ec_q#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc i.rheart i.rcancr rshlt  PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==1 & BIRTHYR <=1960  & ec_q <. , mle 
   est sto EC_SR
  outreg2 using "$results/gcm_sr_01282021.xls", append dec(3)  alpha(0.001, 0.01, 0.05)
 
// Walter model 4 + controls + EC X PGS , black
mixed rbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.ec_q#c.bmi_giant15    ib1.ec_q#c.age_mc ib1.ec_q#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc  i.rheart i.rcancr rshlt PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==2 & BIRTHYR <=1960  & ec_q <. , mle  
   outreg2 using "$results/gcm_sr_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	
//log close
  /****************************************************/
  
   
  /**************************************************************
  4.1 counterfactual BMI if cohort born after 1945 were born 1935-1944
  *****************************************************************/
   
    mixed rbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.cohort#c.bmi_giant15    ib1.cohort#c.age_mc ib1.cohort#c.age_mc2 ///
 c.bmi_giant15#c.age_mc  ///
     i.raeduc  i.GENDER   i.GENDER#c.age_mc i.cohort#i.GENDER   PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==1 & BIRTHYR <=1960 & ec_q <. , mle 
   
   
   // model predicted BMI 
   predict pdct_bmi if e(sample) & cohort ==4 
    label var pdct_bmi "model predicted BMI" 
   
   drop cohort2 
	 tab cohort, gen(cohort)
   tab ec_q, gen(ec_q)
   tab raeduc, gen(raeduc)
   tab GENDER, gen(GENDER)
    recode bmi_giant15 (min/-2=1) (-2/-1=2) (-1/1=3) (1/2=4) (2/5=5), gen(pgs_cat)
	
   // if cohort 4 were born in cohort 3 
  g cf1_bmi = _b[bmi_giant15]*bmi_giant15 + ///
   _b[3.cohort]*cohort4 + ///
   _b[age_mc]*age_mc + _b[age_mc2]*age_mc2 + _b[1.ec_q]*ec_q1 + _b[2.ec_q]*ec_q2 + _b[3.ec_q]*ec_q3 + ///
   _b[3.cohort#c.bmi_giant15]*cohort4*bmi_giant15 + ///
   _b[3.cohort#age_mc]*cohort4*age_mc + _b[3.cohort#age_mc2]*cohort4*age_mc2 + ///
   _b[bmi_giant15#age_mc]*bmi_giant15*age_mc + ///
   _b[1.raeduc]*raeduc1 + _b[2.raeduc]*raeduc2 +_b[3.raeduc]*raeduc3 + ///
   _b[1.GENDER]*GENDER1 + _b[2.GENDER]*GENDER2 + ///
   _b[3.cohort#2.GENDER]*cohort4*GENDER2 + ///
   _b[2.GENDER#age_mc]*GENDER2*age_mc + ///
   _b[PC1_5A]*PC1_5A + _b[PC1_5B]*PC1_5B + _b[PC1_5C]*PC1_5C + _b[PC1_5D]*PC1_5D  + _b[PC1_5E]*PC1_5E + ///
   _b[PC6_10A]*PC6_10A + _b[PC6_10B]*PC6_10B + _b[PC6_10C]*PC6_10C + _b[PC6_10D]*PC6_10D + _b[PC6_10E]*PC6_10E + ///
   26.43285  if cohort ==4 & e(sample) ==1 
   
   label var cf1_bmi "if cohort4 were born in 1935-1944" 
   
   // if cohort 4 had the impact of PGS as cohort 3 
  g cf2_bmi = _b[bmi_giant15]*bmi_giant15 + ///
   _b[4.cohort]*cohort4 + ///
   _b[age_mc]*age_mc + _b[age_mc2]*age_mc2 + _b[1.ec_q]*ec_q1 + _b[2.ec_q]*ec_q2 + _b[3.ec_q]*ec_q3 + ///
   _b[3.cohort#c.bmi_giant15]*cohort4*bmi_giant15 + ///
   _b[4.cohort#age_mc]*cohort4*age_mc + _b[4.cohort#age_mc2]*cohort4*age_mc2 + ///
   _b[bmi_giant15#age_mc]*bmi_giant15*age_mc + ///
   _b[1.raeduc]*raeduc1 + _b[2.raeduc]*raeduc2 +_b[3.raeduc]*raeduc3 + ///
   _b[1.GENDER]*GENDER1 + _b[2.GENDER]*GENDER2 + ///
   _b[4.cohort#2.GENDER]*cohort4*GENDER2 + ///
   _b[2.GENDER#age_mc]*GENDER2*age_mc + ///
   _b[PC1_5A]*PC1_5A + _b[PC1_5B]*PC1_5B + _b[PC1_5C]*PC1_5C + _b[PC1_5D]*PC1_5D  + _b[PC1_5E]*PC1_5E + ///
   _b[PC6_10A]*PC6_10A + _b[PC6_10B]*PC6_10B + _b[PC6_10C]*PC6_10C + _b[PC6_10D]*PC6_10D + _b[PC6_10E]*PC6_10E + ///
   26.43285  if cohort ==4 & e(sample) ==1 
   
   label var cf2_bmi "if cohort4 has the coefficient on PGS as 1935-1944" 
	

   label define gender 1"male" 2"female" ,modify
    label values GENDER  gender
   
    recode cf1_bmi (min/25 = 0) (25/max =1), gen(cf1_ow) 
	recode cf1_bmi (min/30=0) (30/max=1), gen(cf1_obe)
	recode cf2_bmi (min/25=0) (25/max=1), gen(cf2_ow)
	recode cf2_bmi (min/30=0) (30/max=1), gen(cf2_obe)
	recode pdct_bmi (min/25=0) (25/max=1), gen(pdct_ow)
	recode pdct_bmi (min/30=0) (30/max=1), gen(pdct_obe)
		  
		  
		 tab GENDER cf1_ow if RACE==1 & cohort ==4,row
		 tab GENDER cf2_ow if RACE==1 & cohort==4, row
		 tab GENDER pdct_ow if RACE ==1 & cohort==4, row
		 
		 tab GENDER cf1_obe if RACE==1 & cohort ==4,row
		 tab GENDER cf2_obe  if RACE==1 & cohort==4, row
		 tab GENDER pdct_obe if RACE ==1 & cohort==4, row

  
   // if cohort 4 were born in cohort 1
  g cf3_bmi = _b[bmi_giant15]*bmi_giant15 + ///
   _b[1.cohort]*cohort4 + ///
   _b[age_mc]*age_mc + _b[age_mc2]*age_mc2 + _b[1.ec_q]*ec_q1 + _b[2.ec_q]*ec_q2 + _b[3.ec_q]*ec_q3 + ///
   _b[1.cohort#c.bmi_giant15]*cohort4*bmi_giant15 + ///
   _b[1.cohort#age_mc]*cohort4*age_mc + _b[1.cohort#age_mc2]*cohort4*age_mc2 + ///
   _b[bmi_giant15#age_mc]*bmi_giant15*age_mc + ///
   _b[1.raeduc]*raeduc1 + _b[2.raeduc]*raeduc2 +_b[3.raeduc]*raeduc3 + ///
   _b[1.GENDER]*GENDER1 + _b[2.GENDER]*GENDER2 + ///
   _b[1.cohort#2.GENDER]*cohort4*GENDER2 + ///
   _b[2.GENDER#age_mc]*GENDER2*age_mc + ///
   _b[PC1_5A]*PC1_5A + _b[PC1_5B]*PC1_5B + _b[PC1_5C]*PC1_5C + _b[PC1_5D]*PC1_5D  + _b[PC1_5E]*PC1_5E + ///
   _b[PC6_10A]*PC6_10A + _b[PC6_10B]*PC6_10B + _b[PC6_10C]*PC6_10C + _b[PC6_10D]*PC6_10D + _b[PC6_10E]*PC6_10E + ///
   26.43285  if cohort ==4 & e(sample) ==1 
   
   label var cf3_bmi "if cohort4 were born in 1900-1924" 
   
   // if cohort 4 had the impact of PGS as cohort 1 
  g cf4_bmi = _b[bmi_giant15]*bmi_giant15 + ///
   _b[4.cohort]*cohort4 + ///
   _b[age_mc]*age_mc + _b[age_mc2]*age_mc2 + _b[1.ec_q]*ec_q1 + _b[2.ec_q]*ec_q2 + _b[3.ec_q]*ec_q3 + ///
   _b[1.cohort#c.bmi_giant15]*cohort4*bmi_giant15 + ///
   _b[4.cohort#age_mc]*cohort4*age_mc + _b[4.cohort#age_mc2]*cohort4*age_mc2 + ///
   _b[bmi_giant15#age_mc]*bmi_giant15*age_mc + ///
   _b[1.raeduc]*raeduc1 + _b[2.raeduc]*raeduc2 +_b[3.raeduc]*raeduc3 + ///
   _b[1.GENDER]*GENDER1 + _b[2.GENDER]*GENDER2 + ///
   _b[4.cohort#2.GENDER]*cohort4*GENDER2 + ///
   _b[2.GENDER#age_mc]*GENDER2*age_mc + ///
   _b[PC1_5A]*PC1_5A + _b[PC1_5B]*PC1_5B + _b[PC1_5C]*PC1_5C + _b[PC1_5D]*PC1_5D  + _b[PC1_5E]*PC1_5E + ///
   _b[PC6_10A]*PC6_10A + _b[PC6_10B]*PC6_10B + _b[PC6_10C]*PC6_10C + _b[PC6_10D]*PC6_10D + _b[PC6_10E]*PC6_10E + ///
   26.43285  if cohort ==4 & e(sample) ==1 
   
   label var cf4_bmi "if cohort4 has the coefficient on PGS as 1900-1924" 
  
     
/* 4.2 extrapolation *********************************************/
  // assuming the main effect for cohort 5 is 4.54, the interaction effect is .55
	 
 g cohort5 = 1 if BIRTHYR >1960 & BIRTHYR <=1970
 
 g extra_bmi = 	 _b[bmi_giant15]*bmi_giant15 + ///
   4.54*cohort5 + ///
   _b[age_mc]*age_mc + _b[age_mc2]*age_mc2 + _b[1.ec_q]*ec_q1 + _b[2.ec_q]*ec_q2 + _b[3.ec_q]*ec_q3 + ///
   .55*cohort5*bmi_giant15 + ///
   -.14*cohort5*age_mc +  ///
   _b[bmi_giant15#age_mc]*bmi_giant15*age_mc + ///
   _b[1.raeduc]*raeduc1 + _b[2.raeduc]*raeduc2 +_b[3.raeduc]*raeduc3 + ///
   _b[1.GENDER]*GENDER1 + _b[2.GENDER]*GENDER2 + ///
   .758*cohort5*GENDER2 + ///
   _b[2.GENDER#age_mc]*GENDER2*age_mc + ///
   _b[PC1_5A]*PC1_5A + _b[PC1_5B]*PC1_5B + _b[PC1_5C]*PC1_5C + _b[PC1_5D]*PC1_5D  + _b[PC1_5E]*PC1_5E + ///
   _b[PC6_10A]*PC6_10A + _b[PC6_10B]*PC6_10B + _b[PC6_10C]*PC6_10C + _b[PC6_10D]*PC6_10D + _b[PC6_10E]*PC6_10E + ///
   26.43285  
   
   sum extra_bmi if cohort5 ==1 & ec_q<. & RACE ==1 & GENDER ==1 
   
   sort age RACE GENDER
      bys age RACE GENDER: egen extra_bmimean = mean(extra_bmi) if cohort5==1
   bys age RACE GENDER: egen real_bmimean = mean(rbmi) if  cohort5==1

   
   twoway line real_bmimean age if RACE==1  || ///
          line extra_bmimean age if RACE ==1   ||, ///
		  by(GENDER, title("HRS respondents born 1960-70, white") ) ///
		  legend(order(1 "real BMI" 2 "extrapolated BMI" ))  
		 
 
	 
/**** 5. do everything with measured BMI ************************/
 log using "$results/models_pm.log", replace
  /* average sr bmi */ 
 reg pmbmi  bmi_giant15  ib1.cohort    age_mc age_mc2 i.ec_q  ///
  ib1.cohort#c.bmi_giant15 ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==1 & BIRTHYR <=1960  & ec_q <. ,vce(cluster HHID)
   outreg2 using "$results/ols_pm_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
   
 reg pmbmi  bmi_giant15  ib1.cohort    age_mc age_mc2 i.ec_q  ///
  ib1.cohort#c.bmi_giant15 ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==2 & BIRTHYR <=1960  & ec_q <. ,vce(cluster HHID)
      outreg2 using "$results/ols_pm_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	  
   /* average sr bmi, EC X PGS */ 
 reg pmbmi  bmi_giant15  ib1.cohort    age_mc age_mc2 i.ec_q  ///
  ib1.ec_q#c.bmi_giant15 ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==1 & BIRTHYR <=1960  & ec_q <. ,vce(cluster HHID)
      outreg2 using "$results/ols_pm_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	  
 reg pmbmi  bmi_giant15  ib1.cohort    age_mc age_mc2 i.ec_q  ///
  ib1.ec_q#c.bmi_giant15 ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==2 & BIRTHYR <=1960  & ec_q <. ,vce(cluster HHID)
      outreg2 using "$results/ols_pm_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	  
  
  /****************************************************/
 // Walter model 4 +controls , white 
  mixed rmbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.cohort#c.bmi_giant15    ib1.cohort#c.age_mc ib1.cohort#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc  i.rheart i.rcancr rshlt i.rpmhghts PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==1 & BIRTHYR <=1960  & ec_q <. , mle 
   est sto GC_M
   
   coefplot GC_SR , keep(*bmi_giant15 1.*)  label ///
    title("Coefficients on cohort X PGS") plotlabels("Self-reported BMI")
   
   coefplot GC_SR GC_M, keep(*bmi_giant15 cohort)  label ///
    title("Coefficients on cohort X PGS") ///
	plotlabels("Self-reported BMI" "Measured BMI")
   
   twoway qfit rmbmi age if e(sample) & cohort ==1  || ///
          qfit rmbmi age if e(sample) & cohort ==2 || ///
		  qfit rmbmi age if e(sample) & cohort ==3 || ///
		  qfit rmbmi age if e(sample) & cohort ==4 , legen(order (1 "before 1924" 2 "1925-34"  3 "1935-44" 4 "1945-60"))
   
   margins , at( bmi_giant15 = (-3 -2 -1 0 1 2 3) cohort = (1 2 3 4) )
   marginsplot 
   
   // Walter model 4 +controls , white , linear growth curve only 
  mixed rmbmi  bmi_giant15 ib1.cohort age_mc   i.ec_q  ///
   ib1.cohort#c.bmi_giant15    ib1.cohort#c.age_mc ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc i.rheart i.rcancr rshlt i.rpmhghts  PC* ///
   ||hhidpn: age_mc  if RACE ==1 & BIRTHYR <=1960  & ec_q <. , mle
   
     margins , at( bmi_giant15 = (-4 -2 0 2 4) cohort = (1 2 3 4) )
   marginsplot 
   
  outreg2 using "$results/gcm_pm_01282021.xls", append dec(3)  alpha(0.001, 0.01, 0.05)
 
// Walter model 4 + controls, black
  mixed rmbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.cohort#c.bmi_giant15    ib1.cohort#c.age_mc ib1.cohort#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc i.rheart i.rcancr rshlt i.rpmhghts  PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==2 & BIRTHYR <=1960  & ec_q <. , mle 
   
    
   outreg2 using "$results/gcm_pm_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	

	
  /****************************************************/
  // Walter model 4 + control + EC X PGS, white 
mixed rmbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.ec_q#c.bmi_giant15    ib1.ec_q#c.age_mc ib1.ec_q#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc i.rheart i.rcancr rshlt i.rpmhghts   PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==1 & BIRTHYR <=1960  & ec_q <. , mle 
   est sto EC_M
   
   coefplot EC_SR , keep(*bmi_giant15) drop(bmi_giant15) label ///
    title("Coefficients on early condition X PGS") 
   
   
   coefplot EC_SR EC_M, keep(*bmi_giant15) drop(bmi_giant15) label ///
    title("Coefficients on early condition X PGS") ///
	plotlabels("Self-reported BMI" "Measured BMI")
   
    /*mark hasmbmi if e(sample) ==1 
  drop if hasmbmi ==0 
   bys hhidpn: g n2= _n 
   bys hhidpn: g N2 = _N
   
     tab N2 if n2==1  */ 
   
  outreg2 using "$results/gcm_pm_01282021.xls", append dec(3)  alpha(0.001, 0.01, 0.05)
 
// Walter model 4 + controls + EC X PGS , black
mixed rmbmi  bmi_giant15 ib1.cohort age_mc   age_mc2 i.ec_q  ///
   ib1.ec_q#c.bmi_giant15    ib1.ec_q#c.age_mc ib1.ec_q#c.age_mc2 ///
 c.bmi_giant15#c.age_mc ///
     i.raeduc  i.GENDER  ib1.cohort#i.GENDER  i.GENDER#c.age_mc  i.rheart i.rcancr rshlt i.rpmhghts  PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==2 & BIRTHYR <=1960  & ec_q <. , mle 
   

   outreg2 using "$results/gcm_pm_01282021.xls", append dec(3) alpha(0.001, 0.01, 0.05)
	
log close


/*************** depression cohort *********************************************/

mixed rbmi  bmi_giant15 i.cohort2 age_mc   age_mc2  ///
   i.cohort2#c.bmi_giant15    i.cohort2#c.age_mc i.cohort2#c.age_mc2 ///
 c.bmi_giant15#c.age_mc  ///
     i.raeduc  i.GENDER   i.GENDER#c.age_mc i.cohort2#i.GENDER  PC* ///
   ||hhidpn: age_mc age_mc2 if RACE ==1 & BIRTHYR <=1960 & ec_q <. , mle 
   
bys hhidpn: egen diab = max(rdiab)

 
reg diab t2b_pgs i.cohort2  ///
   i.cohort2#c.t2b_pgs    i.raeduc  i.GENDER   i.cohort2#i.GENDER  PC* if n ==1 & RACE==1 
   
reg diab t2b_pgs i.cohort2 BIRTHYR   ///
   i.cohort2#c.t2b_pgs    i.raeduc  i.GENDER   i.cohort2#i.GENDER  PC* if n ==1 & RACE==2 
  
/************* instrumental variables ? *******/ 
// first stage 
oprobit ec_q i.rabplace i.rcendiv bmi_giant15  ib1.cohort  age_mc age_mc2 i.heart i.cancer i.srh   ///
  i.raeduc i.GENDER PC* ///
   if n==1 & RACE ==1 & BIRTHYR <=1960  

   predict ec_q1p ec_q2p ec_q3p 
   
   tab ec_q, gen(ec_q)

ivregress 2sls srbmi  bmi_giant15    ib1.ec_q#c.bmi_giant15  /// 
ib1.cohort  age_mc age_mc2 i.heart i.cancer i.srh  ///
i.raeduc i.GENDER PC* (ec_q2 ec_q3 =ec_q2p ec_q3p) ///
   if n==1 & RACE ==1 & BIRTHYR <=1960  & ec_q <. ,vce(cluster HHID)
   
   
