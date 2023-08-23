;; MIBCOM - Metabolic Individual-Based COMmunity model: Dynamic, allometric, metabolic, individual-based, spatial community model for cental-place forager mammals in heterogenous landscapes
;; copyright Leonna Szangolies, Florian Jeltsch, Cara Gallagher, University of Potsdam

extensions [profiler]

globals [
  ;;species definition with mass, shelter need, foraging type and food preference
  s0mass s0massdev s1mass s1massdev  s2mass s2massdev s3mass s3massdev s4mass s4massdev s5mass s5massdev s6mass s6massdev s7mass s7massdev s8mass s8massdev s9mass s9massdev
  s0shelt s1shelt s2shelt s3shelt s4shelt s5shelt s6shelt s7shelt s8shelt s9shelt
  s0food s1food s2food s3food s4food s5food s6food s7food s8food s9food
  s0fortype s1fortype s2fortype s3fortype s4fortype s5fortype s6fortype s7fortype s8fortype s9fortype

  ;;landscape variables
  feed-matrix feed-struct season-var ;;feed-var
  actcov cover_large cover_small

  ;;helper variables for the functions
  day r-search feed shelter-need maxrad minfeed spec-color movecost foodshare

  ;;model parameters (need sensitivity test)
  spec_num mass_order n_immi max_nat_disp growthcost patch-move cover var bold_prob

  ;; scenario definition
  density_dep small_habitat disturb high-season

  ;;output variables
  focal_spec                                 ;;mainly for plotting home ranges of only one species
  fail                                       ;;juveniles that did not find a home range (per time step)
  fail3                                      ;;individuals that dy because of no storage
  failed-lact                                ;;offspring that gets lost during lactation due to shortage
  repro                                      ;;weaned juveniles that can search their own home range
  suc_juv                                    ;;successful number of juvenile over whole time per species
  suc_immi
  rep_success_0 rep_success_1 rep_success_2 rep_success_3 rep_success_4 rep_success_5 rep_success_6 rep_success_7 rep_success_8 rep_success_9  ;;total weaned juveniles per female
  rep_success_hr_0 rep_success_hr_1 rep_success_hr_2 rep_success_hr_3 rep_success_hr_4 rep_success_hr_5 rep_success_hr_6 rep_success_hr_7 rep_success_hr_8 rep_success_hr_9  ;;total weaned juveniles per female
  life_juv0 life_juv1 life_juv2 life_juv3 life_juv4 life_juv5 life_juv6 life_juv7 life_juv8 life_juv9                                          ;;mean rep_success
  life_juv_hr0 life_juv_hr1 life_juv_hr2 life_juv_hr3 life_juv_hr4 life_juv_hr5 life_juv_hr6 life_juv_hr7 life_juv_hr8 life_juv_hr9                                          ;;mean rep_success
  patches0 patches1 patches2 patches3 patches4 patches5 patches6 patches7 patches8 patches9
  compet0 compet1 compet2 compet3 compet4 compet5 compet6 compet7 compet8 compet9
  order0 order1 order2 order3 order4 order5 order6 order7 order8 order9
  vis_before0 vis_before1 vis_before2 vis_before3 vis_before4 vis_before5 vis_before6 vis_before7 vis_before8 vis_before9
  number number0 number1 number2 number3 number4 number5 number6 number7 number8 number9                                                              ;;number of individuals
  hr hr0 hr1 hr2 hr3 hr4 hr5 hr6 hr7 hr8 hr9                                                                                                   ;;mean home range size
  maxhrs maxhr0 maxhr1 maxhr2 maxhr3 maxhr4 maxhr5 maxhr6 maxhr7 maxhr8 maxhr9
  countmaxhr countmaxhr0 countmaxhr1 countmaxhr2 countmaxhr3 countmaxhr4 countmaxhr5 countmaxhr6 countmaxhr7 countmaxhr8 countmaxhr9
  ages age0 age1 age2 age3 age4 age5 age6 age7 age8 age9                                                                                                   ;;mean home range size
  stors stor0 stor1 stor2 stor3 stor4 stor5 stor6 stor7 stor8 stor9                                                                            ;;mean storage
  abs_stor0 abs_stor1 abs_stor2 abs_stor3 abs_stor4 abs_stor5 abs_stor6 abs_stor7 abs_stor8 abs_stor9                                                                            ;;mean absolute storage of non-pregnant
  rm rm0 rm1 rm2 rm3 rm4 rm5 rm6 rm7 rm8 rm9                                                                                                   ;;mean real mass
  hrpreg hrpreg0 hrpreg1 hrpreg2 hrpreg3 hrpreg4 hrpreg5 hrpreg6 hrpreg7 hrpreg8 hrpreg9                                                       ;;mean home range size during pregnancy and lactation
  storpreg storpreg0 storpreg1 storpreg2 storpreg3 storpreg4 storpreg5 storpreg6 storpreg7 storpreg8 storpreg9 stor_pregstart stor_pregend     ;;mean storage during pregnancy and lactation
  comp5 comp10 comp15 comp20 comp25 compmore                                                                                                   ;;mean number of individuals with 5 or less / 5 - 10 / 10 - 15 / 15 - 20 / 20 - 25 / 25 or more competitors per food patch
  comp5_stor comp5_fmr comp10_stor comp10_fmr comp15_stor comp15_fmr comp20_stor comp20_fmr comp25_stor comp25_fmr compmore_stor compmore_fmr  ;;mean storage and field metabolic rate of turtles with respective competition
  abs_juv juv juv0 juv1 juv2 juv3 juv4 juv5 juv6 juv7 juv8 juv9                                                                                        ;;mean percent of weaned juveniles vs total they get pregnant with
  field_young0 field_young1 field_young2 field_young3 field_young4 field_young5 field_young6 field_young7 field_young8 field_young9            ;;field metabolic rate of young non-reproducing individuals
  field_mal0 field_mal1 field_mal2 field_mal3 field_mal4 field_mal5 field_mal6 field_mal7 field_mal8 field_mal9                                ;;field metabolic rate of male reproducing individuals
  field_preg0 field_preg1 field_preg2 field_preg3 field_preg4 field_preg5 field_preg6 field_preg7 field_preg8 field_preg9                      ;;field metabolic rate of female reproducing individuals
  bas_young0 bas_young1 bas_young2 bas_young3 bas_young4 bas_young5 bas_young6 bas_young7 bas_young8 bas_young9                                ;;basal metabolic rate of young non-reproducing individuals
  bas_mal0 bas_mal1 bas_mal2 bas_mal3 bas_mal4 bas_mal5 bas_mal6 bas_mal7 bas_mal8 bas_mal9                                                    ;;basal metabolic rate of male reproducing individuals
  bas_preg0 bas_preg1 bas_preg2 bas_preg3 bas_preg4 bas_preg5 bas_preg6 bas_preg7 bas_preg8 bas_preg9                                          ;;basal metabolic rate of female reproducing individuals
  repro_young0 repro_young1 repro_young2 repro_young3 repro_young4 repro_young5 repro_young6 repro_young7 repro_young8 repro_young9            ;;metabolic rate for reproduction of young non-reproducing individuals
  repro_mal0 repro_mal1 repro_mal2 repro_mal3 repro_mal4 repro_mal5 repro_mal6 repro_mal7 repro_mal8 repro_mal9                                ;;metabolic rate for reproduction of male reproducing individuals
  repro_preg0 repro_preg1 repro_preg2 repro_preg3 repro_preg4 repro_preg5 repro_preg6 repro_preg7 repro_preg8 repro_preg9                      ;;metabolic rate for reproduction of female reproducing individuals
  loco_young0 loco_young1 loco_young2 loco_young3 loco_young4 loco_young5 loco_young6 loco_young7 loco_young8 loco_young9 loco_fem0            ;;metabolic rate for locomotion of young non-reproducing individuals
  loco_mal0 loco_mal1 loco_mal2 loco_mal3 loco_mal4 loco_mal5 loco_mal6 loco_mal7 loco_mal8 loco_mal9                                          ;;metabolic rate for locomotion of male reproducing individuals
  loco_preg0 loco_preg1 loco_preg2 loco_preg3 loco_preg4 loco_preg5 loco_preg6 loco_preg7 loco_preg8 loco_preg9                                ;;metabolic rate for locomotion of female reproducing individuals
  grow_young0 grow_young1 grow_young2 grow_young3 grow_young4 grow_young5 grow_young6 grow_young7 grow_young8 grow_young9                      ;;metabolic rate for growth of young non-reproducing individuals
  grow_mal0 grow_mal1 grow_mal2 grow_mal3 grow_mal4 grow_mal5 grow_mal6 grow_mal7 grow_mal8 grow_mal9                                          ;;metabolic rate for growth of male reproducing individuals
  grow_preg0 grow_preg1 grow_preg2 grow_preg3 grow_preg4 grow_preg5 grow_preg6 grow_preg7 grow_preg8 grow_preg9                                ;;metabolic rate for growth of female reproducing individuals
  dig_young0 dig_young1 dig_young2 dig_young3 dig_young4 dig_young5 dig_young6 dig_young7 dig_young8 dig_young9                                ;;metabolic rate for digestion/heat increment of feeding of young non-reproducing individuals
  dig_mal0 dig_mal1 dig_mal2 dig_mal3 dig_mal4 dig_mal5 dig_mal6 dig_mal7 dig_mal8 dig_mal9                                                    ;;metabolic rate for digestion/heat increment of feeding of male reproducing individuals
  dig_preg0 dig_preg1 dig_preg2 dig_preg3 dig_preg4 dig_preg5 dig_preg6 dig_preg7 dig_preg8 dig_preg9                                          ;;metabolic rate for digestion/heat increment of feeding of female reproducing individuals
  move0 move1 move2 move3 move4 move5 move6 move7 move8 move9                                                                                  ;;mean movement distance per day
  field_max_preg0 field_max_preg1 field_max_preg2 field_max_preg3 field_max_preg4 field_max_preg5 field_max_preg6 field_max_preg7 field_max_preg8 field_max_preg9   ;;field metabolic rate of female individuals at the end of lactation
  mort_spec mort_food mort_stor mort_order mort_hr
  patchgrowths patchfeeds
  fail3_spec juv_fail juv_fail_spec feed_fail feed_fail_spec0 feed_fail_spec1 feed_fail_spec2 feed_fail_spec3 feed_fail_spec4 feed_fail_spec5 feed_fail_spec6 feed_fail_spec7 feed_fail_spec8 feed_fail_spec9
  failedpreg_spec failedlact_spec failage_spec failbold_spec pregs_spec

  fmr0 fmr1 fmr2 fmr3 fmr4 fmr5 fmr6 fmr7 fmr8 fmr9 loco0 loco1 loco2 loco3 loco4 loco5 loco6 loco7 loco8 loco9 repro0 repro1 repro2 repro3 repro4 repro5 repro6 repro7 repro8 repro9
  grow0 grow1 grow2 grow3 grow4 grow5 grow6 grow7 grow8 grow9 basal0 basal1 basal2 basal3 basal4 basal5 basal6 basal7 basal8 basal9 digest0 digest1 digest2 digest3 digest4 digest5 digest6 digest7 digest8 digest9
  prod0 prod1 prod2 prod3 prod4 prod5 prod6 prod7 prod8 prod9 in0 in1 in2 in3 in4 in5 in6 in7 in8 in9 balance0 balance1 balance2 balance3 balance4 balance5 balance6 balance7 balance8 balance9
  prodbalance0 prodbalance1 prodbalance2 prodbalance3 prodbalance4 prodbalance5 prodbalance6 prodbalance7 prodbalance8 prodbalance9 prod_sum
  inpatch0 inpatch1 inpatch2 inpatch3 inpatch4 inpatch5 inpatch6 inpatch7 inpatch8 inpatch9 inhr0 inhr1 inhr2 inhr3 inhr4 inhr5 inhr6 inhr7 inhr8 inhr9

  ;;cohort based output (only in long output)
  field_spec0 field_spec1 field_spec2 field_spec3 field_spec4 field_spec5 field_spec6 field_spec7 field_spec8 field_spec9
  bas_spec0 bas_spec1 bas_spec2 bas_spec3 bas_spec4 bas_spec5 bas_spec6 bas_spec7 bas_spec8 bas_spec9
  grow_spec0 grow_spec1 grow_spec2 grow_spec3 grow_spec4 grow_spec5 grow_spec6 grow_spec7 grow_spec8 grow_spec9
  repro_spec0 repro_spec1 repro_spec2 repro_spec3 repro_spec4 repro_spec5 repro_spec6 repro_spec7 repro_spec8 repro_spec9
  gest_spec0 gest_spec1 gest_spec2 gest_spec3 gest_spec4 gest_spec5 gest_spec6 gest_spec7 gest_spec8 gest_spec9
  loco_spec0 loco_spec1 loco_spec2 loco_spec3 loco_spec4 loco_spec5 loco_spec6 loco_spec7 loco_spec8 loco_spec9
  lact_spec0 lact_spec1 lact_spec2 lact_spec3 lact_spec4 lact_spec5 lact_spec6 lact_spec7 lact_spec8 lact_spec9
  digest_spec0 digest_spec1 digest_spec2 digest_spec3 digest_spec4 digest_spec5 digest_spec6 digest_spec7 digest_spec8 digest_spec9
  digest_mother_spec0 digest_mother_spec1 digest_mother_spec2 digest_mother_spec3 digest_mother_spec4 digest_mother_spec5 digest_mother_spec6 digest_mother_spec7 digest_mother_spec8 digest_mother_spec9
  mass_spec0 mass_spec1 mass_spec2 mass_spec3 mass_spec4 mass_spec5 mass_spec6 mass_spec7 mass_spec8 mass_spec9
  stor_spec0 stor_spec1 stor_spec2 stor_spec3 stor_spec4 stor_spec5 stor_spec6 stor_spec7 stor_spec8 stor_spec9
  fmr_pup_gest_spec0 fmr_pup_gest_spec1 fmr_pup_gest_spec2 fmr_pup_gest_spec3 fmr_pup_gest_spec4 fmr_pup_gest_spec5 fmr_pup_gest_spec6 fmr_pup_gest_spec7 fmr_pup_gest_spec8 fmr_pup_gest_spec9
  fmr_pup_lact_spec0 fmr_pup_lact_spec1 fmr_pup_lact_spec2 fmr_pup_lact_spec3 fmr_pup_lact_spec4 fmr_pup_lact_spec5 fmr_pup_lact_spec6 fmr_pup_lact_spec7 fmr_pup_lact_spec8 fmr_pup_lact_spec9
  young_mass_spec0 young_mass_spec1 young_mass_spec2 young_mass_spec3 young_mass_spec4 young_mass_spec5 young_mass_spec6 young_mass_spec7 young_mass_spec8 young_mass_spec9
]

patches-own [habitat patchfeed patchquali spec-id spec-list eaten save visited patchgrowth]
turtles-own [mass species shelter food-pref for-type maxhr hrsize lococost feedrate age average_age preg
             pregcount age_first_repro gest_period lact_period sex order stime hunger core forage_large forage_small storage max_storage
             young mass-neonate mass-weaning mass-embryo real-mass real-young focal share competition forage-patches spec_ind_list
             fmr dailydist younggrowth daymoved fmr_basal fmr_growth fmr_repro fmr_loco juvs juvshr preg_next fmr_digest postural postural-time feedspeed vis_before move-factor-allometric dom just-immi
             production intake balance balanceprod]

to setup
  clear-all
  reset-ticks
  set-parameters

  setup-patches-algo

  setup-turtles

  ask patches
  [
    set pcolor scale-color black habitat 1 0
    set patchgrowth 0
    set patchgrowths []
    set patchfeeds []
  ]
  find-hr

  ask one-of turtles [set focal 1]
  init-out
end

;---------------------------------------------

to set-parameters                                      ;; characterize species und basic parameters
  let addmass 0
  ;set maxmass 0.1                                      ;; upper limit of mean body mass
  set s0mass 0.1 * maxmass + addmass                             ;; mean body mass of species 0 in kg - normal distributed
  set s0massdev 0.2 * s0mass                           ;; std dev of body mass species 0 - normal distribution
  set s1mass 0.2 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s1massdev 0.2 * s1mass                           ;; std dev of body mass species 0 - normal distribution
  set s2mass 0.3 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s2massdev 0.2 * s2mass                           ;; std dev of body mass species 0 - normal distribution
  set s3mass 0.4 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s3massdev 0.2 * s3mass                           ;; std dev of body mass species 0 - normal distribution
  set s4mass 0.5 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s4massdev 0.2 * s4mass                           ;; std dev of body mass species 0 - normal distribution
  set s5mass 0.6 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s5massdev 0.2 * s5mass                           ;; std dev of body mass species 0 - normal distribution
  set s6mass 0.7 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s6massdev 0.2 * s6mass                           ;; std dev of body mass species 0 - normal distribution
  set s7mass 0.8 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s7massdev 0.2 * s7mass                           ;; std dev of body mass species 0 - normal distribution
  set s8mass 0.9 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s8massdev 0.2 * s8mass                           ;; std dev of body mass species 0 - normal distribution
  set s9mass maxmass + addmass                                    ;; mean body mass of species 0 - normal distributed
  set s9massdev 0.2 * s9mass                           ;; std dev of body mass species 0 - normal distribution

  set s0shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s1shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s2shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s3shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s4shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s5shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s6shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s7shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s8shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s9shelt 1;0.5                                      ;; index 0..1 of shelter use/need

  set s0food "omni"                                    ;; food preference; not distinguished yet
  set s1food "omni"                                    ;; food preference
  set s2food "omni"                                    ;; food preference
  set s3food "omni"                                    ;; food preference
  set s4food "omni"                                    ;; food preference
  set s5food "omni"                                    ;; food preference
  set s6food "omni"                                    ;; food preference
  set s7food "omni"                                    ;; food preference
  set s8food "omni"                                    ;; food preference
  set s9food "omni"                                    ;; food preference

  set s0fortype "central"                              ;; type of forage movement; not distinguished yet
  set s1fortype "central"                              ;; type of forage movement
  set s2fortype "central"                              ;; type of forage movement
  set s3fortype "central"                              ;; type of forage movement
  set s4fortype "central"                              ;; type of forage movement
  set s5fortype "central"                              ;; type of forage movement
  set s6fortype "central"                              ;; type of forage movement
  set s7fortype "central"                              ;; type of forage movement
  set s8fortype "central"                              ;; type of forage movement
  set s9fortype "central"                              ;; type of forage movement

  set actcov cover                                     ;; actual cover after changes
  set fail 0                                           ;; counting failing attempts by offspring to establish homerange
  set fail3 0                                          ;; counting death events caused by food shortage in attempts to re-establish homerange next day
  set failed-lact 0
  set repro 0                                          ;; counting overall annual reproduction

;;parameters to be changed
  set cover cover_percentage * 101 * 101                ;; xx % of shrub cover * nr of patches; xx * 100 * 100
  set feed-matrix 0.0                                   ;; no food prod in matrix
  set feed-struct 1.7127 * feed-amount ;6.85 ;1.7       ;; 10% of 68.5g dry mass/ grid cell*day Buchmann et al. 2011
  set mass_order prio_fact                              ;; % of turtles that are ordered acording to body mass, remaining turtles are ordered randomly
  set n_immi 0                                          ;; nr of random immigrants per timestep
  set focal_spec 1                                      ;; fcoal species of which daily homeranges are shown (cumulative)
  set patch-move 0
  set bold_prob 1
  set spec_num 10
  set high-season True
  set growthcost (7 + 6) * 1000 / 10
  set var 0
end
;---------------------------------------------

to setup-patches-algo                                  ;; creates a clumped landscape with two types of patches
  let cov 0
  ask patches   [
      set habitat 0
      set patchquali 0
      set spec-id -1
      set spec-list [ ]
      set eaten 0
      set visited 0
   ]
  while [cov < cover]
  [
   ask patches
   [
    ifelse ((habitat > 0) and (any? neighbors with [habitat = 0]) and (cov < cover))
    [ask one-of neighbors with [habitat = 0] [set habitat 1 set cov cov + 1]]
    [if random-float 1 < (1 - clump) and habitat = 0 and (cov < cover) [set habitat 1 set cov cov + 1]]
   ]
  ]
  ask patches[set save (sum [habitat] of neighbors > 6)]

  resource-update
end
;..............................................

to resource-update                                    ;;daily update of food resources, random variation normal distributed
  ask patches
   [
      ifelse habitat = 0 [set patchfeed feed-matrix] [set patchfeed (random-normal feed-struct feed-var) ;/ ((day + 100) / 100) ;random-normal feed-struct feed-var;0.7
      set spec-list [ ]]
   ]
end
;---------------------------------------------

to setup-turtles                                       ;; create xx initial individuals, location at 0,0; identifies species, mass,
  crt starting_indivs
                                                       ;; scenario 10 species, species identity 0 .. 9
  let slope -1.5                                       ;; allometric frequency distribution with mass^-1.5, see Buchmann et al. 2012 Ecography
  let s0m s0mass ^ (slope)
  let s1m s1mass ^ (slope)
  let s2m s2mass ^ (slope)
  let s3m s3mass ^ (slope)
  let s4m s4mass ^ (slope)
  let s5m s5mass ^ (slope)
  let s6m s6mass ^ (slope)
  let s7m s7mass ^ (slope)
  let s8m s8mass ^ (slope)
  let s9m s9mass ^ (slope)
  let mass-sum s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m + s8m + s9m

  ask turtles
  [
   let zuf random-float 1.0

   ifelse zuf < s0m / mass-sum
   [ set species 0 ]
    [ ifelse zuf < (s0m + s1m) / mass-sum [set species 1 ]
      [ ifelse zuf < (s0m + s1m + s2m) / mass-sum [set species 2 ]
        [ ifelse zuf < (s0m + s1m + s2m + s3m) / mass-sum [set species 3 ]
          [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m) / mass-sum [set species 4 ]
            [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m) / mass-sum [set species 5 ]
              [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m) / mass-sum [set species 6 ]
                [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m) / mass-sum [set species 7 ]
                  [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m + s8m) / mass-sum [set species 8 ]
                    [set species 9] ;;last else
                  ];;ifelse8
                ];;ifelse7
              ];;ifelse6
            ];;ifelse5
          ];; ifelse4
        ] ;; ifelse3
      ] ;; ifelse2
    ] ;;ifelse1
    if specs-included != "all" and specs-included != "two" [set species specs-included]
    if specs-included = "two" [if species != species1 and species != species2 [die]]
    specify_turtle
    if starting_indivs = 1 [set sex "female"]
    ;set sex "female"
  ] ;;ask turtles
end
;.............................................

to setup-turtles-random
  crt starting_indivs
  [
    set age random 100
    set species random 10
    if specs-included != "all" [set species specs-included]
    specify_turtle
    if starting_indivs = 1 [set sex "female"]
  ]
end
;.............................................

to specify_turtle                                               ;;characterize and parameterize individuals
  set hrsize 0
  let shy 0
  let zuf1 random-float 1.0
  set color (species * 10 + 5)
  ifelse (zuf1 < bold_prob) [set shy 0][set shy 1]              ;; shy=0 means bold individuals

    if (species = 0) [set mass random-normal s0mass s0massdev   ;; body mass normal distributed
      set shelter s0shelt * shy
      set food-pref s0food
      set for-type s0fortype]
    if (species = 1) [set mass random-normal s1mass s1massdev
      set shelter s1shelt * shy
      set food-pref s1food
      set for-type s1fortype]
    if (species = 2) [set mass random-normal s2mass s2massdev
      set shelter s2shelt * shy
      set food-pref s2food
      set for-type s2fortype]
    if (species = 3) [set mass random-normal s3mass s3massdev
      set shelter s3shelt * shy
      set food-pref s3food
      set for-type s3fortype]
    if (species = 4) [set mass random-normal s4mass s4massdev
      set shelter s4shelt * shy
      set food-pref s4food
      set for-type s4fortype]
    if (species = 5) [set mass random-normal s5mass s5massdev
      set shelter s5shelt * shy
      set food-pref s5food
      set for-type s5fortype]
    if (species = 6) [set mass random-normal s6mass s6massdev
      set shelter s6shelt * shy
      set food-pref s6food
      set for-type s6fortype]
    if (species = 7) [set mass random-normal s7mass s7massdev
      set shelter s7shelt * shy
      set food-pref s7food
      set for-type s7fortype]
    if (species = 8) [set mass random-normal s8mass s8massdev
      set shelter s8shelt * shy
      set food-pref s8food
      set for-type s8fortype]
    if (species = 9) [set mass random-normal s9mass s9massdev
      set shelter s9shelt * shy
      set food-pref s9food
      set for-type s9fortype]

    if mass < 0 [set mass 0 die]
    ;set shelter (1 - shelter)

    set preg 0
    set pregcount 0

    set average_age 1766.53 * (mass ^ 0.21)                      ;; average lifespan: allometric formula of mammals after Hamilton et al 2011 [g, days]
    set age_first_repro 293.17 * (mass ^ 0.27)                   ;; age of first reproduction: allometric fromula after Hamilton et al 2011
    set gest_period 64.14 * (mass ^ 0.24)                        ;; gestation period allometric after Hamilton et al 2011
    set lact_period  57.16 * (mass ^ 0.22)                       ;; lactation period allometric after Hamilton et al 2011
    set hunger 0
    set forage_small 0
    set forage_large 0

    set move-factor-allometric move-factor ;(mass ^ -0.3)

    set age random (average_age - 180) + 180
    let zuf random-float 1.0
    ifelse zuf < 0.5 [set sex "male"] [set sex "female"]
    let zuf# random-float 1.0
    ifelse zuf# < mass_order [set dom 1] [set dom 0]

    life-history
    set real-mass mass
    set storage 0
    set mass-embryo 0
    set spec_ind_list []
    set juvs 0
    set juvshr 0

    calc-maxhr
    calc-distance
    calc-lococost
    calc-feedrate real-mass
    calc-energetics
    calc-distance
    ;life-history

  if age > age_first_repro                                       ;; initialization with some already pregnant individuals
  [
  if sex = "female"
  [
  set preg random 2
  if preg = 1
  [
    set real-young young
    set pregcount random gest_period + lact_period
    repeat min (list gest_period pregcount)
    [
      embryo
    ]
    if pregcount > gest_period
    [
      repeat pregcount - gest_period
      [
        juvenile
      ]
    ]
  ]
    let young-mass real-young * mass-embryo
    calc-feedrate (young-mass + real-mass)
    set max_storage storage_addition * (294.8) * ((real-mass + young-mass) ^ 1.19) + feedrate
  ]
  ]

end
;---------------------------------------------
to calc-energetics
    set max_storage storage_addition * (294.8) * (real-mass ^ 1.19) + feedrate ;; allometric storage after Lindstedt and Boyce 1985
    set share (real-mass / 0.001) ^ (-0.25)                           ;; allometric share of available food see Buchmann 2011
end
;----------------------------------------------

to calc-maxhr                                                      ;; calculate max home range after Kelt & Van Vuren 2001 in ha, see Buchmann et al 2011
  let maxhrx1 56.23 * (real-mass ^ 0.91)                                 ;; max hr for herbivores and omnivores, larger one is used
  let maxhrx2 47.863 * (real-mass ^ 1.18)
  ifelse maxhrx2 > maxhrx1 [set maxhr maxhrx2] [set maxhr maxhrx1]
  set maxhr maxhr * 10000                                          ;; in m2
  set maxhr sqrt (maxhr / pi)                                      ;; radius of maxhr in m
  set maxhr maxhr / 10                                             ;; radius in patch length (=10m)
end
;---------------------------------------------

to calc-lococost                                                   ;; calculate movement costs, specific for food type
  calc-postural
  set lococost 10.7 * (real-mass ^ 0.68) + postural                                ;; costs for mammals in J/m; mass in kg after Calder 1996
  set lococost lococost / 10000                                    ;; costs in g dry biomass/m; after Nagy'99 p.263 Buchmann 2011
  set lococost lococost * 10                                       ;; lococost in patchlength (= 10m)
end
;---------------------------------------------

to calc-postural
  calc-postural-time
  set postural postural-time / (0.01398 * ((real-mass * 1000) ^ 0.217) * (e ^ (-0.002 * ((log (real-mass * 1000) e) ^ 2 )))) ;J/m
end
;----------------------------------------------

to calc-postural-time
  set postural-time (6.03 * (real-mass ^ 0.697) - 2.963 * (real-mass ^ 0.737)) ;J/s
end
;-----------------------------------------------

to calc-feed-postural
  set feedspeed postural-time * (0.71 * (real-mass ^ 0.7) / 60) ;J/g - (g/s Shipley 1994)
  set feedspeed feedspeed * 100 ; g/g
end
;----------------------------------------------

to calc-feedrate [m]                                                   ; basal metabolic rate from Savage et al. 2004  ;; caluclate daily feeding rate, specific for food type
 set feedrate 25.6 * (m ^ 0.737);
end
;----------------------------------------------

to calc-distance
  let dailydist_mean 1.038 * (real-mass ^ 0.25) * 100              ;; allometric daily movement distance Garland 1983 scaled to 10m patches
  set dailydist random-normal dailydist_mean (dailydist_mean * 0.1)
end
;---------------------------------------------

to life-history
  set young round (2.24 * (mass ^ (-0.13)))                       ;; allomatric number of offspring after Hamilton et al. 2011
  set mass-neonate 47.86 * (mass ^ 0.93)                          ;; allometric neonate mass after Hamilton 2011
  set mass-weaning 295.12 * (mass ^ 0.91)                         ;; allometric weaning mass after Hamilton 2011
end
;---------------------------------------------

to embryo                                                         ;; embryo growth after Rickleffs 2010
  let A exp(0.865 + 1.006 * log  mass-neonate e)
  let M0 0.0001
  let b log (A / M0) e
  let k exp(0.627 - 0.905 * log gest_period e)
  set mass-embryo A * exp(- b * exp(- k * pregcount)) / 1000
end
;-----------------------------------------------

to juvenile                                                        ;; lactation growth after Sibly 2013
  let m mass * 1000
  let jm mass-embryo * 1000
  let factor ((m / jm) ^ (1 / 3) - 1)
  let prod factor * 3 / lact_period * log ((1 - ((mass-neonate / m) ^ (1 / 3))) / (1 - ((mass-weaning / m) ^ (1 / 3)))) e
  set jm prod * jm
  set mass-embryo mass-embryo + jm / 1000
end

;-----------------------------------------------

to growth                                                          ;; ongoing growth after Sibly 2013
  let m mass * 1000
  let realm real-mass * 1000
  let factor ((m / realm) ^ (1 / 3) - 1)
  let prod factor * 3 / lact_period * log ((1 - ((mass-neonate / m) ^ (1 / 3))) / (1 - ((mass-weaning / m) ^ (1 / 3)))) e
  set realm prod * realm
  set real-mass min list mass (real-mass + realm / 1000)

  calc-distance
  calc-lococost
  calc-feedrate real-mass
  calc-energetics
end

;-----------------------------------------------

to change_order_density_radius                                      ;; slight ordering of individuals for feeding
  ask turtles
  [
     let zuf random-float 1.0
     let sp species
     let mhr maxhr
     ifelse zuf < mass_order [set order mass * (1 - count turtles in-radius mhr with [species = sp] / count turtles in-radius mhr)] [set order random-float maxmass] ;;
  ]
end
;-----------------------------------------------

to one-step                                                      ;; one time step

  if ( debug = 1 ) and ( ticks > 360 ) [
    profiler:start
  ]

  set day ticks mod 360
  set fail 0
  set fail3 0
  set fail3_spec [0 0 0 0 0 0 0 0 0 0]
  set repro 0

  ask patches
  [
    set pcolor scale-color green habitat 1 0
    set spec-list [ ]
    set eaten 0
    set visited 0
  ]

  ask turtles
  [
   set age age + 1
   set fmr 0
   set preg_next false
  ]

  resource-update                                                  ;; function: update resources

  change_order_density_radius                                      ;; function: order mass_order% of turtles according to body mass
  check-hr                                                         ;; function: check if old homerange is still O.K., adapt if possible

  ask turtles
  [
    set competition sum [length spec-list] of patches in-radius (hrsize / 2)
    set forage-patches count (patches in-radius (hrsize / 2) with [habitat = 1])
  ]

  maintenance                                                      ;; function: energy allocation to maintenance, growth, juveniles
  foreach sort-on [ (- order)] turtles                             ;; offspring of heaviest turtles check hr first for mass_oder% of turtles
  [ ?1 -> ask ?1                                                   ;; ask turtles
   [
   if sex = "female"
    [
     ifelse preg = 1
      [set pregcount pregcount + 1]
        [let spec species ;
          if (real-mass >= maturity * mass) [set preg 1 set real-young young set mass-embryo 0]] ; and day > 91 and day < 315  and any? turtles in-radius maxhr with [species = spec and sex = "male"];(age > age_first_repro) [set preg 1 set real-young young] ];; and storage > 0.2 * feedrate) [ set preg 1 set real-young young] ]                   ;; getting pregnant deterministically => age allometric
     if pregcount >  gest_period + lact_period                     ;; days of gestation/pregnancy PLUS days of lactation==> allometric after Hamilton et al 2011
      [ offspring set pregcount 0 set preg_next true ]
    ] ;; end if sex = female
   mort                                                            ;; mortality function
   ] ;;end ?1 -> ask ?1
   ] ;; end foreach sort

  out

  output
  ask turtles with [preg_next = true] [set preg 0]

  if ( debug = 1 ) and ( ticks > 360 ) [
    profiler:stop
    print profiler:report
    profiler:reset
  ]
end
;-----------------------------------------------

to go                                                              ;; to go function, multiple time steps
  one-step
  tick
end

;---------------------------------------------

to lifetime_success  ;;if a female dies, report how many successfull juveniles it had over its lifetime
  if sex = "female"
  [
  if species = 0 [set rep_success_0 lput (juvs) rep_success_0 set rep_success_hr_0 lput (juvshr) rep_success_hr_0]
  if species = 1 [set rep_success_1 lput (juvs) rep_success_1 set rep_success_hr_1 lput (juvshr) rep_success_hr_1]
  if species = 2 [set rep_success_2 lput (juvs) rep_success_2 set rep_success_hr_2 lput (juvshr) rep_success_hr_2]
  if species = 3 [set rep_success_3 lput (juvs) rep_success_3 set rep_success_hr_3 lput (juvshr) rep_success_hr_3]
  if species = 4 [set rep_success_4 lput (juvs) rep_success_4 set rep_success_hr_4 lput (juvshr) rep_success_hr_4]
  if species = 5 [set rep_success_5 lput (juvs) rep_success_5 set rep_success_hr_5 lput (juvshr) rep_success_hr_5]
  if species = 6 [set rep_success_6 lput (juvs) rep_success_6 set rep_success_hr_6 lput (juvshr) rep_success_hr_6]
  if species = 7 [set rep_success_7 lput (juvs) rep_success_7 set rep_success_hr_7 lput (juvshr) rep_success_hr_7]
  if species = 8 [set rep_success_8 lput (juvs) rep_success_8 set rep_success_hr_8 lput (juvshr) rep_success_hr_8]
  if species = 9 [set rep_success_9 lput (juvs) rep_success_9 set rep_success_hr_9 lput (juvshr) rep_success_hr_9]
  ]
end

;---------------------------------------------


to offspring                                                         ;; offspring, inherits everything from mother
  set repro repro + real-young                                            ;; to count overall annual reproduction
  let f focal
  let stor storage / (real-young + 1)
  set storage stor                                                   ;; mothers storage equally divided to mother and offspring
  set juvs juvs + real-young
  set pregs_spec replace-item species pregs_spec (item species pregs_spec + 1)
  let juvs-found-hr 0
  hatch real-young
    [
     specify_turtle
     set age 0
     set preg 0
     set pregcount 0
     set real-mass (mass-weaning / 1000)
     calc-lococost                                                   ;; note: has to be calculated again according to real mass
     calc-feedrate real-mass
     calc-distance
     calc-energetics
     set storage stor
     set focal f
     find-hr-offspring                                               ;; note: to offspring takes place after gestation AND lactation!
     set juvs-found-hr juvs-found-hr + 1
    ]
  set juvshr juvshr + juvs-found-hr
  calc-feedrate real-mass                           ;; adapt mothers feedrate and max storage without pregnancy
  set max_storage storage_addition * (294.8) * (real-mass ^ 1.19) + feedrate

end
;.............................................

to mort                                                              ;; mortality function related to life span, normally distributed and boldness
  if shelter = 0 and random-float 1 < mortbold [ lifetime_success set failbold_spec replace-item species failbold_spec (item species failbold_spec + 1) die ]    ;; additional mortality for bold indivdiuals (bold=> shelter= 0)
  if age > random-normal average_age ( 0.1 * average_age ) [lifetime_success set failage_spec replace-item species failage_spec (item species failage_spec + 1) die ]   ;; average life span is allometric, see above
end
;.............................................

to find-hr                                                         ;; find suitable homerange for initial distribution
  foreach sort-on [ (- order)] turtles
  [ ?1 -> ask ?1
  [
  set maxrad maxhr
  set minfeed feedrate
  set foodshare share
  set shelter-need shelter
  set spec-color (species * 10 + 5)
  set movecost lococost
  set move-factor move-factor-allometric
  ;set feedtime feedspeed
  let spec species
  let success 0
  let moved 0
  let try 0
  let maxmove dailydist
  let stor storage
  while [success = 0 and try < hr_try_init]                                ;; 100 attempts for each mammal of inital community to find suitable hr
    [
      set try try + 1
      set r-search 0
      set feed 0
      set moved 0
      move-to one-of patches with [habitat > 0]                   ;; random search for potential hr core cell
      set core patchquali
      set feed digestion * patchfeed * foodshare - patch-move * movecost                               ;; calculate food of core cell
      set eaten 1
      if (storage + feed) >= max_storage ;if feed >= minfeed                                           ;; if enough food in core cell - end search
       [
         set success 1
         set eaten 1                                               ;; 1 indicates patch as hr-patch forlater food reduction
         set spec-id spec                                          ;; identify patch as part of hr of species
         set spec-list lput spec-id spec-list                      ;; add species to patch-specific species list
       ]
      while [r-search < maxhr and (storage + feed) < max_storage]                       ;; if not enough food in core cell - search in neighborhood
      [
       set r-search r-search + 1
       ask patches in-radius r-search with [patchfeed > 0 and eaten = 0]
        [
          if (stor + feed) >  (movecost * distance myself) and eaten = 0
          [
           set eaten 2                                             ;; 2 indicates potential use as hr-patch for later food reduction
           ifelse save                  ;; edge effects: patches 'in the open' are less frequently visited or shorter time => reduced food intake
            [set moved moved + 2 * distance myself * move-factor + patch-move
              set feed feed + digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost] ; - feedtime * patchfeed * foodshare]
            [set moved moved + (2 * distance myself * move-factor + patch-move) * (1 - shelter-need)
              set feed feed + (digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost) * (1 - shelter-need)] ;- feedtime * patchfeed * foodshare) * (1 - shelter-need)]
        ]  ;; end ask-patches in r
       ]  ;; end while r-search
      ]

     ifelse storage + feed >= minfeed
       [
          set success 1
          set storage min (list max_storage (feed - minfeed))
          set patchfeed patchfeed * (1 - foodshare)
          set eaten 0
          ask patches in-radius (r-search + 1) with [eaten = 2]
           [
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]        ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter-need))]
              set spec-id spec                                     ;; identify patch as part of hr of species
              set spec-list lput spec-id spec-list                 ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec-color]                   ;; show only hr of focal species
              set eaten 0
           ] ;; end ask parches in-radius
          let patchafter sum [patchfeed] of patches
       ] ;; end ifelse feed >=minfeed cond1
       [
          ask patches in-radius (r-search) [set eaten 0]            ;; set back to not-eaten
       ] ;; end ifelse feed >=minfeed cond2

    ] ;; end while success = 0

    if success = 1 [
      set hrsize 2 * r-search
      set daymoved moved
      set fmr_digest 0.1 * (feed + moved * movecost)
      set fmr_loco moved * movecost
      set storage min (list max_storage (storage + feed))]
    if success = 0 [die]

  ] ;; end ask turtles
  ]
end
;----------------------------------------------

to mother-feedrate
ask turtles
  [
    if age > 0                                                         ;; only for non-offspring
     [
     let young-mass 0
     if preg = 1                                                       ;; gestating or lactating females
      [
       let ym mass-embryo
       ifelse (pregcount < gest_period)                                  ;; gestating
        [
         embryo
         set young-mass real-young * mass-embryo
        ]
      [
        if (pregcount < gest_period + lact_period)
         [
          juvenile
          set young-mass real-young * mass-embryo
         ]
      ] ;;end ifelse pregcount
          set younggrowth (mass-embryo - ym)
          calc-feedrate (young-mass + real-mass)
          set max_storage storage_addition * (294.8) * ((real-mass + young-mass) ^ 1.19) + feedrate
     ]
    ]
  ]
end

;----------------------------------------------
to check-hr                                                              ;; check if existing homerange is still sufficiant => food is reduced!
 mother-feedrate
 foreach sort-on [ (- order)] turtles
  [ ?1 -> ask ?1
   [
    if age > 0                                                           ;; only for non-offspring
    [
     let m real-mass
     set real-mass real-mass + storage / 3930
     calc-lococost
     set move-factor move-factor-allometric
     set real-mass m
     let spec species
     set maxrad maxhr
    set foodshare share
    set shelter-need shelter
    set movecost lococost
    let success 0
    set r-search 0
    set feed 0
    let stor storage
    set feed digestion * patchfeed * foodshare - patch-move * movecost   ;; calculate food of core cell
    set patchfeed patchfeed * (1 - foodshare)                            ;; reduce remaining food in core patch
    let in patchfeed * foodshare
    set eaten 1
    let moved 0
    let before visited
    let maxmove dailydist
    if patchfeed < 0 [set patchfeed 0]
    set visited visited + 1
    if (storage + feed) >= max_storage                                                   ;; if enough food in core cell - end
       [
         set success 1
         set hunger 0
         set pcolor spec-color
         set spec-id spec
         set spec-list lput spec-id spec-list
       ]
        while [r-search < maxhr and (storage + feed) < max_storage]                         ;; if not enough food in core cell - search in neighborhood
      [
       set r-search r-search + 1
         ask patches in-radius r-search with [eaten = 0 and patchfeed > 0]
        [
         if (stor + feed) >  (movecost * distance myself) and eaten = 0
         [
         if spec = focal_spec [set pcolor spec-color]                    ;; show only hr of focal species
         set eaten 2
         set spec-id spec
         set spec-list lput spec-id spec-list
         ifelse save                           ;; edge effect: patches 'in the open' are potentially less frequently or for shorter time visited
           [
             set moved moved + 2 * distance myself * move-factor + patch-move
             set feed feed + digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost
             set in in + patchfeed * foodshare
             set patchfeed patchfeed * (1 - foodshare)
             set before before + visited
             set visited visited + 1
           ]
           [
             set moved moved + (2 * distance myself * move-factor + patch-move) * (1 - shelter-need)
             set feed feed + (digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost) * (1 - shelter-need)
             set in in + (patchfeed * foodshare) * (1 - shelter-need)
             set patchfeed patchfeed * (1 - foodshare * (1 - shelter-need))  ;;reduce food in patch
             set before before + visited * (1 - shelter-need)
             set visited visited + 1 * (1 - shelter-need)
           ]
         if patchfeed < 0 [set patchfeed 0]
        ] ;; end ask patches in-radius
        ]
       ] ;; end while r-search
        set fmr_digest 0.1 * (feed / digestion + moved * movecost)
        set fmr_loco moved * movecost
        set fmr fmr + fmr_loco + fmr_digest
        set daymoved moved
        set storage min (list max_storage (storage + feed))
        set hrsize 2 * r-search
        set vis_before before
        set intake in
        ask patches in-radius r-search [set eaten 0]

      ]
    ]
  ]
end
;----------------------------------------------------------
to maintenance
ask turtles
  [
    set production 0
    set fmr_repro 0
    set fmr_growth 0
    if age > 0                                                         ;; only for non-offspring
     [
     let young-mass 0
     let synth 0
     if preg = 1                                                      ;; gestating or lactating females
      [
       set young-mass real-young * mass-embryo
       set synth (real-young * younggrowth * growthcost)
      ]
      let old feedrate
      calc-feedrate (young-mass + real-mass)
      let act_feedrate feedrate
      set feedrate old
      ifelse storage > (act_feedrate + synth)                            ;; if generally enough energy, try to grow
        [
          let m real-mass
          growth
          let diff real-mass - m
          let diffcost diff * growthcost
          set storage storage - act_feedrate - synth
          set production production + (synth / growthcost)
          ifelse storage > diffcost
          [set storage storage - diffcost
           set production production + diff
           calc-feedrate (real-mass + young-mass)
           set act_feedrate feedrate
           calc-feedrate (real-mass)
           set fmr fmr + act_feedrate + diffcost + synth
           set fmr_basal feedrate
           set fmr_repro synth + (act_feedrate - feedrate)
           set fmr_growth diffcost]
          [set real-mass m
              calc-distance
              calc-lococost
              calc-feedrate (young-mass + real-mass)
              calc-energetics
              set act_feedrate feedrate
              calc-feedrate (real-mass)
              set fmr_basal feedrate
              set fmr_repro synth + (act_feedrate - feedrate)
              set fmr fmr + act_feedrate + synth
            ]
        ]
      [
        ifelse preg = 1
        [
          ifelse pregcount > gest_period
            [
                while [storage < act_feedrate + (real-young * younggrowth * growthcost) and real-young > 0]  ;;if not enough energy and in lactation period, loose one juvenile after another
              [
              set real-young real-young - 1
              set failed-lact failed-lact + 1
              set failedlact_spec replace-item species failedlact_spec (item species failedlact_spec + 1)
              set young-mass real-young * mass-embryo
                set old feedrate
                calc-feedrate (young-mass + real-mass)
                set act_feedrate feedrate
                set feedrate old
              ]
              ifelse storage > act_feedrate + (real-young * younggrowth * growthcost)
              [set storage storage - act_feedrate - (real-young * younggrowth * growthcost)
               set production production + real-young * younggrowth
               set fmr fmr + act_feedrate + (real-young * younggrowth * growthcost)
               set act_feedrate feedrate
               calc-feedrate (real-mass)
               set fmr_basal feedrate
               set fmr_repro (real-young * younggrowth * growthcost) + (act_feedrate - feedrate)
               if real-young = 0 [set preg 0 set pregcount 0 set failedpreg_spec replace-item species failedpreg_spec (item species failedpreg_spec + 1)]
              ]
              [                                                                                                 ;;if not enough energy after loosing all juveniles, die
                set mort_order replace-item species mort_order (sentence order (item species mort_order))
                set mort_stor replace-item species mort_stor (sentence (storage / max_storage) (item species mort_stor))
                set mort_food replace-item species mort_food (sentence (competition / forage-patches) (item species mort_food))
                set mort_hr replace-item species mort_hr (sentence hrsize (item species mort_hr))
                set fail3_spec replace-item species fail3_spec (item species fail3_spec + 1)
                set preg 0
                set pregcount 0
                set storage 0
                set fail3 fail3 + 1
                lifetime_success
                die
              ]
            ]
        [                                                                                                      ;;if not enough energy and in gestation, loose pregnancy
          set preg 0
          set pregcount 0
          set failedpreg_spec replace-item species failedpreg_spec (item species failedpreg_spec + 1)
          calc-feedrate real-mass
          ifelse storage > feedrate
          [set storage storage - feedrate
           set fmr fmr + feedrate
           set fmr_basal feedrate]
          [                                                                                                     ;;if not enough energy after loosing pregnancy, die
           set mort_order replace-item species mort_order (sentence order (item species mort_order))
           set mort_stor replace-item species mort_stor (sentence (storage / max_storage) (item species mort_stor))
           set mort_food replace-item species mort_food (sentence (competition / forage-patches) (item species mort_food))
           set mort_hr replace-item species mort_hr (sentence hrsize (item species mort_hr))
           set fail3_spec replace-item species fail3_spec (item species fail3_spec + 1)
           set storage 0
           set fail3 fail3 + 1
           lifetime_success
           die
          ]
        ]
        ]
        [                                                                                                        ;;if not enough energy, die
           set mort_order replace-item species mort_order (sentence order (item species mort_order))
           set mort_stor replace-item species mort_stor (sentence (storage / max_storage) (item species mort_stor))
           set mort_food replace-item species mort_food (sentence (competition / forage-patches) (item species mort_food))
           set mort_hr replace-item species mort_hr (sentence hrsize (item species mort_hr))
           set fail3_spec replace-item species fail3_spec (item species fail3_spec + 1)
           set storage 0
           set fail3 fail3 + 1
           lifetime_success
           die
        ]

      ]
    if intake > 0
        [
        set balance fmr / (intake * digestion)
        set balanceprod production / intake
        ]
    ]
  ]
  ;]
end
;----------------------------------------------

to find-hr-offspring                                                         ;; offspring searches for new homerange
  set maxrad maxhr
  set minfeed feedrate;;
  set foodshare share
  set shelter-need shelter
  set spec-color (species * 10 + 5)
  set movecost lococost
  set move-factor move-factor-allometric
  set max_nat_disp 3.31 * (mass ^ 0.65)                                      ;; allometric maximum natal dispersal after Sutherland et al (2000) in km
  let spec species
  let success 0
  let try 0
  let xx xcor
  let yy ycor
  let moved 0
  let maxmove dailydist
  let stor storage
  while [success = 0 and try < hr_try_juv]                                   ;; hr_try_juv attempts for each mammal to find suitable hr
    [
      setxy xx yy
      set try try + 1
      let patch-before patch-here
      move-to one-of patches in-radius (100 * max_nat_disp) with [habitat > 0]
      set core patchquali
      set r-search 0
      set moved 0
      set feed digestion * patchfeed * foodshare - lococost * distance patch-before - patch-move * movecost;- feedtime * patchfeed * foodshare                                        ;; take food of core cell
      set eaten 1
      if (storage + feed) >= max_storage                                                    ;; if enough food in core cell - end
       [
         set success 1
         set pcolor spec-color
         set eaten 1                                                         ;; 1 indicates patch as hr-patch for later food reduction
         set spec-id spec                                                    ;; identify patch as part of hr of species
         set spec-list lput spec-id spec-list                                ;; add species to patch-specific species list
       ]
      while [r-search < maxhr and (storage + feed) < max_storage] ;[r-search < maxrad]
     [
       set r-search r-search + 1
       ask patches in-radius (r-search) with [patchfeed > 0 and eaten = 0]
        [
          if (stor + feed) >  (movecost * distance myself) and eaten = 0
          [
          set eaten 2                                                        ;;  2 indicates potential use as hr-patch for later food reduction
          ifelse save                              ;; edge effect: patches 'in the open' are less frequently visited or fot shorter time
            [set moved moved + 2 * distance myself * move-factor + patch-move
              set feed feed + digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost]
            [set moved moved + (2 * distance myself * move-factor + patch-move) * (1 - shelter-need)
              set feed feed + (digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost) * (1 - shelter-need)]
        ] ;; end ask patches in radius
       ] ;; end while r-search
      ]

      ifelse (storage + feed) >= minfeed
       [
         set success 1
         set storage  min ( list max_storage (storage + feed - minfeed))
         set patchfeed patchfeed * (1 - foodshare)
         set eaten 0
         set daymoved moved
         ask patches in-radius (r-search + 1) with [eaten = 2]
           [
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]                  ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter-need))]
              set spec-id spec                                               ;; identify patch as part of hr of species
              set spec-list lput spec-id spec-list                           ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec-color]                   ;; show only hr of focal species
              set eaten 0
             ;; end if eaten = 1
           ] ;; end ask patches in-radius
         ] ;; end ifelse feed >=minfeed cond 1
         [
          ask patches in-radius (r-search) [set eaten 0]                     ;; set back to not-eaten
         ] ;;end ifelse feed >=minfeed cond 2
    ] ;; end while success = 0
   if success = 1 [set hrsize 2 * r-search set suc_juv replace-item species suc_juv (item species suc_juv + 1)]
   if success = 0 [set fail (fail + 1) set juv_fail_spec replace-item species juv_fail_spec (item species juv_fail_spec + 1) die]                                                    ;; unsuccessful turtles die or emigrate from area
end
;----------------------------------------------

to output  ;; write 'critical time to diversity loss' and shannion-diversity in file
  let specs [species] of turtles
  let specs_unique remove-duplicates specs
  set spec_num length specs_unique
end
;----------------------------------------------
;initialization of the output
to init-out
  set suc_juv [0 0 0 0 0 0 0 0 0 0]
  set suc_immi [0 0 0 0 0 0 0 0 0 0]
  set prod_sum [0 0 0 0 0 0 0 0 0 0]
  set juv_fail_spec [0 0 0 0 0 0 0 0 0 0]
  set failage_spec [0 0 0 0 0 0 0 0 0 0]
  set failbold_spec [0 0 0 0 0 0 0 0 0 0]
  set failedpreg_spec [0 0 0 0 0 0 0 0 0 0]
  set failedlact_spec [0 0 0 0 0 0 0 0 0 0]
  set pregs_spec [0 0 0 0 0 0 0 0 0 0]

  set rep_success_0 []
  set rep_success_1 []
  set rep_success_2 []
  set rep_success_3 []
  set rep_success_4 []
  set rep_success_5 []
  set rep_success_6 []
  set rep_success_7 []
  set rep_success_8 []
  set rep_success_9 []

  set rep_success_hr_0 []
  set rep_success_hr_1 []
  set rep_success_hr_2 []
  set rep_success_hr_3 []
  set rep_success_hr_4 []
  set rep_success_hr_5 []
  set rep_success_hr_6 []
  set rep_success_hr_7 []
  set rep_success_hr_8 []
  set rep_success_hr_9 []

  set number []
  set number0 []
  set number1 []
  set number2 []
  set number3 []
  set number4 []
  set number5 []
  set number6 []
  set number7 []
  set number8 []
  set number9 []

  set juv_fail []
  set feed_fail []

  set feed_fail_spec0 []
  set feed_fail_spec1 []
  set feed_fail_spec2 []
  set feed_fail_spec3 []
  set feed_fail_spec4 []
  set feed_fail_spec5 []
  set feed_fail_spec6 []
  set feed_fail_spec7 []
  set feed_fail_spec8 []
  set feed_fail_spec9 []

  set patches0 []
  set patches1 []
  set patches2 []
  set patches3 []
  set patches4 []
  set patches5 []
  set patches6 []
  set patches7 []
  set patches8 []
  set patches9 []

  set compet0 []
  set compet1 []
  set compet2 []
  set compet3 []
  set compet4 []
  set compet5 []
  set compet6 []
  set compet7 []
  set compet8 []
  set compet9 []

  set order0 []
  set order1 []
  set order2 []
  set order3 []
  set order4 []
  set order5 []
  set order6 []
  set order7 []
  set order8 []
  set order9 []

  set vis_before0 []
  set vis_before1 []
  set vis_before2 []
  set vis_before3 []
  set vis_before4 []
  set vis_before5 []
  set vis_before6 []
  set vis_before7 []
  set vis_before8 []
  set vis_before9 []

  set hr []
  set hr0 []
  set hr1 []
  set hr2 []
  set hr3 []
  set hr4 []
  set hr5 []
  set hr6 []
  set hr7 []
  set hr8 []
  set hr9 []

  set maxhrs []
  set maxhr0 []
  set maxhr1 []
  set maxhr2 []
  set maxhr3 []
  set maxhr4 []
  set maxhr5 []
  set maxhr6 []
  set maxhr7 []
  set maxhr8 []
  set maxhr9 []

  set countmaxhr []
  set countmaxhr0 []
  set countmaxhr1 []
  set countmaxhr2 []
  set countmaxhr3 []
  set countmaxhr4 []
  set countmaxhr5 []
  set countmaxhr6 []
  set countmaxhr7 []
  set countmaxhr8 []
  set countmaxhr9 []

  set ages []
  set age0 []
  set age1 []
  set age2 []
  set age3 []
  set age4 []
  set age5 []
  set age6 []
  set age7 []
  set age8 []
  set age9 []

  set stors []
  set stor0 []
  set stor1 []
  set stor2 []
  set stor3 []
  set stor4 []
  set stor5 []
  set stor6 []
  set stor7 []
  set stor8 []
  set stor9 []

  set abs_stor0 []
  set abs_stor1 []
  set abs_stor2 []
  set abs_stor3 []
  set abs_stor4 []
  set abs_stor5 []
  set abs_stor6 []
  set abs_stor7 []
  set abs_stor8 []
  set abs_stor9 []

  set rm []
  set rm0 []
  set rm1 []
  set rm2 []
  set rm3 []
  set rm4 []
  set rm5 []
  set rm6 []
  set rm7 []
  set rm8 []
  set rm9 []

  set fmr0 []
  set fmr1 []
  set fmr2 []
  set fmr3 []
  set fmr4 []
  set fmr5 []
  set fmr6 []
  set fmr7 []
  set fmr8 []
  set fmr9 []

  set loco0 []
  set loco1 []
  set loco2 []
  set loco3 []
  set loco4 []
  set loco5 []
  set loco6 []
  set loco7 []
  set loco8 []
  set loco9 []

  set repro0 []
  set repro1 []
  set repro2 []
  set repro3 []
  set repro4 []
  set repro5 []
  set repro6 []
  set repro7 []
  set repro8 []
  set repro9 []

  set grow0 []
  set grow1 []
  set grow2 []
  set grow3 []
  set grow4 []
  set grow5 []
  set grow6 []
  set grow7 []
  set grow8 []
  set grow9 []

  set basal0 []
  set basal1 []
  set basal2 []
  set basal3 []
  set basal4 []
  set basal5 []
  set basal6 []
  set basal7 []
  set basal8 []
  set basal9 []

  set digest0 []
  set digest1 []
  set digest2 []
  set digest3 []
  set digest4 []
  set digest5 []
  set digest6 []
  set digest7 []
  set digest8 []
  set digest9 []

  set prod0 []
  set prod1 []
  set prod2 []
  set prod3 []
  set prod4 []
  set prod5 []
  set prod6 []
  set prod7 []
  set prod8 []
  set prod9 []

  set in0 []
  set in1 []
  set in2 []
  set in3 []
  set in4 []
  set in5 []
  set in6 []
  set in7 []
  set in8 []
  set in9 []

  set inpatch0 []
  set inpatch1 []
  set inpatch2 []
  set inpatch3 []
  set inpatch4 []
  set inpatch5 []
  set inpatch6 []
  set inpatch7 []
  set inpatch8 []
  set inpatch9 []

  set inhr0 []
  set inhr1 []
  set inhr2 []
  set inhr3 []
  set inhr4 []
  set inhr5 []
  set inhr6 []
  set inhr7 []
  set inhr8 []
  set inhr9 []

  set balance0 []
  set balance1 []
  set balance2 []
  set balance3 []
  set balance4 []
  set balance5 []
  set balance6 []
  set balance7 []
  set balance8 []
  set balance9 []

  set prodbalance0 []
  set prodbalance1 []
  set prodbalance2 []
  set prodbalance3 []
  set prodbalance4 []
  set prodbalance5 []
  set prodbalance6 []
  set prodbalance7 []
  set prodbalance8 []
  set prodbalance9 []

  set hrpreg []
  set hrpreg0 []
  set hrpreg1 []
  set hrpreg2 []
  set hrpreg3 []
  set hrpreg4 []
  set hrpreg5 []
  set hrpreg6 []
  set hrpreg7 []
  set hrpreg8 []
  set hrpreg9 []

  set storpreg []
  set storpreg0 []
  set storpreg1 []
  set storpreg2 []
  set storpreg3 []
  set storpreg4 []
  set storpreg5 []
  set storpreg6 []
  set storpreg7 []
  set storpreg8 []
  set storpreg9 []

  set comp5 []
  set comp10 []
  set comp15 []
  set comp20 []
  set comp25 []
  set compmore []

  set comp5_stor []
  set comp10_stor []
  set comp15_stor []
  set comp20_stor []
  set comp25_stor []
  set compmore_stor []

  set comp5_fmr []
  set comp10_fmr []
  set comp15_fmr []
  set comp20_fmr []
  set comp25_fmr []
  set compmore_fmr []

  set stor_pregstart []
  set stor_pregend []

  set juv []
  set juv0 []
  set juv1 []
  set juv2 []
  set juv3 []
  set juv4 []
  set juv5 []
  set juv6 []
  set juv7 []
  set juv8 []
  set juv9 []
  set abs_juv []

set move0 []
set move1 []
set move2 []
set move3 []
set move4 []
set move5 []
set move6 []
set move7 []
set move8 []
set move9 []

set field_max_preg0 []
set field_max_preg1 []
set field_max_preg2 []
set field_max_preg3 []
set field_max_preg4 []
set field_max_preg5 []
set field_max_preg6 []
set field_max_preg7 []
set field_max_preg8 []
set field_max_preg9 []

  if stages_output
  [

  set field_young0 []
set field_young1 []
set field_young2 []
set field_young3 []
set field_young4 []
set field_young5 []
set field_young6 []
set field_young7 []
set field_young8 []
set field_young9 []

set field_mal0 []
set field_mal1 []
set field_mal2 []
set field_mal3 []
set field_mal4 []
set field_mal5 []
set field_mal6 []
set field_mal7 []
set field_mal8 []
set field_mal9 []

set field_preg0 []
set field_preg1 []
set field_preg2 []
set field_preg3 []
set field_preg4 []
set field_preg5 []
set field_preg6 []
set field_preg7 []
set field_preg8 []
set field_preg9 []


set bas_young0 []
set bas_young1 []
set bas_young2 []
set bas_young3 []
set bas_young4 []
set bas_young5 []
set bas_young6 []
set bas_young7 []
set bas_young8 []
set bas_young9 []

set bas_mal0 []
set bas_mal1 []
set bas_mal2 []
set bas_mal3 []
set bas_mal4 []
set bas_mal5 []
set bas_mal6 []
set bas_mal7 []
set bas_mal8 []
set bas_mal9 []

set bas_preg0 []
set bas_preg1 []
set bas_preg2 []
set bas_preg3 []
set bas_preg4 []
set bas_preg5 []
set bas_preg6 []
set bas_preg7 []
set bas_preg8 []
set bas_preg9 []

set loco_young0 []
set loco_young1 []
set loco_young2 []
set loco_young3 []
set loco_young4 []
set loco_young5 []
set loco_young6 []
set loco_young7 []
set loco_young8 []
set loco_young9 []

set loco_mal0 []
set loco_mal1 []
set loco_mal2 []
set loco_mal3 []
set loco_mal4 []
set loco_mal5 []
set loco_mal6 []
set loco_mal7 []
set loco_mal8 []
set loco_mal9 []

set loco_preg0 []
set loco_preg1 []
set loco_preg2 []
set loco_preg3 []
set loco_preg4 []
set loco_preg5 []
set loco_preg6 []
set loco_preg7 []
set loco_preg8 []
set loco_preg9 []

set repro_young0 []
set repro_young1 []
set repro_young2 []
set repro_young3 []
set repro_young4 []
set repro_young5 []
set repro_young6 []
set repro_young7 []
set repro_young8 []
set repro_young9 []

set repro_mal0 []
set repro_mal1 []
set repro_mal2 []
set repro_mal3 []
set repro_mal4 []
set repro_mal5 []
set repro_mal6 []
set repro_mal7 []
set repro_mal8 []
set repro_mal9 []

set repro_preg0 []
set repro_preg1 []
set repro_preg2 []
set repro_preg3 []
set repro_preg4 []
set repro_preg5 []
set repro_preg6 []
set repro_preg7 []
set repro_preg8 []
set repro_preg9 []

set grow_young0 []
set grow_young1 []
set grow_young2 []
set grow_young3 []
set grow_young4 []
set grow_young5 []
set grow_young6 []
set grow_young7 []
set grow_young8 []
set grow_young9 []

set grow_mal0 []
set grow_mal1 []
set grow_mal2 []
set grow_mal3 []
set grow_mal4 []
set grow_mal5 []
set grow_mal6 []
set grow_mal7 []
set grow_mal8 []
set grow_mal9 []

set grow_preg0 []
set grow_preg1 []
set grow_preg2 []
set grow_preg3 []
set grow_preg4 []
set grow_preg5 []
set grow_preg6 []
set grow_preg7 []
set grow_preg8 []
set grow_preg9 []

set dig_young0 []
set dig_young1 []
set dig_young2 []
set dig_young3 []
set dig_young4 []
set dig_young5 []
set dig_young6 []
set dig_young7 []
set dig_young8 []
set dig_young9 []

set dig_mal0 []
set dig_mal1 []
set dig_mal2 []
set dig_mal3 []
set dig_mal4 []
set dig_mal5 []
set dig_mal6 []
set dig_mal7 []
set dig_mal8 []
set dig_mal9 []

set dig_preg0 []
set dig_preg1 []
set dig_preg2 []
set dig_preg3 []
set dig_preg4 []
set dig_preg5 []
set dig_preg6 []
set dig_preg7 []
set dig_preg8 []
set dig_preg9 []
  ]
set life_juv0 []
set life_juv1 []
set life_juv2 []
set life_juv3 []
set life_juv4 []
set life_juv5 []
set life_juv6 []
set life_juv7 []
set life_juv8 []
set life_juv9 []

set life_juv_hr0 []
set life_juv_hr1 []
set life_juv_hr2 []
set life_juv_hr3 []
set life_juv_hr4 []
set life_juv_hr5 []
set life_juv_hr6 []
set life_juv_hr7 []
set life_juv_hr8 []
set life_juv_hr9 []

  set mort_spec []
  set mort_order []
  set mort_stor []
  set mort_food []
  set mort_hr []

  if long_output
  [
set field_spec0 []
set field_spec1 []
set field_spec2 []
set field_spec3 []
set field_spec4 []
set field_spec5 []
set field_spec6 []
set field_spec7 []
set field_spec8 []
set field_spec9 []

  set bas_spec0 []
set bas_spec1 []
set bas_spec2 []
set bas_spec3 []
set bas_spec4 []
set bas_spec5 []
set bas_spec6 []
set bas_spec7 []
set bas_spec8 []
set bas_spec9 []

  set grow_spec0 []
set grow_spec1 []
set grow_spec2 []
set grow_spec3 []
set grow_spec4 []
set grow_spec5 []
set grow_spec6 []
set grow_spec7 []
set grow_spec8 []
set grow_spec9 []

  set repro_spec0 []
set repro_spec1 []
set repro_spec2 []
set repro_spec3 []
set repro_spec4 []
set repro_spec5 []
set repro_spec6 []
set repro_spec7 []
set repro_spec8 []
set repro_spec9 []

  set gest_spec0 []
set gest_spec1 []
set gest_spec2 []
set gest_spec3 []
set gest_spec4 []
set gest_spec5 []
set gest_spec6 []
set gest_spec7 []
set gest_spec8 []
set gest_spec9 []

  set lact_spec0 []
set lact_spec1 []
set lact_spec2 []
set lact_spec3 []
set lact_spec4 []
set lact_spec5 []
set lact_spec6 []
set lact_spec7 []
set lact_spec8 []
set lact_spec9 []

  set loco_spec0 []
set loco_spec1 []
set loco_spec2 []
set loco_spec3 []
set loco_spec4 []
set loco_spec5 []
set loco_spec6 []
set loco_spec7 []
set loco_spec8 []
set loco_spec9 []

    set digest_spec0 []
set digest_spec1 []
set digest_spec2 []
set digest_spec3 []
set digest_spec4 []
set digest_spec5 []
set digest_spec6 []
set digest_spec7 []
set digest_spec8 []
set digest_spec9 []

      set digest_mother_spec0 []
set digest_mother_spec1 []
set digest_mother_spec2 []
set digest_mother_spec3 []
set digest_mother_spec4 []
set digest_mother_spec5 []
set digest_mother_spec6 []
set digest_mother_spec7 []
set digest_mother_spec8 []
set digest_mother_spec9 []

  set mass_spec0 []
    set mass_spec1 []
    set mass_spec2 []
    set mass_spec3 []
    set mass_spec4 []
    set mass_spec5 []
    set mass_spec6 []
    set mass_spec7 []
    set mass_spec8 []
    set mass_spec9 []

    set young_mass_spec0 []
    set young_mass_spec1 []
    set young_mass_spec2 []
    set young_mass_spec3 []
    set young_mass_spec4 []
    set young_mass_spec5 []
    set young_mass_spec6 []
    set young_mass_spec7 []
    set young_mass_spec8 []
    set young_mass_spec9 []

    set stor_spec0 []
    set stor_spec1 []
    set stor_spec2 []
    set stor_spec3 []
    set stor_spec4 []
    set stor_spec5 []
    set stor_spec6 []
    set stor_spec7 []
    set stor_spec8 []
    set stor_spec9 []

  repeat 1766.53 * (maxmass ^ 0.21)
  [
    set field_spec0 lput [] field_spec0
    set field_spec1 lput [] field_spec1
    set field_spec2 lput [] field_spec2
    set field_spec3 lput [] field_spec3
    set field_spec4 lput [] field_spec4
    set field_spec5 lput [] field_spec5
    set field_spec6 lput [] field_spec6
    set field_spec7 lput [] field_spec7
    set field_spec8 lput [] field_spec8
    set field_spec9 lput [] field_spec9

    set bas_spec0 lput [] bas_spec0
    set bas_spec1 lput [] bas_spec1
    set bas_spec2 lput [] bas_spec2
    set bas_spec3 lput [] bas_spec3
    set bas_spec4 lput [] bas_spec4
    set bas_spec5 lput [] bas_spec5
    set bas_spec6 lput [] bas_spec6
    set bas_spec7 lput [] bas_spec7
    set bas_spec8 lput [] bas_spec8
    set bas_spec9 lput [] bas_spec9

    set grow_spec0 lput [] grow_spec0
    set grow_spec1 lput [] grow_spec1
    set grow_spec2 lput [] grow_spec2
    set grow_spec3 lput [] grow_spec3
    set grow_spec4 lput [] grow_spec4
    set grow_spec5 lput [] grow_spec5
    set grow_spec6 lput [] grow_spec6
    set grow_spec7 lput [] grow_spec7
    set grow_spec8 lput [] grow_spec8
    set grow_spec9 lput [] grow_spec9

    set repro_spec0 lput [] repro_spec0
    set repro_spec1 lput [] repro_spec1
    set repro_spec2 lput [] repro_spec2
    set repro_spec3 lput [] repro_spec3
    set repro_spec4 lput [] repro_spec4
    set repro_spec5 lput [] repro_spec5
    set repro_spec6 lput [] repro_spec6
    set repro_spec7 lput [] repro_spec7
    set repro_spec8 lput [] repro_spec8
    set repro_spec9 lput [] repro_spec9

    set gest_spec0 lput [] gest_spec0
    set gest_spec1 lput [] gest_spec1
    set gest_spec2 lput [] gest_spec2
    set gest_spec3 lput [] gest_spec3
    set gest_spec4 lput [] gest_spec4
    set gest_spec5 lput [] gest_spec5
    set gest_spec6 lput [] gest_spec6
    set gest_spec7 lput [] gest_spec7
    set gest_spec8 lput [] gest_spec8
    set gest_spec9 lput [] gest_spec9

    set lact_spec0 lput [] lact_spec0
    set lact_spec1 lput [] lact_spec1
    set lact_spec2 lput [] lact_spec2
    set lact_spec3 lput [] lact_spec3
    set lact_spec4 lput [] lact_spec4
    set lact_spec5 lput [] lact_spec5
    set lact_spec6 lput [] lact_spec6
    set lact_spec7 lput [] lact_spec7
    set lact_spec8 lput [] lact_spec8
    set lact_spec9 lput [] lact_spec9

    set loco_spec0 lput [] loco_spec0
    set loco_spec1 lput [] loco_spec1
    set loco_spec2 lput [] loco_spec2
    set loco_spec3 lput [] loco_spec3
    set loco_spec4 lput [] loco_spec4
    set loco_spec5 lput [] loco_spec5
    set loco_spec6 lput [] loco_spec6
    set loco_spec7 lput [] loco_spec7
    set loco_spec8 lput [] loco_spec8
    set loco_spec9 lput [] loco_spec9

    set digest_spec0 lput [] digest_spec0
    set digest_spec1 lput [] digest_spec1
    set digest_spec2 lput [] digest_spec2
    set digest_spec3 lput [] digest_spec3
    set digest_spec4 lput [] digest_spec4
    set digest_spec5 lput [] digest_spec5
    set digest_spec6 lput [] digest_spec6
    set digest_spec7 lput [] digest_spec7
    set digest_spec8 lput [] digest_spec8
    set digest_spec9 lput [] digest_spec9

    set mass_spec0 lput [] mass_spec0
    set mass_spec1 lput [] mass_spec1
    set mass_spec2 lput [] mass_spec2
    set mass_spec3 lput [] mass_spec3
    set mass_spec4 lput [] mass_spec4
    set mass_spec5 lput [] mass_spec5
    set mass_spec6 lput [] mass_spec6
    set mass_spec7 lput [] mass_spec7
    set mass_spec8 lput [] mass_spec8
    set mass_spec9 lput [] mass_spec9

    set young_mass_spec0 lput [] young_mass_spec0
    set young_mass_spec1 lput [] young_mass_spec1
    set young_mass_spec2 lput [] young_mass_spec2
    set young_mass_spec3 lput [] young_mass_spec3
    set young_mass_spec4 lput [] young_mass_spec4
    set young_mass_spec5 lput [] young_mass_spec5
    set young_mass_spec6 lput [] young_mass_spec6
    set young_mass_spec7 lput [] young_mass_spec7
    set young_mass_spec8 lput [] young_mass_spec8
    set young_mass_spec9 lput [] young_mass_spec9

    set stor_spec0 lput [] stor_spec0
    set stor_spec1 lput [] stor_spec1
    set stor_spec2 lput [] stor_spec2
    set stor_spec3 lput [] stor_spec3
    set stor_spec4 lput [] stor_spec4
    set stor_spec5 lput [] stor_spec5
    set stor_spec6 lput [] stor_spec6
    set stor_spec7 lput [] stor_spec7
    set stor_spec8 lput [] stor_spec8
    set stor_spec9 lput [] stor_spec9
  ]
    repeat 75
  [
    set digest_mother_spec0 lput [] digest_mother_spec0
    set digest_mother_spec1 lput [] digest_mother_spec1
    set digest_mother_spec2 lput [] digest_mother_spec2
    set digest_mother_spec3 lput [] digest_mother_spec3
    set digest_mother_spec4 lput [] digest_mother_spec4
    set digest_mother_spec5 lput [] digest_mother_spec5
    set digest_mother_spec6 lput [] digest_mother_spec6
    set digest_mother_spec7 lput [] digest_mother_spec7
    set digest_mother_spec8 lput [] digest_mother_spec8
    set digest_mother_spec9 lput [] digest_mother_spec9
  ]
  ]

    set fmr_pup_gest_spec0 []
    set fmr_pup_gest_spec1 []
    set fmr_pup_gest_spec2 []
    set fmr_pup_gest_spec3 []
    set fmr_pup_gest_spec4 []
    set fmr_pup_gest_spec5 []
    set fmr_pup_gest_spec6 []
    set fmr_pup_gest_spec7 []
    set fmr_pup_gest_spec8 []
    set fmr_pup_gest_spec9 []

    set fmr_pup_lact_spec0 []
    set fmr_pup_lact_spec1 []
    set fmr_pup_lact_spec2 []
    set fmr_pup_lact_spec3 []
    set fmr_pup_lact_spec4 []
    set fmr_pup_lact_spec5 []
    set fmr_pup_lact_spec6 []
    set fmr_pup_lact_spec7 []
    set fmr_pup_lact_spec8 []
    set fmr_pup_lact_spec9 []

    repeat 5
    [
      set fmr_pup_gest_spec0 lput [] fmr_pup_gest_spec0
      set fmr_pup_gest_spec1 lput [] fmr_pup_gest_spec1
      set fmr_pup_gest_spec2 lput [] fmr_pup_gest_spec2
      set fmr_pup_gest_spec3 lput [] fmr_pup_gest_spec3
      set fmr_pup_gest_spec4 lput [] fmr_pup_gest_spec4
      set fmr_pup_gest_spec5 lput [] fmr_pup_gest_spec5
      set fmr_pup_gest_spec6 lput [] fmr_pup_gest_spec6
      set fmr_pup_gest_spec7 lput [] fmr_pup_gest_spec7
      set fmr_pup_gest_spec8 lput [] fmr_pup_gest_spec8
      set fmr_pup_gest_spec9 lput [] fmr_pup_gest_spec9

      set fmr_pup_lact_spec0 lput [] fmr_pup_lact_spec0
      set fmr_pup_lact_spec1 lput [] fmr_pup_lact_spec1
      set fmr_pup_lact_spec2 lput [] fmr_pup_lact_spec2
      set fmr_pup_lact_spec3 lput [] fmr_pup_lact_spec3
      set fmr_pup_lact_spec4 lput [] fmr_pup_lact_spec4
      set fmr_pup_lact_spec5 lput [] fmr_pup_lact_spec5
      set fmr_pup_lact_spec6 lput [] fmr_pup_lact_spec6
      set fmr_pup_lact_spec7 lput [] fmr_pup_lact_spec7
      set fmr_pup_lact_spec8 lput [] fmr_pup_lact_spec8
      set fmr_pup_lact_spec9 lput [] fmr_pup_lact_spec9
    ]

  repeat 10
  [
    set mort_order lput [] mort_order
    set mort_stor lput [] mort_stor
    set mort_food lput [] mort_food
    set mort_hr lput [] mort_hr
  ]
  ;]
end

;------------------------------------------------------------------------------------------------------------------------
;output in form of lists with daily values

to out

set number lput count turtles number
set number0 lput count turtles with [species = 0] number0
set number1 lput count turtles with [species = 1] number1
set number2 lput count turtles with [species = 2] number2
set number3 lput count turtles with [species = 3] number3
set number4 lput count turtles with [species = 4] number4
set number5 lput count turtles with [species = 5] number5
set number6 lput count turtles with [species = 6] number6
set number7 lput count turtles with [species = 7] number7
set number8 lput count turtles with [species = 8] number8
set number9 lput count turtles with [species = 9] number9

set juv_fail lput fail juv_fail
set feed_fail lput fail3 feed_fail

let spec 0
repeat 10
  [
  set suc_immi replace-item spec suc_immi (item spec suc_immi + count turtles with [species = spec and just-immi = 1])
  set prod_sum replace-item spec prod_sum (item spec prod_sum + sum [production] of turtles with [species = spec])
  set spec spec + 1
  ]

set feed_fail_spec0 lput item 0 fail3_spec feed_fail_spec0
set feed_fail_spec1 lput item 1 fail3_spec feed_fail_spec1
set feed_fail_spec2 lput item 2 fail3_spec feed_fail_spec2
set feed_fail_spec3 lput item 3 fail3_spec feed_fail_spec3
set feed_fail_spec4 lput item 4 fail3_spec feed_fail_spec4
set feed_fail_spec5 lput item 5 fail3_spec feed_fail_spec5
set feed_fail_spec6 lput item 6 fail3_spec feed_fail_spec6
set feed_fail_spec7 lput item 7 fail3_spec feed_fail_spec7
set feed_fail_spec8 lput item 8 fail3_spec feed_fail_spec8
set feed_fail_spec9 lput item 9 fail3_spec feed_fail_spec9

  if any? turtles with [age > 1]
  [
set hr lput mean [hrsize] of turtles with [age > 1] hr
set maxhrs lput mean [hrsize / maxhr] of turtles with [age > 1] maxhrs
set countmaxhr lput count turtles with [age > 1 and hrsize / maxhr > 0.9] countmaxhr
set stors lput mean [storage / max_storage] of turtles with [age > 1] stors
set ages lput mean [age] of turtles with [age > 1] ages
  ]

  if any? turtles with [species = 0 and age > 1]
  [
set hr0 lput mean [hrsize] of turtles with [species = 0 and age > 1] hr0
set maxhr0 lput mean [hrsize / maxhr] of turtles with [species = 0 and age > 1] maxhr0
set countmaxhr0 lput count turtles with [species = 0 and age > 1 and hrsize / maxhr > 0.9] countmaxhr0
set age0 lput mean [age] of turtles with [species = 0 and age > 1] age0
set stor0 lput mean [storage / max_storage] of turtles with [species = 0 and age > 1] stor0
set move0 lput mean [daymoved] of turtles with [species = 0 and age > 1] move0
set patches0 lput mean [forage-patches] of turtles with [species = 0 and age > 1] patches0
set compet0 lput mean [competition / forage-patches] of turtles with [species = 0 and age > 1] compet0
set order0 lput mean [order] of turtles with [species = 0 and age > 1] order0
set vis_before0 lput mean [vis_before / forage-patches] of turtles with [species = 0 and age > 1] vis_before0
set fmr0 lput mean [fmr] of turtles with [species = 0 and age > 1] fmr0
set loco0 lput mean [fmr_loco / fmr] of turtles with [species = 0 and age > 1] loco0
set repro0 lput mean [fmr_repro / fmr] of turtles with [species = 0 and age > 1] repro0
set grow0 lput mean [fmr_growth / fmr] of turtles with [species = 0 and age > 1] grow0
set basal0 lput mean [fmr_basal / fmr] of turtles with [species = 0 and age > 1] basal0
set digest0 lput mean [fmr_digest / fmr] of turtles with [species = 0 and age > 1] digest0
set prod0 lput mean [production] of turtles with [species = 0 and age > 1] prod0
set in0 lput mean [intake] of turtles with [species = 0 and age > 1] in0
set inpatch0 lput mean [intake / forage-patches] of turtles with [species = 0 and age > 1] inpatch0
set inhr0 lput mean [intake / (hrsize + 1)] of turtles with [species = 0 and age > 1] inhr0
set balance0 lput median [balance] of turtles with [species = 0 and age > 1] balance0
set prodbalance0 lput median [balanceprod] of turtles with [species = 0 and age > 1] prodbalance0
  ]
  if any? turtles with [species = 0 and preg = 0]
  [
set abs_stor0 lput mean [storage] of turtles with [species = 0 and preg = 0] abs_stor0
  ]

  if any? turtles with [species = 1 and age > 1]
  [
set hr1 lput mean [hrsize] of turtles with [species = 1 and age > 1] hr1
set maxhr1 lput mean [hrsize / maxhr] of turtles with [species = 1 and age > 1] maxhr1
set countmaxhr1 lput count turtles with [species = 1 and age > 1 and hrsize / maxhr > 0.9] countmaxhr1
set age1 lput mean [age] of turtles with [species = 1 and age > 1] age1
set stor1 lput mean [storage / max_storage] of turtles with [species = 1 and age > 1] stor1
set move1 lput mean [daymoved] of turtles with [species = 1 and age > 1] move1
set patches1 lput mean [forage-patches] of turtles with [species = 1 and age > 1] patches1
set compet1 lput mean [competition / forage-patches] of turtles with [species = 1 and age > 1] compet1
set order1 lput mean [order] of turtles with [species = 1 and age > 1] order1
set vis_before1 lput mean [vis_before / forage-patches] of turtles with [species = 1 and age > 1] vis_before1
set fmr1 lput mean [fmr] of turtles with [species = 1 and age > 1] fmr1
set loco1 lput mean [fmr_loco / fmr] of turtles with [species = 1 and age > 1] loco1
set repro1 lput mean [fmr_repro / fmr] of turtles with [species = 1 and age > 1] repro1
set grow1 lput mean [fmr_growth / fmr] of turtles with [species = 1 and age > 1] grow1
set basal1 lput mean [fmr_basal / fmr] of turtles with [species = 1 and age > 1] basal1
set digest1 lput mean [fmr_digest / fmr] of turtles with [species = 1 and age > 1] digest1
set prod1 lput mean [production] of turtles with [species = 1 and age > 1] prod1
set in1 lput mean [intake] of turtles with [species = 1 and age > 1] in1
set inpatch1 lput mean [intake / forage-patches] of turtles with [species = 1 and age > 1] inpatch1
set inhr1 lput mean [intake / (hrsize + 1)] of turtles with [species = 1 and age > 1] inhr1
set balance1 lput median [balance] of turtles with [species = 1 and age > 1] balance1
set prodbalance1 lput median [balanceprod] of turtles with [species = 1 and age > 1] prodbalance1
  ]
  if any? turtles with [species = 1 and preg = 0]
  [
set abs_stor1 lput mean [storage] of turtles with [species = 1 and preg = 0] abs_stor1
  ]

  if any? turtles with [species = 2 and age > 1]
  [
set hr2 lput mean [hrsize] of turtles with [species = 2 and age > 1] hr2
set maxhr2 lput mean [hrsize / maxhr] of turtles with [species = 2 and age > 1] maxhr2
set countmaxhr2 lput count turtles with [species = 2 and age > 1 and hrsize / maxhr > 0.9] countmaxhr2
set age2 lput mean [age] of turtles with [species = 2 and age > 1] age2
set stor2 lput mean [storage / max_storage] of turtles with [species = 2 and age > 1] stor2
set move2 lput mean [daymoved] of turtles with [species = 2 and age > 1] move2
set patches2 lput mean [forage-patches] of turtles with [species = 2 and age > 1] patches2
set compet2 lput mean [competition / forage-patches] of turtles with [species = 2 and age > 1] compet2
set order2 lput mean [order] of turtles with [species = 2 and age > 1] order2
set vis_before2 lput mean [vis_before / forage-patches] of turtles with [species = 2 and age > 1] vis_before2
set fmr2 lput mean [fmr] of turtles with [species = 2 and age > 1] fmr2
set loco2 lput mean [fmr_loco / fmr] of turtles with [species = 2 and age > 1] loco2
set repro2 lput mean [fmr_repro / fmr] of turtles with [species = 2 and age > 1] repro2
set grow2 lput mean [fmr_growth / fmr] of turtles with [species = 2 and age > 1] grow2
set basal2 lput mean [fmr_basal / fmr] of turtles with [species = 2 and age > 1] basal2
set digest2 lput mean [fmr_digest / fmr] of turtles with [species = 2 and age > 1] digest2
set prod2 lput mean [production] of turtles with [species = 2 and age > 1] prod2
set in2 lput mean [intake] of turtles with [species = 2 and age > 1] in2
set inpatch2 lput mean [intake / forage-patches] of turtles with [species = 2 and age > 1] inpatch2
set inhr2 lput mean [intake / (hrsize + 1)] of turtles with [species = 2 and age > 1] inhr2
set balance2 lput median [balance] of turtles with [species = 2 and age > 1] balance2
set prodbalance2 lput median [balanceprod] of turtles with [species = 2 and age > 1] prodbalance2
  ]
  if any? turtles with [species = 2 and preg = 0]
  [
set abs_stor2 lput mean [storage] of turtles with [species = 2 and preg = 0] abs_stor2
  ]

  if any? turtles with [species = 3 and age > 1]
  [
set hr3 lput mean [hrsize] of turtles with [species = 3 and age > 1] hr3
set maxhr3 lput mean [hrsize / maxhr] of turtles with [species = 3 and age > 1] maxhr3
set countmaxhr3 lput count turtles with [species = 3 and age > 1 and hrsize / maxhr > 0.9] countmaxhr3
set age3 lput mean [age] of turtles with [species = 3 and age > 1] age3
set stor3 lput mean [storage / max_storage] of turtles with [species = 3 and age > 1] stor3
set move3 lput mean [daymoved] of turtles with [species = 3 and age > 1] move3
set patches3 lput mean [forage-patches] of turtles with [species = 3 and age > 1] patches3
set compet3 lput mean [competition / forage-patches] of turtles with [species = 3 and age > 1] compet3
set order3 lput mean [order] of turtles with [species = 3 and age > 1] order3
set vis_before3 lput mean [vis_before / forage-patches] of turtles with [species = 3 and age > 1] vis_before3
set fmr3 lput mean [fmr] of turtles with [species = 3 and age > 1] fmr3
set loco3 lput mean [fmr_loco / fmr] of turtles with [species = 3 and age > 1] loco3
set repro3 lput mean [fmr_repro / fmr] of turtles with [species = 3 and age > 1] repro3
set grow3 lput mean [fmr_growth / fmr] of turtles with [species = 3 and age > 1] grow3
set basal3 lput mean [fmr_basal / fmr] of turtles with [species = 3 and age > 1] basal3
set digest3 lput mean [fmr_digest / fmr] of turtles with [species = 3 and age > 1] digest3
set prod3 lput mean [production] of turtles with [species = 3 and age > 1] prod3
set in3 lput mean [intake] of turtles with [species = 3 and age > 1] in3
set inpatch3 lput mean [intake / forage-patches] of turtles with [species = 3 and age > 1] inpatch3
set inhr3 lput mean [intake / (hrsize + 1)] of turtles with [species = 3 and age > 1] inhr3
set balance3 lput median [balance] of turtles with [species = 3 and age > 1] balance3
set prodbalance3 lput median [balanceprod] of turtles with [species = 3 and age > 1] prodbalance3
  ]
  if any? turtles with [species = 3 and preg = 0]
  [
set abs_stor3 lput mean [storage] of turtles with [species = 3 and preg = 0] abs_stor3
  ]

  if any? turtles with [species = 4 and age > 1]
  [
set hr4 lput mean [hrsize] of turtles with [species = 4 and age > 1] hr4
set maxhr4 lput mean [hrsize / maxhr] of turtles with [species = 4 and age > 1] maxhr4
set countmaxhr4 lput count turtles with [species = 4 and age > 1 and hrsize / maxhr > 0.9] countmaxhr4
set age4 lput mean [age] of turtles with [species = 4 and age > 1] age4
set stor4 lput mean [storage / max_storage] of turtles with [species = 4 and age > 1] stor4
set move4 lput mean [daymoved] of turtles with [species = 4 and age > 1] move4
set patches4 lput mean [forage-patches] of turtles with [species = 4 and age > 1] patches4
set compet4 lput mean [competition / forage-patches] of turtles with [species = 4 and age > 1] compet4
set order4 lput mean [order] of turtles with [species = 4 and age > 1] order4
set vis_before4 lput mean [vis_before / forage-patches] of turtles with [species = 4 and age > 1] vis_before4
set fmr4 lput mean [fmr] of turtles with [species = 4 and age > 1] fmr4
set loco4 lput mean [fmr_loco / fmr] of turtles with [species = 4 and age > 1] loco4
set repro4 lput mean [fmr_repro / fmr] of turtles with [species = 4 and age > 1] repro4
set grow4 lput mean [fmr_growth / fmr] of turtles with [species = 4 and age > 1] grow4
set basal4 lput mean [fmr_basal / fmr] of turtles with [species = 4 and age > 1] basal4
set digest4 lput mean [fmr_digest / fmr] of turtles with [species = 4 and age > 1] digest4
set prod4 lput mean [production] of turtles with [species = 4 and age > 1] prod4
set in4 lput mean [intake] of turtles with [species = 4 and age > 1] in4
set inpatch4 lput mean [intake / forage-patches] of turtles with [species = 4 and age > 1] inpatch4
set inhr4 lput mean [intake / (hrsize + 1)] of turtles with [species = 4 and age > 1] inhr4
set balance4 lput median [balance] of turtles with [species = 4 and age > 1] balance4
set prodbalance4 lput median [balanceprod] of turtles with [species = 4 and age > 1] prodbalance4
  ]
  if any? turtles with [species = 4 and preg = 0]
  [
set abs_stor4 lput mean [storage] of turtles with [species = 4 and preg = 0] abs_stor4
  ]

  if any? turtles with [species = 5 and age > 1]
  [
set hr5 lput mean [hrsize] of turtles with [species = 5 and age > 1] hr5
set maxhr5 lput mean [hrsize / maxhr] of turtles with [species = 5 and age > 1] maxhr5
set countmaxhr5 lput count turtles with [species = 5 and age > 1 and hrsize / maxhr > 0.9] countmaxhr5
set age5 lput mean [age] of turtles with [species = 5 and age > 1] age5
set stor5 lput mean [storage / max_storage] of turtles with [species = 5 and age > 1] stor5
set move5 lput mean [daymoved] of turtles with [species = 5 and age > 1] move5
set patches5 lput mean [forage-patches] of turtles with [species = 5 and age > 1] patches5
set compet5 lput mean [competition / forage-patches] of turtles with [species = 5 and age > 1] compet5
set order5 lput mean [order] of turtles with [species = 5 and age > 1] order5
set vis_before5 lput mean [vis_before / forage-patches] of turtles with [species = 5 and age > 1] vis_before5
set fmr5 lput mean [fmr] of turtles with [species = 5 and age > 1] fmr5
set loco5 lput mean [fmr_loco / fmr] of turtles with [species = 5 and age > 1] loco5
set repro5 lput mean [fmr_repro / fmr] of turtles with [species = 5 and age > 1] repro5
set grow5 lput mean [fmr_growth / fmr] of turtles with [species = 5 and age > 1] grow5
set basal5 lput mean [fmr_basal / fmr] of turtles with [species = 5 and age > 1] basal5
set digest5 lput mean [fmr_digest / fmr] of turtles with [species = 5 and age > 1] digest5
set prod5 lput mean [production] of turtles with [species = 5 and age > 1] prod5
set in5 lput mean [intake] of turtles with [species = 5 and age > 1] in5
set inpatch5 lput mean [intake / forage-patches] of turtles with [species = 5 and age > 1] inpatch5
set inhr5 lput mean [intake / (hrsize + 1)] of turtles with [species = 5 and age > 1] inhr5
set balance5 lput median [balance] of turtles with [species = 5 and age > 1] balance5
set prodbalance5 lput median [balanceprod] of turtles with [species = 5 and age > 1] prodbalance5
  ]
  if any? turtles with [species = 5 and preg = 0]
  [
set abs_stor5 lput mean [storage] of turtles with [species = 5 and preg = 0] abs_stor5
  ]

  if any? turtles with [species = 6 and age > 1]
  [
set hr6 lput mean [hrsize] of turtles with [species = 6 and age > 1] hr6
set maxhr6 lput mean [hrsize / maxhr] of turtles with [species = 6 and age > 1] maxhr6
set countmaxhr6 lput count turtles with [species = 6 and age > 1 and hrsize / maxhr > 0.9] countmaxhr6
set age6 lput mean [age] of turtles with [species = 6 and age > 1] age6
set stor6 lput mean [storage / max_storage] of turtles with [species = 6 and age > 1] stor6
set move6 lput mean [daymoved] of turtles with [species = 6 and age > 1] move6
set patches6 lput mean [forage-patches] of turtles with [species = 6 and age > 1] patches6
set compet6 lput mean [competition / forage-patches] of turtles with [species = 6 and age > 1] compet6
set order6 lput mean [order] of turtles with [species = 6 and age > 1] order6
set vis_before6 lput mean [vis_before / forage-patches] of turtles with [species = 6 and age > 1] vis_before6
set fmr6 lput mean [fmr] of turtles with [species = 6 and age > 1] fmr6
set loco6 lput mean [fmr_loco / fmr] of turtles with [species = 6 and age > 1] loco6
set repro6 lput mean [fmr_repro / fmr] of turtles with [species = 6 and age > 1] repro6
set grow6 lput mean [fmr_growth / fmr] of turtles with [species = 6 and age > 1] grow6
set basal6 lput mean [fmr_basal / fmr] of turtles with [species = 6 and age > 1] basal6
set digest6 lput mean [fmr_digest / fmr] of turtles with [species = 6 and age > 1] digest6
set prod6 lput mean [production] of turtles with [species = 6 and age > 1] prod6
set in6 lput mean [intake] of turtles with [species = 6 and age > 1] in6
set inpatch6 lput mean [intake / forage-patches] of turtles with [species = 6 and age > 1] inpatch6
set inhr6 lput mean [intake / (hrsize + 1)] of turtles with [species = 6 and age > 1] inhr6
set balance6 lput median [balance] of turtles with [species = 6 and age > 1] balance6
set prodbalance6 lput median [balanceprod] of turtles with [species = 6 and age > 1] prodbalance6
  ]
  if any? turtles with [species = 6 and preg = 0]
  [
set abs_stor6 lput mean [storage] of turtles with [species = 6 and preg = 0] abs_stor6
  ]

  if any? turtles with [species = 7 and age > 1]
  [
set hr7 lput mean [hrsize] of turtles with [species = 7 and age > 1] hr7
set maxhr7 lput mean [hrsize / maxhr] of turtles with [species = 7 and age > 1] maxhr7
set countmaxhr7 lput count turtles with [species = 7 and age > 1 and hrsize / maxhr > 0.9] countmaxhr7
set age7 lput mean [age] of turtles with [species = 7 and age > 1] age7
set move7 lput mean [daymoved] of turtles with [species = 7 and age > 1] move7
set stor7 lput mean [storage / max_storage] of turtles with [species = 7 and age > 1] stor7
set patches7 lput mean [forage-patches] of turtles with [species = 7 and age > 1] patches7
set compet7 lput mean [competition / forage-patches] of turtles with [species = 7 and age > 1] compet7
set order7 lput mean [order] of turtles with [species = 7 and age > 1] order7
set vis_before7 lput mean [vis_before / forage-patches] of turtles with [species = 7 and age > 1] vis_before7
set fmr7 lput mean [fmr] of turtles with [species = 7 and age > 1] fmr7
set loco7 lput mean [fmr_loco / fmr] of turtles with [species = 7 and age > 1] loco7
set repro7 lput mean [fmr_repro / fmr] of turtles with [species = 7 and age > 1] repro7
set grow7 lput mean [fmr_growth / fmr] of turtles with [species = 7 and age > 1] grow7
set basal7 lput mean [fmr_basal / fmr] of turtles with [species = 7 and age > 1] basal7
set digest7 lput mean [fmr_digest / fmr] of turtles with [species = 7 and age > 1] digest7
set prod7 lput mean [production] of turtles with [species = 7 and age > 1] prod7
set in7 lput mean [intake] of turtles with [species = 7 and age > 1] in7
set inpatch7 lput mean [intake / forage-patches] of turtles with [species = 7 and age > 1] inpatch7
set inhr7 lput mean [intake / (hrsize + 1)] of turtles with [species = 7 and age > 1] inhr7
set balance7 lput median [balance] of turtles with [species = 7 and age > 1] balance7
set prodbalance7 lput median [balanceprod] of turtles with [species = 7 and age > 1] prodbalance7
  ]
  if any? turtles with [species = 7 and preg = 0]
  [
set abs_stor7 lput mean [storage] of turtles with [species = 7 and preg = 0] abs_stor7
  ]

  if any? turtles with [species = 8 and age > 1]
  [
set hr8 lput mean [hrsize] of turtles with [species = 8 and age > 1] hr8
set maxhr8 lput mean [hrsize / maxhr] of turtles with [species = 8 and age > 1] maxhr8
set countmaxhr8 lput count turtles with [species = 8 and age > 1 and hrsize / maxhr > 0.9] countmaxhr8
set age8 lput mean [age] of turtles with [species = 8 and age > 1] age8
set stor8 lput mean [storage / max_storage] of turtles with [species = 8 and age > 1] stor8
set move8 lput mean [daymoved] of turtles with [species = 8 and age > 1] move8
set patches8 lput mean [forage-patches] of turtles with [species = 8 and age > 1] patches8
set compet8 lput mean [competition / forage-patches] of turtles with [species = 8 and age > 1] compet8
set order8 lput mean [order] of turtles with [species = 8 and age > 1] order8
set vis_before8 lput mean [vis_before / forage-patches] of turtles with [species = 8 and age > 1] vis_before8
set fmr8 lput mean [fmr] of turtles with [species = 8 and age > 1] fmr8
set loco8 lput mean [fmr_loco / fmr] of turtles with [species = 8 and age > 1] loco8
set repro8 lput mean [fmr_repro / fmr] of turtles with [species = 8 and age > 1] repro8
set grow8 lput mean [fmr_growth / fmr] of turtles with [species = 8 and age > 1] grow8
set basal8 lput mean [fmr_basal / fmr] of turtles with [species = 8 and age > 1] basal8
set digest8 lput mean [fmr_digest / fmr] of turtles with [species = 8 and age > 1] digest8
set prod8 lput mean [production] of turtles with [species = 8 and age > 1] prod8
set in8 lput mean [intake] of turtles with [species = 8 and age > 1] in8
set inpatch8 lput mean [intake / forage-patches] of turtles with [species = 8 and age > 1] inpatch8
set inhr8 lput mean [intake / (hrsize + 1)] of turtles with [species = 8 and age > 1] inhr8
set balance8 lput median [balance] of turtles with [species = 8 and age > 1] balance8
set prodbalance8 lput median [balanceprod] of turtles with [species = 8 and age > 1] prodbalance8
  ]
  if any? turtles with [species = 8 and preg = 0]
  [
set abs_stor8 lput mean [storage] of turtles with [species = 8 and preg = 0] abs_stor8
  ]

  if any? turtles with [species = 9 and age > 1]
  [
set hr9 lput mean [hrsize] of turtles with [species = 9 and age > 1] hr9
set maxhr9 lput mean [hrsize / maxhr] of turtles with [species = 9 and age > 1] maxhr9
set countmaxhr9 lput count turtles with [species = 9 and age > 1 and hrsize / maxhr > 0.9] countmaxhr9
set age9 lput mean [age] of turtles with [species = 9 and age > 1] age9
set stor9 lput mean [storage / max_storage] of turtles with [species = 9 and age > 1] stor9
set move9 lput mean [daymoved] of turtles with [species = 9 and age > 1] move9
set patches9 lput mean [forage-patches] of turtles with [species = 9 and age > 1] patches9
set compet9 lput mean [competition / forage-patches] of turtles with [species = 9 and age > 1] compet9
set order9 lput mean [order] of turtles with [species = 9 and age > 1] order9
set vis_before9 lput mean [vis_before / forage-patches] of turtles with [species = 9 and age > 1] vis_before9
set fmr9 lput mean [fmr] of turtles with [species = 9 and age > 1] fmr9
set loco9 lput mean [fmr_loco / fmr] of turtles with [species = 9 and age > 1] loco9
set repro9 lput mean [fmr_repro / fmr] of turtles with [species = 9 and age > 1] repro9
set grow9 lput mean [fmr_growth / fmr] of turtles with [species = 9 and age > 1] grow9
set basal9 lput mean [fmr_basal / fmr] of turtles with [species = 9 and age > 1] basal9
set digest9 lput mean [fmr_digest / fmr] of turtles with [species = 9 and age > 1] digest9
set prod9 lput mean [production] of turtles with [species = 9 and age > 1] prod9
set in9 lput mean [intake] of turtles with [species = 9 and age > 1] in9
set inpatch9 lput mean [intake / forage-patches] of turtles with [species = 9 and age > 1] inpatch9
set inhr9 lput mean [intake / (hrsize + 1)] of turtles with [species = 9 and age > 1] inhr9
set balance9 lput median [balance] of turtles with [species = 9 and age > 1] balance9
set prodbalance9 lput median [balanceprod] of turtles with [species = 9 and age > 1] prodbalance9
  ]
  if any? turtles with [species = 9 and preg = 0]
  [
set abs_stor9 lput mean [storage] of turtles with [species = 9 and preg = 0] abs_stor9
  ]

  if any? turtles with [age > age_first_repro]
  [
set rm lput mean [real-mass] of turtles with [age > age_first_repro] rm
    if any? turtles with [species = 0 and age > age_first_repro]
  [
set rm0 lput mean [real-mass] of turtles with [species = 0 and age > age_first_repro] rm0
  ]
    if any? turtles with [species = 1 and age > age_first_repro]
  [
set rm1 lput mean [real-mass] of turtles with [species = 1 and age > age_first_repro] rm1
  ]
    if any? turtles with [species = 2 and age > age_first_repro]
  [
set rm2 lput mean [real-mass] of turtles with [species = 2 and age > age_first_repro] rm2
  ]
    if any? turtles with [species = 3 and age > age_first_repro]
  [
set rm3 lput mean [real-mass] of turtles with [species = 3 and age > age_first_repro] rm3
  ]
    if any? turtles with [species = 4 and age > age_first_repro]
  [
set rm4 lput mean [real-mass] of turtles with [species = 4 and age > age_first_repro] rm4
  ]
    if any? turtles with [species = 5 and age > age_first_repro]
  [
set rm5 lput mean [real-mass] of turtles with [species = 5 and age > age_first_repro] rm5
  ]
    if any? turtles with [species = 6 and age > age_first_repro]
  [
set rm6 lput mean [real-mass] of turtles with [species = 6 and age > age_first_repro] rm6
  ]
    if any? turtles with [species = 7 and age > age_first_repro]
  [
set rm7 lput mean [real-mass] of turtles with [species = 7 and age > age_first_repro] rm7
  ]
    if any? turtles with [species = 8 and age > age_first_repro]
  [
set rm8 lput mean [real-mass] of turtles with [species = 8 and age > age_first_repro] rm8
  ]
    if any? turtles with [species = 9 and age > age_first_repro]
  [
set rm9 lput mean [real-mass] of turtles with [species = 9 and age > age_first_repro] rm9
]
]

  if any? turtles with [preg = 1]
  [
set hrpreg lput mean [hrsize] of turtles with [preg = 1] hrpreg
set storpreg lput mean [storage / max_storage] of turtles with [preg = 1] storpreg
    if any? turtles with [species = 0 and preg = 1]
    [
set hrpreg0 lput mean [hrsize] of turtles with [species = 0 and preg = 1] hrpreg0
set storpreg0 lput mean [storage / max_storage] of turtles with [species = 0 and preg = 1] storpreg0
    ]
    if any? turtles with [species = 1 and preg = 1]
    [
set hrpreg1 lput mean [hrsize] of turtles with [species = 1 and preg = 1] hrpreg1
set storpreg1 lput mean [storage / max_storage] of turtles with [species = 1 and preg = 1] storpreg1
      ]
    if any? turtles with [species = 2 and preg = 1]
    [
set hrpreg2 lput mean [hrsize] of turtles with [species = 2 and preg = 1] hrpreg2
set storpreg2 lput mean [storage / max_storage] of turtles with [species = 2 and preg = 1] storpreg2
      ]
    if any? turtles with [species = 3 and preg = 1]
    [
set hrpreg3 lput mean [hrsize] of turtles with [species = 3 and preg = 1] hrpreg3
set storpreg3 lput mean [storage / max_storage] of turtles with [species = 3 and preg = 1] storpreg3
      ]
    if any? turtles with [species = 4 and preg = 1]
    [
set hrpreg4 lput mean [hrsize] of turtles with [species = 4 and preg = 1] hrpreg4
set storpreg4 lput mean [storage / max_storage] of turtles with [species = 4 and preg = 1] storpreg4
      ]
    if any? turtles with [species = 5 and preg = 1]
    [
set hrpreg5 lput mean [hrsize] of turtles with [species = 5 and preg = 1] hrpreg5
set storpreg5 lput mean [storage / max_storage] of turtles with [species = 5 and preg = 1] storpreg5
      ]
    if any? turtles with [species = 6 and preg = 1]
    [
set hrpreg6 lput mean [hrsize] of turtles with [species = 6 and preg = 1] hrpreg6
set storpreg6 lput mean [storage / max_storage] of turtles with [species = 6 and preg = 1] storpreg6
      ]
    if any? turtles with [species = 7 and preg = 1]
    [
set hrpreg7 lput mean [hrsize] of turtles with [species = 7 and preg = 1] hrpreg7
set storpreg7 lput mean [storage / max_storage] of turtles with [species = 7 and preg = 1] storpreg7
      ]
    if any? turtles with [species = 8 and preg = 1]
    [
set hrpreg8 lput mean [hrsize] of turtles with [species = 8 and preg = 1] hrpreg8
set storpreg8 lput mean [storage / max_storage] of turtles with [species = 8 and preg = 1] storpreg8
      ]
    if any? turtles with [species = 9 and preg = 1]
    [
set hrpreg9 lput mean [hrsize] of turtles with [species = 9 and preg = 1] hrpreg9
set storpreg9 lput mean [storage / max_storage] of turtles with [species = 9 and preg = 1] storpreg9
    ]
    if any? turtles with [preg = 1 and pregcount = 0]
    [
set stor_pregstart lput mean [storage / max_storage] of turtles with [preg = 1 and pregcount = 0] stor_pregstart
    ]
    if any? turtles with [preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set stor_pregend lput mean [storage / max_storage] of turtles with [preg = 1 and pregcount = round(gest_period + lact_period - 1)] stor_pregend
set juv lput mean [real-young / young] of turtles with [preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv
set abs_juv lput mean [real-young] of turtles with [preg = 1 and pregcount = round(gest_period + lact_period - 1)] abs_juv

      if any? turtles with [species = 0 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv0 lput mean [real-young / young] of turtles with [species = 0 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv0
set field_max_preg0 lput mean [fmr] of turtles with [species = 0 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg0
    ]
      if any? turtles with [species = 1 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv1 lput mean [real-young / young] of turtles with [species = 1 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv1
set field_max_preg1 lput mean [fmr] of turtles with [species = 1 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg1
        ]
      if any? turtles with [species = 2 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv2 lput mean [real-young / young] of turtles with [species = 2 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv2
set field_max_preg2 lput mean [fmr] of turtles with [species = 2 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg2
        ]
      if any? turtles with [species = 3 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv3 lput mean [real-young / young] of turtles with [species = 3 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv3
set field_max_preg3 lput mean [fmr] of turtles with [species = 3 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg3
        ]
      if any? turtles with [species = 4 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv4 lput mean [real-young / young] of turtles with [species = 4 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv4
set field_max_preg4 lput mean [fmr] of turtles with [species = 4 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg4
        ]
      if any? turtles with [species = 5 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv5 lput mean [real-young / young] of turtles with [species = 5 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv5
set field_max_preg5 lput mean [fmr] of turtles with [species = 5 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg5
        ]
      if any? turtles with [species = 6 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv6 lput mean [real-young / young] of turtles with [species = 6 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv6
set field_max_preg6 lput mean [fmr] of turtles with [species = 6 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg6
        ]
      if any? turtles with [species = 7 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv7 lput mean [real-young / young] of turtles with [species = 7 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv7
set field_max_preg7 lput mean [fmr] of turtles with [species = 7 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg7
        ]
      if any? turtles with [species = 8 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv8 lput mean [real-young / young] of turtles with [species = 8 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv8
set field_max_preg8 lput mean [fmr] of turtles with [species = 8 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg8
        ]
      if any? turtles with [species = 9 and preg = 1 and pregcount = round(gest_period + lact_period - 1)]
    [
set juv9 lput mean [real-young / young] of turtles with [species = 9 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] juv9
set field_max_preg9 lput mean [fmr] of turtles with [species = 9 and preg = 1 and pregcount = round(gest_period + lact_period - 1)] field_max_preg9
      ]
  ]
  ]

  if any? turtles with [competition / forage-patches < 5]
  [
set comp5 lput count turtles with [competition / forage-patches < 5] comp5
set comp5_stor lput mean [storage / max_storage] of turtles with [competition / forage-patches < 5] comp5_stor
set comp5_fmr lput mean [fmr / mass] of turtles with [competition / forage-patches < 5] comp5_fmr
  ]
  if any? turtles with [competition / forage-patches >= 5 and competition / forage-patches < 10]
  [
set comp10 lput count turtles with [competition / forage-patches >= 5 and competition / forage-patches < 10] comp10
set comp10_stor lput mean [storage / max_storage] of turtles with [competition / forage-patches >= 5 and competition / forage-patches < 10] comp10_stor
set comp10_fmr lput mean [fmr / mass] of turtles with [competition / forage-patches >= 5 and competition / forage-patches < 10] comp10_fmr
    ]
  if any? turtles with [competition / forage-patches >= 10 and competition / forage-patches < 15]
  [
set comp15 lput count turtles with [competition / forage-patches >= 10 and competition / forage-patches < 15] comp15
set comp15_stor lput mean [storage / max_storage] of turtles with [competition / forage-patches >= 10 and competition / forage-patches < 15] comp15_stor
set comp15_fmr lput mean [fmr / mass] of turtles with [competition / forage-patches >= 10 and competition / forage-patches < 15] comp15_fmr
    ]
  if any? turtles with [competition / forage-patches >= 15 and competition / forage-patches < 20]
  [
set comp20 lput count turtles with [competition / forage-patches >= 15 and competition / forage-patches < 20] comp20
set comp20_stor lput mean [storage / max_storage] of turtles with [competition / forage-patches >= 15 and competition / forage-patches < 20] comp20_stor
set comp20_fmr lput mean [fmr / mass] of turtles with [competition / forage-patches >= 15 and competition / forage-patches < 20] comp20_fmr
    ]
  if any? turtles with [competition / forage-patches >= 20 and competition / forage-patches < 25]
  [
set comp25 lput count turtles with [competition / forage-patches >= 20 and competition / forage-patches < 25] comp25
set comp25_stor lput mean [storage / max_storage] of turtles with [competition / forage-patches >= 20 and competition / forage-patches < 25] comp25_stor
set comp25_fmr lput mean [fmr / mass] of turtles with [competition / forage-patches >= 20 and competition / forage-patches < 25] comp25_fmr
    ]
  if any? turtles with [competition / forage-patches >= 25]
  [
set compmore lput count turtles with [competition / forage-patches >= 25] compmore
set compmore_stor lput mean [storage / max_storage] of turtles with [competition / forage-patches >= 25] compmore_stor
set compmore_fmr lput mean [fmr / mass] of turtles with [competition / forage-patches >= 25] compmore_fmr
  ]

  if length rep_success_0 > 0
  [
set life_juv0 lput mean rep_success_0 life_juv0
set life_juv_hr0 lput mean rep_success_hr_0 life_juv_hr0
  ]
  if length rep_success_1 > 0
  [
set life_juv1 lput mean rep_success_1 life_juv1
set life_juv_hr1 lput mean rep_success_hr_1 life_juv_hr1
    ]
  if length rep_success_2 > 0
  [
set life_juv2 lput mean rep_success_2 life_juv2
set life_juv_hr2 lput mean rep_success_hr_2 life_juv_hr2
    ]
  if length rep_success_3 > 0
  [
set life_juv3 lput mean rep_success_3 life_juv3
set life_juv_hr3 lput mean rep_success_hr_3 life_juv_hr3
    ]
  if length rep_success_4 > 0
  [
set life_juv4 lput mean rep_success_4 life_juv4
set life_juv_hr4 lput mean rep_success_hr_4 life_juv_hr4
    ]
  if length rep_success_5 > 0
  [
set life_juv5 lput mean rep_success_5 life_juv5
set life_juv_hr5 lput mean rep_success_hr_5 life_juv_hr5
    ]
  if length rep_success_6 > 0
  [
set life_juv6 lput mean rep_success_6 life_juv6
set life_juv_hr6 lput mean rep_success_hr_6 life_juv_hr6
    ]
  if length rep_success_7 > 0
  [
set life_juv7 lput mean rep_success_7 life_juv7
set life_juv_hr7 lput mean rep_success_hr_7 life_juv_hr7
    ]
  if length rep_success_8 > 0
  [
set life_juv8 lput mean rep_success_8 life_juv8
set life_juv_hr8 lput mean rep_success_hr_8 life_juv_hr8
    ]
  if length rep_success_9 > 0
  [
set life_juv9 lput mean rep_success_9 life_juv9
set life_juv_hr9 lput mean rep_success_hr_9 life_juv_hr9
  ]


  if stages_output
  [

  if any? turtles with [species = 0 and age < age_first_repro]
  [
set field_young0 lput mean [fmr] of turtles with [species = 0 and age < age_first_repro] field_young0
set bas_young0 lput mean [fmr_basal / fmr] of turtles with [species = 0 and age < age_first_repro] bas_young0
set loco_young0 lput mean [fmr_loco / fmr] of turtles with [species = 0 and age < age_first_repro] loco_young0
set repro_young0 lput mean [fmr_repro / fmr] of turtles with [species = 0 and age < age_first_repro] repro_young0
set grow_young0 lput mean [fmr_growth / fmr] of turtles with [species = 0 and age < age_first_repro] grow_young0
set dig_young0 lput mean [fmr_digest / fmr] of turtles with [species = 0 and age < age_first_repro] dig_young0
  ]
  if any? turtles with [species = 1 and age < age_first_repro]
  [
set field_young1 lput mean [fmr] of turtles with [species = 1 and age < age_first_repro] field_young1
set bas_young1 lput mean [fmr_basal / fmr] of turtles with [species = 1 and age < age_first_repro] bas_young1
set loco_young1 lput mean [fmr_loco / fmr] of turtles with [species = 1 and age < age_first_repro] loco_young1
set repro_young1 lput mean [fmr_repro / fmr] of turtles with [species = 1 and age < age_first_repro] repro_young1
set grow_young1 lput mean [fmr_growth / fmr] of turtles with [species = 1 and age < age_first_repro] grow_young1
set dig_young1 lput mean [fmr_digest / fmr] of turtles with [species = 1 and age < age_first_repro] dig_young1
    ]
  if any? turtles with [species = 2 and age < age_first_repro]
  [
set field_young2 lput mean [fmr] of turtles with [species = 2 and age < age_first_repro] field_young2
set bas_young2 lput mean [fmr_basal / fmr] of turtles with [species = 2 and age < age_first_repro] bas_young2
set loco_young2 lput mean [fmr_loco / fmr] of turtles with [species = 2 and age < age_first_repro] loco_young2
set repro_young2 lput mean [fmr_repro / fmr] of turtles with [species = 2 and age < age_first_repro] repro_young2
set grow_young2 lput mean [fmr_growth / fmr] of turtles with [species = 2 and age < age_first_repro] grow_young2
set dig_young2 lput mean [fmr_digest / fmr] of turtles with [species = 2 and age < age_first_repro] dig_young2
    ]
  if any? turtles with [species = 3 and age < age_first_repro]
  [
set field_young3 lput mean [fmr] of turtles with [species = 3 and age < age_first_repro] field_young3
set bas_young3 lput mean [fmr_basal / fmr] of turtles with [species = 3 and age < age_first_repro] bas_young3
set loco_young3 lput mean [fmr_loco / fmr] of turtles with [species = 3 and age < age_first_repro] loco_young3
set repro_young3 lput mean [fmr_repro / fmr] of turtles with [species = 3 and age < age_first_repro] repro_young3
set grow_young3 lput mean [fmr_growth / fmr] of turtles with [species = 3 and age < age_first_repro] grow_young3
set dig_young3 lput mean [fmr_digest / fmr] of turtles with [species = 3 and age < age_first_repro] dig_young3
    ]
  if any? turtles with [species = 4 and age < age_first_repro]
  [
set field_young4 lput mean [fmr] of turtles with [species = 4 and age < age_first_repro] field_young4
set bas_young4 lput mean [fmr_basal / fmr] of turtles with [species = 4 and age < age_first_repro] bas_young4
set loco_young4 lput mean [fmr_loco / fmr] of turtles with [species = 4 and age < age_first_repro] loco_young4
set repro_young4 lput mean [fmr_repro / fmr] of turtles with [species = 4 and age < age_first_repro] repro_young4
set grow_young4 lput mean [fmr_growth / fmr] of turtles with [species = 4 and age < age_first_repro] grow_young4
set dig_young4 lput mean [fmr_digest / fmr] of turtles with [species = 4 and age < age_first_repro] dig_young4
    ]
  if any? turtles with [species = 5 and age < age_first_repro]
  [
set field_young5 lput mean [fmr] of turtles with [species = 5 and age < age_first_repro] field_young5
set bas_young5 lput mean [fmr_basal / fmr] of turtles with [species = 5 and age < age_first_repro] bas_young5
set loco_young5 lput mean [fmr_loco / fmr] of turtles with [species = 5 and age < age_first_repro] loco_young5
set repro_young5 lput mean [fmr_repro / fmr] of turtles with [species = 5 and age < age_first_repro] repro_young5
set grow_young5 lput mean [fmr_growth / fmr] of turtles with [species = 5 and age < age_first_repro] grow_young5
set dig_young5 lput mean [fmr_digest / fmr] of turtles with [species = 5 and age < age_first_repro] dig_young5
    ]
  if any? turtles with [species = 6 and age < age_first_repro]
  [
set field_young6 lput mean [fmr] of turtles with [species = 6 and age < age_first_repro] field_young6
set bas_young6 lput mean [fmr_basal / fmr] of turtles with [species = 6 and age < age_first_repro] bas_young6
set loco_young6 lput mean [fmr_loco / fmr] of turtles with [species = 6 and age < age_first_repro] loco_young6
set repro_young6 lput mean [fmr_repro / fmr] of turtles with [species = 6 and age < age_first_repro] repro_young6
set grow_young6 lput mean [fmr_growth / fmr] of turtles with [species = 6 and age < age_first_repro] grow_young6
set dig_young6 lput mean [fmr_digest / fmr] of turtles with [species = 6 and age < age_first_repro] dig_young6
    ]
  if any? turtles with [species = 7 and age < age_first_repro]
  [
set field_young7 lput mean [fmr] of turtles with [species = 7 and age < age_first_repro] field_young7
set bas_young7 lput mean [fmr_basal / fmr] of turtles with [species = 7 and age < age_first_repro] bas_young7
set loco_young7 lput mean [fmr_loco / fmr] of turtles with [species = 7 and age < age_first_repro] loco_young7
set repro_young7 lput mean [fmr_repro / fmr] of turtles with [species = 7 and age < age_first_repro] repro_young7
set grow_young7 lput mean [fmr_growth / fmr] of turtles with [species = 7 and age < age_first_repro] grow_young7
set dig_young7 lput mean [fmr_digest / fmr] of turtles with [species = 7 and age < age_first_repro] dig_young7
    ]
  if any? turtles with [species = 8 and age < age_first_repro]
  [
set field_young8 lput mean [fmr] of turtles with [species = 8 and age < age_first_repro] field_young8
set bas_young8 lput mean [fmr_basal / fmr] of turtles with [species = 8 and age < age_first_repro] bas_young8
set loco_young8 lput mean [fmr_loco / fmr] of turtles with [species = 8 and age < age_first_repro] loco_young8
set repro_young8 lput mean [fmr_repro / fmr] of turtles with [species = 8 and age < age_first_repro] repro_young8
set grow_young8 lput mean [fmr_growth / fmr] of turtles with [species = 8 and age < age_first_repro] grow_young8
set dig_young8 lput mean [fmr_digest / fmr] of turtles with [species = 8 and age < age_first_repro] dig_young8
    ]
  if any? turtles with [species = 9 and age < age_first_repro]
  [
set field_young9 lput mean [fmr] of turtles with [species = 9 and age < age_first_repro] field_young9
set bas_young9 lput mean [fmr_basal / fmr] of turtles with [species = 9 and age < age_first_repro] bas_young9
set loco_young9 lput mean [fmr_loco / fmr] of turtles with [species = 9 and age < age_first_repro] loco_young9
set repro_young9 lput mean [fmr_repro / fmr] of turtles with [species = 9 and age < age_first_repro] repro_young9
set grow_young9 lput mean [fmr_growth / fmr] of turtles with [species = 9 and age < age_first_repro] grow_young9
set dig_young9 lput mean [fmr_digest / fmr] of turtles with [species = 9 and age < age_first_repro] dig_young9
  ]

  if any? turtles with [species = 0 and age > age_first_repro and sex = "male"]
  [
set field_mal0 lput mean [fmr] of turtles with [species = 0 and age > age_first_repro and sex = "male"] field_mal0
set bas_mal0 lput mean [fmr_basal / fmr] of turtles with [species = 0 and age > age_first_repro and sex = "male"] bas_mal0
set loco_mal0 lput mean [fmr_loco / fmr] of turtles with [species = 0 and age > age_first_repro and sex = "male"] loco_mal0
set repro_mal0 lput mean [fmr_repro / fmr] of turtles with [species = 0 and age > age_first_repro and sex = "male"] repro_mal0
set grow_mal0 lput mean [fmr_growth / fmr] of turtles with [species = 0 and age > age_first_repro and sex = "male"] grow_mal0
set dig_mal0 lput mean [fmr_digest / fmr] of turtles with [species = 0 and age > age_first_repro and sex = "male"] dig_mal0
  ]
  if any? turtles with [species = 1 and age > age_first_repro and sex = "male"]
  [
set field_mal1 lput mean [fmr] of turtles with [species = 1 and age > age_first_repro and sex = "male"] field_mal1
set bas_mal1 lput mean [fmr_basal / fmr] of turtles with [species = 1 and age > age_first_repro and sex = "male"] bas_mal1
set loco_mal1 lput mean [fmr_loco / fmr] of turtles with [species = 1 and age > age_first_repro and sex = "male"] loco_mal1
set repro_mal1 lput mean [fmr_repro / fmr] of turtles with [species = 1 and age > age_first_repro and sex = "male"] repro_mal1
set grow_mal1 lput mean [fmr_growth / fmr] of turtles with [species = 1 and age > age_first_repro and sex = "male"] grow_mal1
set dig_mal1 lput mean [fmr_digest / fmr] of turtles with [species = 1 and age > age_first_repro and sex = "male"] dig_mal1
     ]
  if any? turtles with [species = 2 and age > age_first_repro and sex = "male"]
  [
set field_mal2 lput mean [fmr] of turtles with [species = 2 and age > age_first_repro and sex = "male"] field_mal2
set bas_mal2 lput mean [fmr_basal / fmr] of turtles with [species = 2 and age > age_first_repro and sex = "male"] bas_mal2
set loco_mal2 lput mean [fmr_loco / fmr] of turtles with [species = 2 and age > age_first_repro and sex = "male"] loco_mal2
set repro_mal2 lput mean [fmr_repro / fmr] of turtles with [species = 2 and age > age_first_repro and sex = "male"] repro_mal2
set grow_mal2 lput mean [fmr_growth / fmr] of turtles with [species = 2 and age > age_first_repro and sex = "male"] grow_mal2
set dig_mal2 lput mean [fmr_digest / fmr] of turtles with [species = 2 and age > age_first_repro and sex = "male"] dig_mal2
     ]
  if any? turtles with [species = 3 and age > age_first_repro and sex = "male"]
  [
set field_mal3 lput mean [fmr] of turtles with [species = 3 and age > age_first_repro and sex = "male"] field_mal3
set bas_mal3 lput mean [fmr_basal / fmr] of turtles with [species = 3 and age > age_first_repro and sex = "male"] bas_mal3
set loco_mal3 lput mean [fmr_loco / fmr] of turtles with [species = 3 and age > age_first_repro and sex = "male"] loco_mal3
set repro_mal3 lput mean [fmr_repro / fmr] of turtles with [species = 3 and age > age_first_repro and sex = "male"] repro_mal3
set grow_mal3 lput mean [fmr_growth / fmr] of turtles with [species = 3 and age > age_first_repro and sex = "male"] grow_mal3
set dig_mal3 lput mean [fmr_digest / fmr] of turtles with [species = 3 and age > age_first_repro and sex = "male"] dig_mal3
     ]
  if any? turtles with [species = 4 and age > age_first_repro and sex = "male"]
  [
set field_mal4 lput mean [fmr] of turtles with [species = 4 and age > age_first_repro and sex = "male"] field_mal4
set bas_mal4 lput mean [fmr_basal / fmr] of turtles with [species = 4 and age > age_first_repro and sex = "male"] bas_mal4
set loco_mal4 lput mean [fmr_loco / fmr] of turtles with [species = 4 and age > age_first_repro and sex = "male"] loco_mal4
set repro_mal4 lput mean [fmr_repro / fmr] of turtles with [species = 4 and age > age_first_repro and sex = "male"] repro_mal4
set grow_mal4 lput mean [fmr_growth / fmr] of turtles with [species = 4 and age > age_first_repro and sex = "male"] grow_mal4
set dig_mal4 lput mean [fmr_digest / fmr] of turtles with [species = 4 and age > age_first_repro and sex = "male"] dig_mal4
     ]
  if any? turtles with [species = 5 and age > age_first_repro and sex = "male"]
  [
set field_mal5 lput mean [fmr] of turtles with [species = 5 and age > age_first_repro and sex = "male"] field_mal5
set bas_mal5 lput mean [fmr_basal / fmr] of turtles with [species = 5 and age > age_first_repro and sex = "male"] bas_mal5
set loco_mal5 lput mean [fmr_loco / fmr] of turtles with [species = 5 and age > age_first_repro and sex = "male"] loco_mal5
set repro_mal5 lput mean [fmr_repro / fmr] of turtles with [species = 5 and age > age_first_repro and sex = "male"] repro_mal5
set grow_mal5 lput mean [fmr_growth / fmr] of turtles with [species = 5 and age > age_first_repro and sex = "male"] grow_mal5
set dig_mal5 lput mean [fmr_digest / fmr] of turtles with [species = 5 and age > age_first_repro and sex = "male"] dig_mal5
     ]
  if any? turtles with [species = 6 and age > age_first_repro and sex = "male"]
  [
set field_mal6 lput mean [fmr] of turtles with [species = 6 and age > age_first_repro and sex = "male"] field_mal6
set bas_mal6 lput mean [fmr_basal / fmr] of turtles with [species = 6 and age > age_first_repro and sex = "male"] bas_mal6
set loco_mal6 lput mean [fmr_loco / fmr] of turtles with [species = 6 and age > age_first_repro and sex = "male"] loco_mal6
set repro_mal6 lput mean [fmr_repro / fmr] of turtles with [species = 6 and age > age_first_repro and sex = "male"] repro_mal6
set grow_mal6 lput mean [fmr_growth / fmr] of turtles with [species = 6 and age > age_first_repro and sex = "male"] grow_mal6
set dig_mal6 lput mean [fmr_digest / fmr] of turtles with [species = 6 and age > age_first_repro and sex = "male"] dig_mal6
     ]
  if any? turtles with [species = 7 and age > age_first_repro and sex = "male"]
  [
set field_mal7 lput mean [fmr] of turtles with [species = 7 and age > age_first_repro and sex = "male"] field_mal7
set bas_mal7 lput mean [fmr_basal / fmr] of turtles with [species = 7 and age > age_first_repro and sex = "male"] bas_mal7
set loco_mal7 lput mean [fmr_loco / fmr] of turtles with [species = 7 and age > age_first_repro and sex = "male"] loco_mal7
set repro_mal7 lput mean [fmr_repro / fmr] of turtles with [species = 7 and age > age_first_repro and sex = "male"] repro_mal7
set grow_mal7 lput mean [fmr_growth / fmr] of turtles with [species = 7 and age > age_first_repro and sex = "male"] grow_mal7
set dig_mal7 lput mean [fmr_digest / fmr] of turtles with [species = 7 and age > age_first_repro and sex = "male"] dig_mal7
     ]
  if any? turtles with [species = 8 and age > age_first_repro and sex = "male"]
  [
set field_mal8 lput mean [fmr] of turtles with [species = 8 and age > age_first_repro and sex = "male"] field_mal8
set bas_mal8 lput mean [fmr_basal / fmr] of turtles with [species = 8 and age > age_first_repro and sex = "male"] bas_mal8
set loco_mal8 lput mean [fmr_loco / fmr] of turtles with [species = 8 and age > age_first_repro and sex = "male"] loco_mal8
set repro_mal8 lput mean [fmr_repro / fmr] of turtles with [species = 8 and age > age_first_repro and sex = "male"] repro_mal8
set grow_mal8 lput mean [fmr_growth / fmr] of turtles with [species = 8 and age > age_first_repro and sex = "male"] grow_mal8
set dig_mal8 lput mean [fmr_digest / fmr] of turtles with [species = 8 and age > age_first_repro and sex = "male"] dig_mal8
     ]
  if any? turtles with [species = 9 and age > age_first_repro and sex = "male"]
  [
set field_mal9 lput mean [fmr] of turtles with [species = 9 and age > age_first_repro and sex = "male"] field_mal9
set bas_mal9 lput mean [fmr_basal / fmr] of turtles with [species = 9 and age > age_first_repro and sex = "male"] bas_mal9
set loco_mal9 lput mean [fmr_loco / fmr] of turtles with [species = 9 and age > age_first_repro and sex = "male"] loco_mal9
set repro_mal9 lput mean [fmr_repro / fmr] of turtles with [species = 9 and age > age_first_repro and sex = "male"] repro_mal9
set grow_mal9 lput mean [fmr_growth / fmr] of turtles with [species = 9 and age > age_first_repro and sex = "male"] grow_mal9
set dig_mal9 lput mean [fmr_digest / fmr] of turtles with [species = 9 and age > age_first_repro and sex = "male"] dig_mal9
  ]

  if any? turtles with [species = 0 and age > age_first_repro and preg = 1]
  [
set field_preg0 lput mean [fmr] of turtles with [species = 0 and age > age_first_repro and preg = 1] field_preg0
set bas_preg0 lput mean [fmr_basal / fmr] of turtles with [species = 0 and age > age_first_repro and preg = 1] bas_preg0
set loco_preg0 lput mean [fmr_loco / fmr] of turtles with [species = 0 and age > age_first_repro and preg = 1] loco_preg0
set repro_preg0 lput mean [fmr_repro / fmr] of turtles with [species = 0 and age > age_first_repro and preg = 1] repro_preg0
set grow_preg0 lput mean [fmr_growth / fmr] of turtles with [species = 0 and age > age_first_repro and preg = 1] grow_preg0
set dig_preg0 lput mean [fmr_digest / fmr] of turtles with [species = 0 and age > age_first_repro and preg = 1] dig_preg0
  ]
  if any? turtles with [species = 1 and age > age_first_repro and preg = 1]
  [
set field_preg1 lput mean [fmr] of turtles with [species = 1 and age > age_first_repro and preg = 1] field_preg1
set bas_preg1 lput mean [fmr_basal / fmr] of turtles with [species = 1 and age > age_first_repro and preg = 1] bas_preg1
set loco_preg1 lput mean [fmr_loco / fmr] of turtles with [species = 1 and age > age_first_repro and preg = 1] loco_preg1
set repro_preg1 lput mean [fmr_repro / fmr] of turtles with [species = 1 and age > age_first_repro and preg = 1] repro_preg1
set grow_preg1 lput mean [fmr_growth / fmr] of turtles with [species = 1 and age > age_first_repro and preg = 1] grow_preg1
set dig_preg1 lput mean [fmr_digest / fmr] of turtles with [species = 1 and age > age_first_repro and preg = 1] dig_preg1
    ]
  if any? turtles with [species = 2 and age > age_first_repro and preg = 1]
  [
set field_preg2 lput mean [fmr] of turtles with [species = 2 and age > age_first_repro and preg = 1] field_preg2
set bas_preg2 lput mean [fmr_basal / fmr] of turtles with [species = 2 and age > age_first_repro and preg = 1] bas_preg2
set loco_preg2 lput mean [fmr_loco / fmr] of turtles with [species = 2 and age > age_first_repro and preg = 1] loco_preg2
set repro_preg2 lput mean [fmr_repro / fmr] of turtles with [species = 2 and age > age_first_repro and preg = 1] repro_preg2
set grow_preg2 lput mean [fmr_growth / fmr] of turtles with [species = 2 and age > age_first_repro and preg = 1] grow_preg2
set dig_preg2 lput mean [fmr_digest / fmr] of turtles with [species = 2 and age > age_first_repro and preg = 1] dig_preg2
    ]
  if any? turtles with [species = 3 and age > age_first_repro and preg = 1]
  [
set field_preg3 lput mean [fmr] of turtles with [species = 3 and age > age_first_repro and preg = 1] field_preg3
set bas_preg3 lput mean [fmr_basal / fmr] of turtles with [species = 3 and age > age_first_repro and preg = 1] bas_preg3
set loco_preg3 lput mean [fmr_loco / fmr] of turtles with [species = 3 and age > age_first_repro and preg = 1] loco_preg3
set repro_preg3 lput mean [fmr_repro / fmr] of turtles with [species = 3 and age > age_first_repro and preg = 1] repro_preg3
set grow_preg3 lput mean [fmr_growth / fmr] of turtles with [species = 3 and age > age_first_repro and preg = 1] grow_preg3
set dig_preg3 lput mean [fmr_digest / fmr] of turtles with [species = 3 and age > age_first_repro and preg = 1] dig_preg3
    ]
  if any? turtles with [species = 4 and age > age_first_repro and preg = 1]
  [
set field_preg4 lput mean [fmr] of turtles with [species = 4 and age > age_first_repro and preg = 1] field_preg4
set bas_preg4 lput mean [fmr_basal / fmr] of turtles with [species = 4 and age > age_first_repro and preg = 1] bas_preg4
set loco_preg4 lput mean [fmr_loco / fmr] of turtles with [species = 4 and age > age_first_repro and preg = 1] loco_preg4
set repro_preg4 lput mean [fmr_repro / fmr] of turtles with [species = 4 and age > age_first_repro and preg = 1] repro_preg4
set grow_preg4 lput mean [fmr_growth / fmr] of turtles with [species = 4 and age > age_first_repro and preg = 1] grow_preg4
set dig_preg4 lput mean [fmr_digest / fmr] of turtles with [species = 4 and age > age_first_repro and preg = 1] dig_preg4
    ]
  if any? turtles with [species = 5 and age > age_first_repro and preg = 1]
  [
set field_preg5 lput mean [fmr] of turtles with [species = 5 and age > age_first_repro and preg = 1] field_preg5
set bas_preg5 lput mean [fmr_basal / fmr] of turtles with [species = 5 and age > age_first_repro and preg = 1] bas_preg5
set loco_preg5 lput mean [fmr_loco / fmr] of turtles with [species = 5 and age > age_first_repro and preg = 1] loco_preg5
set repro_preg5 lput mean [fmr_repro / fmr] of turtles with [species = 5 and age > age_first_repro and preg = 1] repro_preg5
set grow_preg5 lput mean [fmr_growth / fmr] of turtles with [species = 5 and age > age_first_repro and preg = 1] grow_preg5
set dig_preg5 lput mean [fmr_digest / fmr] of turtles with [species = 5 and age > age_first_repro and preg = 1] dig_preg5
    ]
  if any? turtles with [species = 6 and age > age_first_repro and preg = 1]
  [
set field_preg6 lput mean [fmr] of turtles with [species = 6 and age > age_first_repro and preg = 1] field_preg6
set bas_preg6 lput mean [fmr_basal / fmr] of turtles with [species = 6 and age > age_first_repro and preg = 1] bas_preg6
set loco_preg6 lput mean [fmr_loco / fmr] of turtles with [species = 6 and age > age_first_repro and preg = 1] loco_preg6
set repro_preg6 lput mean [fmr_repro / fmr] of turtles with [species = 6 and age > age_first_repro and preg = 1] repro_preg6
set grow_preg6 lput mean [fmr_growth / fmr] of turtles with [species = 6 and age > age_first_repro and preg = 1] grow_preg6
set dig_preg6 lput mean [fmr_digest / fmr] of turtles with [species = 6 and age > age_first_repro and preg = 1] dig_preg6
    ]
  if any? turtles with [species = 7 and age > age_first_repro and preg = 1]
  [
set field_preg7 lput mean [fmr] of turtles with [species = 7 and age > age_first_repro and preg = 1] field_preg7
set bas_preg7 lput mean [fmr_basal / fmr] of turtles with [species = 7 and age > age_first_repro and preg = 1] bas_preg7
set loco_preg7 lput mean [fmr_loco / fmr] of turtles with [species = 7 and age > age_first_repro and preg = 1] loco_preg7
set repro_preg7 lput mean [fmr_repro / fmr] of turtles with [species = 7 and age > age_first_repro and preg = 1] repro_preg7
set grow_preg7 lput mean [fmr_growth / fmr] of turtles with [species = 7 and age > age_first_repro and preg = 1] grow_preg7
set dig_preg7 lput mean [fmr_digest / fmr] of turtles with [species = 7 and age > age_first_repro and preg = 1] dig_preg7
    ]
  if any? turtles with [species = 8 and age > age_first_repro and preg = 1]
  [
set field_preg8 lput mean [fmr] of turtles with [species = 8 and age > age_first_repro and preg = 1] field_preg8
set bas_preg8 lput mean [fmr_basal / fmr] of turtles with [species = 8 and age > age_first_repro and preg = 1] bas_preg8
set loco_preg8 lput mean [fmr_loco / fmr] of turtles with [species = 8 and age > age_first_repro and preg = 1] loco_preg8
set repro_preg8 lput mean [fmr_repro / fmr] of turtles with [species = 8 and age > age_first_repro and preg = 1] repro_preg8
set grow_preg8 lput mean [fmr_growth / fmr] of turtles with [species = 8 and age > age_first_repro and preg = 1] grow_preg8
set dig_preg8 lput mean [fmr_digest / fmr] of turtles with [species = 8 and age > age_first_repro and preg = 1] dig_preg8
    ]
  if any? turtles with [species = 9 and age > age_first_repro and preg = 1]
  [
set field_preg9 lput mean [fmr] of turtles with [species = 9 and age > age_first_repro and preg = 1] field_preg9
set bas_preg9 lput mean [fmr_basal / fmr] of turtles with [species = 9 and age > age_first_repro and preg = 1] bas_preg9
set loco_preg9 lput mean [fmr_loco / fmr] of turtles with [species = 9 and age > age_first_repro and preg = 1] loco_preg9
set repro_preg9 lput mean [fmr_repro / fmr] of turtles with [species = 9 and age > age_first_repro and preg = 1] repro_preg9
set grow_preg9 lput mean [fmr_growth / fmr] of turtles with [species = 9 and age > age_first_repro and preg = 1] grow_preg9
set dig_preg9 lput mean [fmr_digest / fmr] of turtles with [species = 9 and age > age_first_repro and preg = 1] dig_preg9
  ]

  ]



  if long_output and ticks > 365
  [
  foreach range (1766.53 * (maxmass ^ 0.21) - 1)
 [
  ?1 ->
  set field_spec0 replace-item ?1 field_spec0 (sentence (fieldrate 0 ?1) (item ?1 field_spec0))
  set field_spec1 replace-item ?1 field_spec1 (sentence (fieldrate 1 ?1) (item ?1 field_spec1))
  set field_spec2 replace-item ?1 field_spec2 (sentence (fieldrate 2 ?1) (item ?1 field_spec2))
  set field_spec3 replace-item ?1 field_spec3 (sentence (fieldrate 3 ?1) (item ?1 field_spec3))
  set field_spec4 replace-item ?1 field_spec4 (sentence (fieldrate 4 ?1) (item ?1 field_spec4))
  set field_spec5 replace-item ?1 field_spec5 (sentence (fieldrate 5 ?1) (item ?1 field_spec5))
  set field_spec6 replace-item ?1 field_spec6 (sentence (fieldrate 6 ?1) (item ?1 field_spec6))
  set field_spec7 replace-item ?1 field_spec7 (sentence (fieldrate 7 ?1) (item ?1 field_spec7))
  set field_spec8 replace-item ?1 field_spec8 (sentence (fieldrate 8 ?1) (item ?1 field_spec8))
  set field_spec9 replace-item ?1 field_spec9 (sentence (fieldrate 9 ?1) (item ?1 field_spec9))

  set bas_spec0 replace-item ?1 bas_spec0 (sentence (basrate 0 ?1) (item ?1 bas_spec0))
  set bas_spec1 replace-item ?1 bas_spec1 (sentence (basrate 1 ?1) (item ?1 bas_spec1))
  set bas_spec2 replace-item ?1 bas_spec2 (sentence (basrate 2 ?1) (item ?1 bas_spec2))
  set bas_spec3 replace-item ?1 bas_spec3 (sentence (basrate 3 ?1) (item ?1 bas_spec3))
  set bas_spec4 replace-item ?1 bas_spec4 (sentence (basrate 4 ?1) (item ?1 bas_spec4))
  set bas_spec5 replace-item ?1 bas_spec5 (sentence (basrate 5 ?1) (item ?1 bas_spec5))
  set bas_spec6 replace-item ?1 bas_spec6 (sentence (basrate 6 ?1) (item ?1 bas_spec6))
  set bas_spec7 replace-item ?1 bas_spec7 (sentence (basrate 7 ?1) (item ?1 bas_spec7))
  set bas_spec8 replace-item ?1 bas_spec8 (sentence (basrate 8 ?1) (item ?1 bas_spec8))
  set bas_spec9 replace-item ?1 bas_spec9 (sentence (basrate 9 ?1) (item ?1 bas_spec9))

    set grow_spec0 replace-item ?1 grow_spec0 (sentence (growrate 0 ?1) (item ?1 grow_spec0))
  set grow_spec1 replace-item ?1 grow_spec1 (sentence (growrate 1 ?1) (item ?1 grow_spec1))
  set grow_spec2 replace-item ?1 grow_spec2 (sentence (growrate 2 ?1) (item ?1 grow_spec2))
  set grow_spec3 replace-item ?1 grow_spec3 (sentence (growrate 3 ?1) (item ?1 grow_spec3))
  set grow_spec4 replace-item ?1 grow_spec4 (sentence (growrate 4 ?1) (item ?1 grow_spec4))
  set grow_spec5 replace-item ?1 grow_spec5 (sentence (growrate 5 ?1) (item ?1 grow_spec5))
  set grow_spec6 replace-item ?1 grow_spec6 (sentence (growrate 6 ?1) (item ?1 grow_spec6))
  set grow_spec7 replace-item ?1 grow_spec7 (sentence (growrate 7 ?1) (item ?1 grow_spec7))
  set grow_spec8 replace-item ?1 grow_spec8 (sentence (growrate 8 ?1) (item ?1 grow_spec8))
  set grow_spec9 replace-item ?1 grow_spec9 (sentence (growrate 9 ?1) (item ?1 grow_spec9))

    set repro_spec0 replace-item ?1 repro_spec0 (sentence (reprorate 0 ?1) (item ?1 repro_spec0))
  set repro_spec1 replace-item ?1 repro_spec1 (sentence (reprorate 1 ?1) (item ?1 repro_spec1))
  set repro_spec2 replace-item ?1 repro_spec2 (sentence (reprorate 2 ?1) (item ?1 repro_spec2))
  set repro_spec3 replace-item ?1 repro_spec3 (sentence (reprorate 3 ?1) (item ?1 repro_spec3))
  set repro_spec4 replace-item ?1 repro_spec4 (sentence (reprorate 4 ?1) (item ?1 repro_spec4))
  set repro_spec5 replace-item ?1 repro_spec5 (sentence (reprorate 5 ?1) (item ?1 repro_spec5))
  set repro_spec6 replace-item ?1 repro_spec6 (sentence (reprorate 6 ?1) (item ?1 repro_spec6))
  set repro_spec7 replace-item ?1 repro_spec7 (sentence (reprorate 7 ?1) (item ?1 repro_spec7))
  set repro_spec8 replace-item ?1 repro_spec8 (sentence (reprorate 8 ?1) (item ?1 repro_spec8))
  set repro_spec9 replace-item ?1 repro_spec9 (sentence (reprorate 9 ?1) (item ?1 repro_spec9))

    set gest_spec0 replace-item ?1 gest_spec0 (sentence (gestrate 0 ?1) (item ?1 gest_spec0))
  set gest_spec1 replace-item ?1 gest_spec1 (sentence (gestrate 1 ?1) (item ?1 gest_spec1))
  set gest_spec2 replace-item ?1 gest_spec2 (sentence (gestrate 2 ?1) (item ?1 gest_spec2))
  set gest_spec3 replace-item ?1 gest_spec3 (sentence (gestrate 3 ?1) (item ?1 gest_spec3))
  set gest_spec4 replace-item ?1 gest_spec4 (sentence (gestrate 4 ?1) (item ?1 gest_spec4))
  set gest_spec5 replace-item ?1 gest_spec5 (sentence (gestrate 5 ?1) (item ?1 gest_spec5))
  set gest_spec6 replace-item ?1 gest_spec6 (sentence (gestrate 6 ?1) (item ?1 gest_spec6))
  set gest_spec7 replace-item ?1 gest_spec7 (sentence (gestrate 7 ?1) (item ?1 gest_spec7))
  set gest_spec8 replace-item ?1 gest_spec8 (sentence (gestrate 8 ?1) (item ?1 gest_spec8))
  set gest_spec9 replace-item ?1 gest_spec9 (sentence (gestrate 9 ?1) (item ?1 gest_spec9))

    set lact_spec0 replace-item ?1 lact_spec0 (sentence (lactrate 0 ?1) (item ?1 lact_spec0))
  set lact_spec1 replace-item ?1 lact_spec1 (sentence (lactrate 1 ?1) (item ?1 lact_spec1))
  set lact_spec2 replace-item ?1 lact_spec2 (sentence (lactrate 2 ?1) (item ?1 lact_spec2))
  set lact_spec3 replace-item ?1 lact_spec3 (sentence (lactrate 3 ?1) (item ?1 lact_spec3))
  set lact_spec4 replace-item ?1 lact_spec4 (sentence (lactrate 4 ?1) (item ?1 lact_spec4))
  set lact_spec5 replace-item ?1 lact_spec5 (sentence (lactrate 5 ?1) (item ?1 lact_spec5))
  set lact_spec6 replace-item ?1 lact_spec6 (sentence (lactrate 6 ?1) (item ?1 lact_spec6))
  set lact_spec7 replace-item ?1 lact_spec7 (sentence (lactrate 7 ?1) (item ?1 lact_spec7))
  set lact_spec8 replace-item ?1 lact_spec8 (sentence (lactrate 8 ?1) (item ?1 lact_spec8))
  set lact_spec9 replace-item ?1 lact_spec9 (sentence (lactrate 9 ?1) (item ?1 lact_spec9))

    set loco_spec0 replace-item ?1 loco_spec0 (sentence (locorate 0 ?1) (item ?1 loco_spec0))
  set loco_spec1 replace-item ?1 loco_spec1 (sentence (locorate 1 ?1) (item ?1 loco_spec1))
  set loco_spec2 replace-item ?1 loco_spec2 (sentence (locorate 2 ?1) (item ?1 loco_spec2))
  set loco_spec3 replace-item ?1 loco_spec3 (sentence (locorate 3 ?1) (item ?1 loco_spec3))
  set loco_spec4 replace-item ?1 loco_spec4 (sentence (locorate 4 ?1) (item ?1 loco_spec4))
  set loco_spec5 replace-item ?1 loco_spec5 (sentence (locorate 5 ?1) (item ?1 loco_spec5))
  set loco_spec6 replace-item ?1 loco_spec6 (sentence (locorate 6 ?1) (item ?1 loco_spec6))
  set loco_spec7 replace-item ?1 loco_spec7 (sentence (locorate 7 ?1) (item ?1 loco_spec7))
  set loco_spec8 replace-item ?1 loco_spec8 (sentence (locorate 8 ?1) (item ?1 loco_spec8))
  set loco_spec9 replace-item ?1 loco_spec9 (sentence (locorate 9 ?1) (item ?1 loco_spec9))

  set digest_spec0 replace-item ?1 digest_spec0 (sentence (digestrate 0 ?1) (item ?1 digest_spec0))
  set digest_spec1 replace-item ?1 digest_spec1 (sentence (digestrate 1 ?1) (item ?1 digest_spec1))
  set digest_spec2 replace-item ?1 digest_spec2 (sentence (digestrate 2 ?1) (item ?1 digest_spec2))
  set digest_spec3 replace-item ?1 digest_spec3 (sentence (digestrate 3 ?1) (item ?1 digest_spec3))
  set digest_spec4 replace-item ?1 digest_spec4 (sentence (digestrate 4 ?1) (item ?1 digest_spec4))
  set digest_spec5 replace-item ?1 digest_spec5 (sentence (digestrate 5 ?1) (item ?1 digest_spec5))
  set digest_spec6 replace-item ?1 digest_spec6 (sentence (digestrate 6 ?1) (item ?1 digest_spec6))
  set digest_spec7 replace-item ?1 digest_spec7 (sentence (digestrate 7 ?1) (item ?1 digest_spec7))
  set digest_spec8 replace-item ?1 digest_spec8 (sentence (digestrate 8 ?1) (item ?1 digest_spec8))
  set digest_spec9 replace-item ?1 digest_spec9 (sentence (digestrate 9 ?1) (item ?1 digest_spec9))

  set stor_spec0 replace-item ?1 stor_spec0 (sentence (storrate 0 ?1) (item ?1 stor_spec0))
  set stor_spec1 replace-item ?1 stor_spec1 (sentence (storrate 1 ?1) (item ?1 stor_spec1))
  set stor_spec2 replace-item ?1 stor_spec2 (sentence (storrate 2 ?1) (item ?1 stor_spec2))
  set stor_spec3 replace-item ?1 stor_spec3 (sentence (storrate 3 ?1) (item ?1 stor_spec3))
  set stor_spec4 replace-item ?1 stor_spec4 (sentence (storrate 4 ?1) (item ?1 stor_spec4))
  set stor_spec5 replace-item ?1 stor_spec5 (sentence (storrate 5 ?1) (item ?1 stor_spec5))
  set stor_spec6 replace-item ?1 stor_spec6 (sentence (storrate 6 ?1) (item ?1 stor_spec6))
  set stor_spec7 replace-item ?1 stor_spec7 (sentence (storrate 7 ?1) (item ?1 stor_spec7))
  set stor_spec8 replace-item ?1 stor_spec8 (sentence (storrate 8 ?1) (item ?1 stor_spec8))
  set stor_spec9 replace-item ?1 stor_spec9 (sentence (storrate 9 ?1) (item ?1 stor_spec9))

  set mass_spec0 replace-item ?1 mass_spec0 (sentence (massrate 0 ?1) (item ?1 mass_spec0))
  set mass_spec1 replace-item ?1 mass_spec1 (sentence (massrate 1 ?1) (item ?1 mass_spec1))
  set mass_spec2 replace-item ?1 mass_spec2 (sentence (massrate 2 ?1) (item ?1 mass_spec2))
  set mass_spec3 replace-item ?1 mass_spec3 (sentence (massrate 3 ?1) (item ?1 mass_spec3))
  set mass_spec4 replace-item ?1 mass_spec4 (sentence (massrate 4 ?1) (item ?1 mass_spec4))
  set mass_spec5 replace-item ?1 mass_spec5 (sentence (massrate 5 ?1) (item ?1 mass_spec5))
  set mass_spec6 replace-item ?1 mass_spec6 (sentence (massrate 6 ?1) (item ?1 mass_spec6))
  set mass_spec7 replace-item ?1 mass_spec7 (sentence (massrate 7 ?1) (item ?1 mass_spec7))
  set mass_spec8 replace-item ?1 mass_spec8 (sentence (massrate 8 ?1) (item ?1 mass_spec8))
  set mass_spec9 replace-item ?1 mass_spec9 (sentence (massrate 9 ?1) (item ?1 mass_spec9))

  set young_mass_spec0 replace-item ?1 young_mass_spec0 (sentence (youngmassrate 0 ?1) (item ?1 young_mass_spec0))
  set young_mass_spec1 replace-item ?1 young_mass_spec1 (sentence (youngmassrate 1 ?1) (item ?1 young_mass_spec1))
  set young_mass_spec2 replace-item ?1 young_mass_spec2 (sentence (youngmassrate 2 ?1) (item ?1 young_mass_spec2))
  set young_mass_spec3 replace-item ?1 young_mass_spec3 (sentence (youngmassrate 3 ?1) (item ?1 young_mass_spec3))
  set young_mass_spec4 replace-item ?1 young_mass_spec4 (sentence (youngmassrate 4 ?1) (item ?1 young_mass_spec4))
  set young_mass_spec5 replace-item ?1 young_mass_spec5 (sentence (youngmassrate 5 ?1) (item ?1 young_mass_spec5))
  set young_mass_spec6 replace-item ?1 young_mass_spec6 (sentence (youngmassrate 6 ?1) (item ?1 young_mass_spec6))
  set young_mass_spec7 replace-item ?1 young_mass_spec7 (sentence (youngmassrate 7 ?1) (item ?1 young_mass_spec7))
  set young_mass_spec8 replace-item ?1 young_mass_spec8 (sentence (youngmassrate 8 ?1) (item ?1 young_mass_spec8))
  set young_mass_spec9 replace-item ?1 young_mass_spec9 (sentence (youngmassrate 9 ?1) (item ?1 young_mass_spec9))
  ]
    foreach range 75
 [
  ?1 ->
      set digest_mother_spec0 replace-item ?1 digest_mother_spec0 (sentence (digestrate_mother 0 ?1) (item ?1 digest_mother_spec0))
  set digest_mother_spec1 replace-item ?1 digest_mother_spec1 (sentence (digestrate_mother 1 ?1) (item ?1 digest_mother_spec1))
  set digest_mother_spec2 replace-item ?1 digest_mother_spec2 (sentence (digestrate_mother 2 ?1) (item ?1 digest_mother_spec2))
  set digest_mother_spec3 replace-item ?1 digest_mother_spec3 (sentence (digestrate_mother 3 ?1) (item ?1 digest_mother_spec3))
  set digest_mother_spec4 replace-item ?1 digest_mother_spec4 (sentence (digestrate_mother 4 ?1) (item ?1 digest_mother_spec4))
  set digest_mother_spec5 replace-item ?1 digest_mother_spec5 (sentence (digestrate_mother 5 ?1) (item ?1 digest_mother_spec5))
  set digest_mother_spec6 replace-item ?1 digest_mother_spec6 (sentence (digestrate_mother 6 ?1) (item ?1 digest_mother_spec6))
  set digest_mother_spec7 replace-item ?1 digest_mother_spec7 (sentence (digestrate_mother 7 ?1) (item ?1 digest_mother_spec7))
  set digest_mother_spec8 replace-item ?1 digest_mother_spec8 (sentence (digestrate_mother 8 ?1) (item ?1 digest_mother_spec8))
  set digest_mother_spec9 replace-item ?1 digest_mother_spec9 (sentence (digestrate_mother 9 ?1) (item ?1 digest_mother_spec9))
    ]
  ]
    foreach range 5
    [
      ?1 ->
      set fmr_pup_gest_spec0 replace-item ?1 fmr_pup_gest_spec0 (sentence (perpup_gest 0 ?1) (item ?1 fmr_pup_gest_spec0))
      set fmr_pup_gest_spec1 replace-item ?1 fmr_pup_gest_spec1 (sentence (perpup_gest 1 ?1) (item ?1 fmr_pup_gest_spec1))
      set fmr_pup_gest_spec2 replace-item ?1 fmr_pup_gest_spec2 (sentence (perpup_gest 2 ?1) (item ?1 fmr_pup_gest_spec2))
      set fmr_pup_gest_spec3 replace-item ?1 fmr_pup_gest_spec3 (sentence (perpup_gest 3 ?1) (item ?1 fmr_pup_gest_spec3))
      set fmr_pup_gest_spec4 replace-item ?1 fmr_pup_gest_spec4 (sentence (perpup_gest 4 ?1) (item ?1 fmr_pup_gest_spec4))
      set fmr_pup_gest_spec5 replace-item ?1 fmr_pup_gest_spec5 (sentence (perpup_gest 5 ?1) (item ?1 fmr_pup_gest_spec5))
      set fmr_pup_gest_spec6 replace-item ?1 fmr_pup_gest_spec6 (sentence (perpup_gest 6 ?1) (item ?1 fmr_pup_gest_spec6))
      set fmr_pup_gest_spec7 replace-item ?1 fmr_pup_gest_spec7 (sentence (perpup_gest 7 ?1) (item ?1 fmr_pup_gest_spec7))
      set fmr_pup_gest_spec8 replace-item ?1 fmr_pup_gest_spec8 (sentence (perpup_gest 8 ?1) (item ?1 fmr_pup_gest_spec8))
      set fmr_pup_gest_spec9 replace-item ?1 fmr_pup_gest_spec9 (sentence (perpup_gest 9 ?1) (item ?1 fmr_pup_gest_spec9))

      set fmr_pup_lact_spec0 replace-item ?1 fmr_pup_lact_spec0 (sentence (perpup_lact 0 ?1) (item ?1 fmr_pup_lact_spec0))
      set fmr_pup_lact_spec1 replace-item ?1 fmr_pup_lact_spec1 (sentence (perpup_lact 1 ?1) (item ?1 fmr_pup_lact_spec1))
      set fmr_pup_lact_spec2 replace-item ?1 fmr_pup_lact_spec2 (sentence (perpup_lact 2 ?1) (item ?1 fmr_pup_lact_spec2))
      set fmr_pup_lact_spec3 replace-item ?1 fmr_pup_lact_spec3 (sentence (perpup_lact 3 ?1) (item ?1 fmr_pup_lact_spec3))
      set fmr_pup_lact_spec4 replace-item ?1 fmr_pup_lact_spec4 (sentence (perpup_lact 4 ?1) (item ?1 fmr_pup_lact_spec4))
      set fmr_pup_lact_spec5 replace-item ?1 fmr_pup_lact_spec5 (sentence (perpup_lact 5 ?1) (item ?1 fmr_pup_lact_spec5))
      set fmr_pup_lact_spec6 replace-item ?1 fmr_pup_lact_spec6 (sentence (perpup_lact 6 ?1) (item ?1 fmr_pup_lact_spec6))
      set fmr_pup_lact_spec7 replace-item ?1 fmr_pup_lact_spec7 (sentence (perpup_lact 7 ?1) (item ?1 fmr_pup_lact_spec7))
      set fmr_pup_lact_spec8 replace-item ?1 fmr_pup_lact_spec8 (sentence (perpup_lact 8 ?1) (item ?1 fmr_pup_lact_spec8))
      set fmr_pup_lact_spec9 replace-item ?1 fmr_pup_lact_spec9 (sentence (perpup_lact 9 ?1) (item ?1 fmr_pup_lact_spec9))
    ]

end

to-report fieldrate [s a]
  report [fmr] of turtles with [species = s and age = a and sex = "female"]
end

to-report basrate [s a]
  report [fmr_basal] of turtles with [species = s and age = a and sex = "female"]
end

to-report growrate [s a]
  report [fmr_growth] of turtles with [species = s and age = a and sex = "female"]
end

to-report reprorate [s a]
  report [fmr_repro] of turtles with [species = s and age = a and sex = "female"]
end

to-report gestrate [s a]
  let part [fmr_repro] of turtles with [species = s and age = a and sex = "female" and pregcount <= gest_period]
  report sentence part n-values (length [fmr_repro] of turtles with [species = s and age = a and sex = "female"] - length part) [0]
end

to-report lactrate [s a]
  let part [fmr_repro] of turtles with [species = s and age = a and sex = "female" and pregcount > gest_period and pregcount < (gest_period + lact_period)]
  report sentence part n-values (length [fmr_repro] of turtles with [species = s and age = a and sex = "female"] - length part) [0]
end

to-report locorate [s a]
  report [fmr_loco] of turtles with [species = s and age = a and sex = "female"]
end

to-report digestrate [s a]
  report [fmr_digest] of turtles with [species = s and age = a and sex = "female"]
end

to-report digestrate_mother [s a]
  report [fmr_digest] of turtles with [species = s and sex = "female" and preg = 1 and pregcount = a]
end

to-report storrate [s a]
  report [storage] of turtles with [species = s and age = a and sex = "female"]
end

to-report massrate [s a]
  report [real-mass] of turtles with [species = s and age = a and sex = "female"]
end

to-report youngmassrate [s a]
  report [real-young * mass-embryo] of turtles with [species = s and age = a and sex = "female"]
end

to-report perpup_gest [s p]
  report [fmr] of turtles with [species = s and preg = 1 and pregcount = round(gest_period) and real-young = p] ;; < instead of = ?
end

to-report perpup_lact [s p]
  report [fmr] of turtles with [species = s and preg = 1 and pregcount = round(gest_period + lact_period - 1) and real-young = p]
end

to-report nanmean [l]
  ifelse length l > 0
  [report mean l]
  [report 999]
end
@#$#@#$#@
GRAPHICS-WINDOW
212
13
688
490
-1
-1
4.634
1
6
1
1
1
0
1
1
1
0
100
0
100
1
1
1
ticks
30.0

BUTTON
29
17
109
51
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
29
55
111
89
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
1356
15
1609
219
Species richness
species ID
nr individuals
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -12895429 true "histogram [species] of turtles ;with [hrsize > 0]" "histogram [species] of turtles ;with [hrsize > 0]"

TEXTBOX
426
483
593
503
1km
11
0.0
1

PLOT
698
15
1347
220
Species numbers
time
number of ind.
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot count turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot count turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot count turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot count turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot count turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot count turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot count turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot count turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot count turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot count turtles with [species = 9]"

PLOT
1360
438
1609
643
Leftover resources
time
patchfeed
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"patchfeed" 1.0 0 -16777216 true "plot sum [patchfeed] of patches" "plot sum [patchfeed] of patches"

PLOT
699
227
1348
432
Mean energy storage
time
energy storage
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot mean [storage] of turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot mean [storage] of turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot mean [storage] of turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot mean [storage] of turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot mean [storage] of turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot mean [storage] of turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot mean [storage] of turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot mean [storage] of turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot mean [storage] of turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot mean [storage] of turtles with [species = 9]"

PLOT
1357
228
1609
431
Species Number
time
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot spec_num"

CHOOSER
32
194
205
239
clump
clump
0.9 0.99 0.999 0.9999
3

CHOOSER
32
363
170
408
specs-included
specs-included
"all" "two" 0 1 2 3 4 5 6 7 8 9
0

PLOT
700
438
1350
644
Mean field metabolic rate
time
FMR [g food]
0.0
12.0
0.0
12.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot mean [fmr] of turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot mean [fmr] of turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot mean [fmr] of turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot mean [fmr] of turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot mean [fmr] of turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot mean [fmr] of turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot mean [fmr] of turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot mean [fmr] of turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot mean [fmr] of turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot mean [fmr] of turtles with [species = 9]"

SLIDER
32
243
204
276
feed-amount
feed-amount
0
50
16.0
1
1
NIL
HORIZONTAL

SWITCH
30
520
191
553
long_output
long_output
1
1
-1000

CHOOSER
29
561
191
606
debug
debug
0 1
0

INPUTBOX
255
521
323
581
digestion
0.5
1
0
Number

INPUTBOX
330
521
395
581
mortbold
3.33E-5
1
0
Number

INPUTBOX
404
521
466
581
hr_try_juv
10.0
1
0
Number

INPUTBOX
404
588
466
648
maxmass
0.1
1
0
Number

INPUTBOX
256
585
323
645
prio_fact
0.2
1
0
Number

INPUTBOX
330
587
395
647
move-factor
3.0
1
0
Number

INPUTBOX
474
521
537
581
hr_try_init
100.0
1
0
Number

INPUTBOX
32
130
206
190
cover_percentage
0.05
1
0
Number

INPUTBOX
474
588
537
648
starting_indivs
1000.0
1
0
Number

TEXTBOX
29
113
179
131
Landscape Setup:
12
0.0
1

INPUTBOX
48
424
107
484
species1
8.0
1
0
Number

INPUTBOX
112
424
169
484
species2
9.0
1
0
Number

SWITCH
29
610
191
643
stages_output
stages_output
1
1
-1000

INPUTBOX
542
521
609
581
maturity
0.95
1
0
Number

INPUTBOX
32
280
204
340
feed-var
0.7
1
0
Number

INPUTBOX
542
588
609
648
storage_addition
3.0
1
0
Number

TEXTBOX
36
346
186
364
Species Setup:
12
0.0
1

TEXTBOX
51
409
201
427
if two species:
12
0.0
1

TEXTBOX
34
501
184
519
Output:
12
0.0
1

TEXTBOX
213
502
363
520
Parameters:
12
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment_frag_all" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3650"/>
    <metric>spec_num</metric>
    <metric>number0</metric>
    <metric>number1</metric>
    <metric>number2</metric>
    <metric>number3</metric>
    <metric>number4</metric>
    <metric>number5</metric>
    <metric>number6</metric>
    <metric>number7</metric>
    <metric>number8</metric>
    <metric>number9</metric>
    <metric>mean age0</metric>
    <metric>mean age1</metric>
    <metric>mean age2</metric>
    <metric>mean age3</metric>
    <metric>mean age4</metric>
    <metric>mean age5</metric>
    <metric>mean age6</metric>
    <metric>mean age7</metric>
    <metric>mean age8</metric>
    <metric>mean age9</metric>
    <metric>sum [real-mass] of turtles</metric>
    <metric>mean [real-mass] of turtles</metric>
    <metric>mean rep_success_0</metric>
    <metric>mean rep_success_1</metric>
    <metric>mean rep_success_2</metric>
    <metric>mean rep_success_3</metric>
    <metric>mean rep_success_4</metric>
    <metric>mean rep_success_5</metric>
    <metric>mean rep_success_6</metric>
    <metric>mean rep_success_7</metric>
    <metric>mean rep_success_8</metric>
    <metric>mean rep_success_9</metric>
    <metric>mean rep_success_hr_0</metric>
    <metric>mean rep_success_hr_1</metric>
    <metric>mean rep_success_hr_2</metric>
    <metric>mean rep_success_hr_3</metric>
    <metric>mean rep_success_hr_4</metric>
    <metric>mean rep_success_hr_5</metric>
    <metric>mean rep_success_hr_6</metric>
    <metric>mean rep_success_hr_7</metric>
    <metric>mean rep_success_hr_8</metric>
    <metric>mean rep_success_hr_9</metric>
    <metric>mean life_juv_hr0</metric>
    <metric>mean life_juv_hr1</metric>
    <metric>mean life_juv_hr2</metric>
    <metric>mean life_juv_hr3</metric>
    <metric>mean life_juv_hr4</metric>
    <metric>mean life_juv_hr5</metric>
    <metric>mean life_juv_hr6</metric>
    <metric>mean life_juv_hr7</metric>
    <metric>mean life_juv_hr8</metric>
    <metric>mean life_juv_hr9</metric>
    <metric>mean life_juv0</metric>
    <metric>mean life_juv1</metric>
    <metric>mean life_juv2</metric>
    <metric>mean life_juv3</metric>
    <metric>mean life_juv4</metric>
    <metric>mean life_juv5</metric>
    <metric>mean life_juv6</metric>
    <metric>mean life_juv7</metric>
    <metric>mean life_juv8</metric>
    <metric>mean life_juv9</metric>
    <metric>mean hr</metric>
    <metric>mean hr0</metric>
    <metric>mean hr1</metric>
    <metric>mean hr2</metric>
    <metric>mean hr3</metric>
    <metric>mean hr4</metric>
    <metric>mean hr5</metric>
    <metric>mean hr6</metric>
    <metric>mean hr7</metric>
    <metric>mean hr8</metric>
    <metric>mean hr9</metric>
    <metric>variance hr0</metric>
    <metric>variance hr1</metric>
    <metric>variance hr2</metric>
    <metric>variance hr3</metric>
    <metric>variance hr4</metric>
    <metric>variance hr5</metric>
    <metric>variance hr6</metric>
    <metric>variance hr7</metric>
    <metric>variance hr8</metric>
    <metric>variance hr9</metric>
    <metric>mean maxhr0</metric>
    <metric>mean maxhr1</metric>
    <metric>mean maxhr2</metric>
    <metric>mean maxhr3</metric>
    <metric>mean maxhr4</metric>
    <metric>mean maxhr5</metric>
    <metric>mean maxhr6</metric>
    <metric>mean maxhr7</metric>
    <metric>mean maxhr8</metric>
    <metric>mean maxhr9</metric>
    <metric>mean countmaxhr0</metric>
    <metric>mean countmaxhr1</metric>
    <metric>mean countmaxhr2</metric>
    <metric>mean countmaxhr3</metric>
    <metric>mean countmaxhr4</metric>
    <metric>mean countmaxhr5</metric>
    <metric>mean countmaxhr6</metric>
    <metric>mean countmaxhr7</metric>
    <metric>mean countmaxhr8</metric>
    <metric>mean countmaxhr9</metric>
    <metric>mean stor0</metric>
    <metric>mean stor1</metric>
    <metric>mean stor2</metric>
    <metric>mean stor3</metric>
    <metric>mean stor4</metric>
    <metric>mean stor5</metric>
    <metric>mean stor6</metric>
    <metric>mean stor7</metric>
    <metric>mean stor8</metric>
    <metric>mean stor9</metric>
    <metric>mean abs_stor0</metric>
    <metric>mean abs_stor1</metric>
    <metric>mean abs_stor2</metric>
    <metric>mean abs_stor3</metric>
    <metric>mean abs_stor4</metric>
    <metric>mean abs_stor5</metric>
    <metric>mean abs_stor6</metric>
    <metric>mean abs_stor7</metric>
    <metric>mean abs_stor8</metric>
    <metric>mean abs_stor9</metric>
    <metric>mean rm</metric>
    <metric>mean rm0</metric>
    <metric>mean rm1</metric>
    <metric>mean rm2</metric>
    <metric>mean rm3</metric>
    <metric>mean rm4</metric>
    <metric>mean rm5</metric>
    <metric>mean rm6</metric>
    <metric>mean rm7</metric>
    <metric>mean rm8</metric>
    <metric>mean rm9</metric>
    <metric>mean juv</metric>
    <metric>mean juv0</metric>
    <metric>mean juv1</metric>
    <metric>mean juv2</metric>
    <metric>mean juv3</metric>
    <metric>mean juv4</metric>
    <metric>mean juv5</metric>
    <metric>mean juv6</metric>
    <metric>mean juv7</metric>
    <metric>mean juv8</metric>
    <metric>mean juv9</metric>
    <metric>mean fmr0</metric>
    <metric>mean fmr1</metric>
    <metric>mean fmr2</metric>
    <metric>mean fmr3</metric>
    <metric>mean fmr4</metric>
    <metric>mean fmr5</metric>
    <metric>mean fmr6</metric>
    <metric>mean fmr7</metric>
    <metric>mean fmr8</metric>
    <metric>mean fmr9</metric>
    <metric>mean loco0</metric>
    <metric>mean loco1</metric>
    <metric>mean loco2</metric>
    <metric>mean loco3</metric>
    <metric>mean loco4</metric>
    <metric>mean loco5</metric>
    <metric>mean loco6</metric>
    <metric>mean loco7</metric>
    <metric>mean loco8</metric>
    <metric>mean loco9</metric>
    <metric>mean repro0</metric>
    <metric>mean repro1</metric>
    <metric>mean repro2</metric>
    <metric>mean repro3</metric>
    <metric>mean repro4</metric>
    <metric>mean repro5</metric>
    <metric>mean repro6</metric>
    <metric>mean repro7</metric>
    <metric>mean repro8</metric>
    <metric>mean repro9</metric>
    <metric>mean grow0</metric>
    <metric>mean grow1</metric>
    <metric>mean grow2</metric>
    <metric>mean grow3</metric>
    <metric>mean grow4</metric>
    <metric>mean grow5</metric>
    <metric>mean grow6</metric>
    <metric>mean grow7</metric>
    <metric>mean grow8</metric>
    <metric>mean grow9</metric>
    <metric>mean basal0</metric>
    <metric>mean basal1</metric>
    <metric>mean basal2</metric>
    <metric>mean basal3</metric>
    <metric>mean basal4</metric>
    <metric>mean basal5</metric>
    <metric>mean basal6</metric>
    <metric>mean basal7</metric>
    <metric>mean basal8</metric>
    <metric>mean basal9</metric>
    <metric>mean digest0</metric>
    <metric>mean digest1</metric>
    <metric>mean digest2</metric>
    <metric>mean digest3</metric>
    <metric>mean digest4</metric>
    <metric>mean digest5</metric>
    <metric>mean digest6</metric>
    <metric>mean digest7</metric>
    <metric>mean digest8</metric>
    <metric>mean digest9</metric>
    <metric>mean prod0</metric>
    <metric>mean prod1</metric>
    <metric>mean prod2</metric>
    <metric>mean prod3</metric>
    <metric>mean prod4</metric>
    <metric>mean prod5</metric>
    <metric>mean prod6</metric>
    <metric>mean prod7</metric>
    <metric>mean prod8</metric>
    <metric>mean prod9</metric>
    <metric>mean in0</metric>
    <metric>mean in1</metric>
    <metric>mean in2</metric>
    <metric>mean in3</metric>
    <metric>mean in4</metric>
    <metric>mean in5</metric>
    <metric>mean in6</metric>
    <metric>mean in7</metric>
    <metric>mean in8</metric>
    <metric>mean in9</metric>
    <metric>variance in0</metric>
    <metric>variance in1</metric>
    <metric>variance in2</metric>
    <metric>variance in3</metric>
    <metric>variance in4</metric>
    <metric>variance in5</metric>
    <metric>variance in6</metric>
    <metric>variance in7</metric>
    <metric>variance in8</metric>
    <metric>variance in9</metric>
    <metric>mean inpatch0</metric>
    <metric>mean inpatch1</metric>
    <metric>mean inpatch2</metric>
    <metric>mean inpatch3</metric>
    <metric>mean inpatch4</metric>
    <metric>mean inpatch5</metric>
    <metric>mean inpatch6</metric>
    <metric>mean inpatch7</metric>
    <metric>mean inpatch8</metric>
    <metric>mean inpatch9</metric>
    <metric>mean inhr0</metric>
    <metric>mean inhr1</metric>
    <metric>mean inhr2</metric>
    <metric>mean inhr3</metric>
    <metric>mean inhr4</metric>
    <metric>mean inhr5</metric>
    <metric>mean inhr6</metric>
    <metric>mean inhr7</metric>
    <metric>mean inhr8</metric>
    <metric>mean inhr9</metric>
    <metric>mean move0</metric>
    <metric>mean move1</metric>
    <metric>mean move2</metric>
    <metric>mean move3</metric>
    <metric>mean move4</metric>
    <metric>mean move5</metric>
    <metric>mean move6</metric>
    <metric>mean move7</metric>
    <metric>mean move8</metric>
    <metric>mean move9</metric>
    <metric>mean patches0</metric>
    <metric>mean patches1</metric>
    <metric>mean patches2</metric>
    <metric>mean patches3</metric>
    <metric>mean patches4</metric>
    <metric>mean patches5</metric>
    <metric>mean patches6</metric>
    <metric>mean patches7</metric>
    <metric>mean patches8</metric>
    <metric>mean patches9</metric>
    <metric>mean compet0</metric>
    <metric>mean compet1</metric>
    <metric>mean compet2</metric>
    <metric>mean compet3</metric>
    <metric>mean compet4</metric>
    <metric>mean compet5</metric>
    <metric>mean compet6</metric>
    <metric>mean compet7</metric>
    <metric>mean compet8</metric>
    <metric>mean compet9</metric>
    <metric>mean vis_before0</metric>
    <metric>mean vis_before1</metric>
    <metric>mean vis_before2</metric>
    <metric>mean vis_before3</metric>
    <metric>mean vis_before4</metric>
    <metric>mean vis_before5</metric>
    <metric>mean vis_before6</metric>
    <metric>mean vis_before7</metric>
    <metric>mean vis_before8</metric>
    <metric>mean vis_before9</metric>
    <metric>suc_juv</metric>
    <metric>juv_fail_spec</metric>
    <metric>feed_fail_spec0</metric>
    <metric>feed_fail_spec1</metric>
    <metric>feed_fail_spec2</metric>
    <metric>feed_fail_spec3</metric>
    <metric>feed_fail_spec4</metric>
    <metric>feed_fail_spec5</metric>
    <metric>feed_fail_spec6</metric>
    <metric>feed_fail_spec7</metric>
    <metric>feed_fail_spec8</metric>
    <metric>feed_fail_spec9</metric>
    <metric>failedpreg_spec</metric>
    <metric>failbold_spec</metric>
    <metric>failage_spec</metric>
    <metric>failedlact_spec</metric>
    <metric>pregs_spec</metric>
    <metric>mort_order</metric>
    <metric>mort_food</metric>
    <metric>[length spec-list] of patches with [habitat = 1]</metric>
    <metric>[sum [patchfeed] of patches with [habitat = 1] in-radius maxhr] of turtles</metric>
    <enumeratedValueSet variable="fragmentation_mode">
      <value value="&quot;frag&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clump">
      <value value="0.9999"/>
      <value value="0.999"/>
      <value value="0.99"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting_indivs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="feed-amount">
      <value value="16"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_frag_singlespec" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <metric>number</metric>
    <metric>mean ages</metric>
    <metric>sum [real-mass] of turtles</metric>
    <metric>mean [real-mass] of turtles</metric>
    <metric>mean (sentence rep_success_0 rep_success_1 rep_success_2 rep_success_3 rep_success_4 rep_success_5 rep_success_6 rep_success_7 rep_success_8 rep_success_9)</metric>
    <metric>mean (sentence rep_success_hr_0 rep_success_hr_1 rep_success_hr_2 rep_success_hr_3 rep_success_hr_4 rep_success_hr_5 rep_success_hr_6 rep_success_hr_7 rep_success_hr_8 rep_success_hr_9)</metric>
    <metric>mean (sentence life_juv_hr0 life_juv_hr1 life_juv_hr2 life_juv_hr3 life_juv_hr4 life_juv_hr5 life_juv_hr6 life_juv_hr7 life_juv_hr8 life_juv_hr9)</metric>
    <metric>mean (sentence life_juv0 life_juv1 life_juv2 life_juv3 life_juv4 life_juv5 life_juv6 life_juv7 life_juv8 life_juv9)</metric>
    <metric>mean hr</metric>
    <metric>variance hr</metric>
    <metric>mean maxhrs</metric>
    <metric>mean countmaxhr</metric>
    <metric>mean stors</metric>
    <metric>mean (sentence abs_stor0 abs_stor1 abs_stor2 abs_stor3 abs_stor4 abs_stor5 abs_stor6 abs_stor7 abs_stor8 abs_stor9)</metric>
    <metric>mean rm</metric>
    <metric>mean juv</metric>
    <metric>mean (sentence fmr0 fmr1 fmr2 fmr3 fmr4 fmr5 fmr6 fmr7 fmr8 fmr9)</metric>
    <metric>mean (sentence loco0 loco1 loco2 loco3 loco4 loco5 loco6 loco7 loco8 loco9)</metric>
    <metric>mean (sentence repro0 repro1 repro2 repro3 repro4 repro5 repro6 repro7 repro8 repro9)</metric>
    <metric>mean (sentence grow0 grow1 grow2 grow3 grow4 grow5 grow6 grow7 grow8 grow9)</metric>
    <metric>mean (sentence basal0 basal1 basal2 basal3 basal4 basal5 basal6 basal7 basal8 basal9)</metric>
    <metric>mean (sentence digest0 digest1 digest2 digest3 digest4 digest5 digest6 digest7 digest8 digest9)</metric>
    <metric>mean (sentence prod0 prod1 prod2 prod3 prod4 prod5 prod6 prod7 prod8 prod9)</metric>
    <metric>mean (sentence in0 in1 in2 in3 in4 in5 in6 in7 in8 in9)</metric>
    <metric>mean (sentence inpatch0 inpatch1 inpatch2 inpatch3 inpatch4 inpatch5 inpatch6 inpatch7 inpatch8 inpatch9)</metric>
    <metric>mean (sentence inhr0 inhr1 inhr2 inhr3 inhr4 inhr5 inhr6 inhr7 inhr8 inhr9)</metric>
    <metric>mean (sentence move0 move1 move2 move3 move4 move5 move6 move7 move8 move9)</metric>
    <metric>mean (sentence patches0 patches1 patches2 patches3 patches4 patches5 patches6 patches7 patches8 patches9)</metric>
    <metric>sum suc_juv</metric>
    <metric>sum juv_fail_spec</metric>
    <metric>sum (sentence feed_fail_spec0 feed_fail_spec1 feed_fail_spec2 feed_fail_spec3 feed_fail_spec4 feed_fail_spec5 feed_fail_spec6 feed_fail_spec7 feed_fail_spec8 feed_fail_spec9)</metric>
    <metric>sum failedpreg_spec</metric>
    <metric>sum failbold_spec</metric>
    <metric>sum failage_spec</metric>
    <metric>sum failedlact_spec</metric>
    <metric>sum pregs_spec</metric>
    <metric>mort_order</metric>
    <metric>mort_food</metric>
    <metric>[length spec-list] of patches with [habitat = 1]</metric>
    <metric>[sum [patchfeed] of patches with [habitat = 1] in-radius maxhr] of turtles</metric>
    <enumeratedValueSet variable="fragmentation_mode">
      <value value="&quot;frag&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clump">
      <value value="0.9999"/>
      <value value="0.999"/>
      <value value="0.99"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="specs-included">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
      <value value="8"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting_indivs">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="feed-amount">
      <value value="16"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
