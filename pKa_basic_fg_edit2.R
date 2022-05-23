library(rcdk)
library(ggplot2)
library(plotly)
library(tidyverse)


setwd("C:/Users/chem5199/OneDrive - Kruvelab/Chimnaz/Data/ionization efficiencies/WP 2 IE predictions")


dataset = read_delim('pKa_logP_data_for_IE_data.csv',
                     delim = ",",
                     col_names = TRUE,
                     trim_ws = TRUE)



dataset = dataset[-c(397, 398, 399, 400),]



list <- as.vector(dataset$SMILES)


#dataset$pKa_acid[is.na(dataset$pKa_acid)] <- 20
#dataset$pKa_base[is.na(dataset$pKa_base)] <- -5


#basic functional groups


function_is_pyridine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyridine <- 'c1ccncc1'
  match = rcdk::matches(query_pyridine, mols)
  return(match)
}
function_is_pyridine = Vectorize(function_is_pyridine)

function_is_amine_primary = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_amine_primary <- '[C;H2][N;H2]'
  query_amine_primary1 <- '[N;H2][C;H1][C]'
  query_amine_primary2 <- '[N;H2][C;C&C]'
  match = rcdk::matches(query_amine_primary, mols)
  match1 = rcdk::matches(query_amine_primary1, mols)
  match2 = rcdk::matches(query_amine_primary2, mols)
  if (match | match1 | match2) {
    match = TRUE
  }
  return(match)
}
function_is_amine_primary = Vectorize(function_is_amine_primary)


#function_is_amine_primary('CN')

function_is_amine_second = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_amine_second <- '[C;H2]~[N;H1]~[C;H2]'
  query_amine_second1 <- '[C;H2]~[N;H1]~[C;C&H1]'
  query_amine_second2 <- '[C;H2][N;H1][C;C&C]'
  query_amine_second3 <- '[C;C&H1][N;H1][C;H2]'
  query_amine_second4 <- '[C;C&H1][N;H1][C;C&H1]'
  query_amine_second5 <- '[C;C&H1][N;H1][C;C&C]'
  query_amine_second6 <- '[C;C&C][N;H1][C;H2]'
  query_amine_second7 <- '[C;C&C][N;H1][C;C&H1]'
  query_amine_second8 <- '[C;C&C][N;H1][C;C&C]'
  match = rcdk::matches(query_amine_second, mols)
  match1 = rcdk::matches(query_amine_second1, mols)
  match2 = rcdk::matches(query_amine_second2, mols)
  match3 = rcdk::matches(query_amine_second3, mols)
  match4 = rcdk::matches(query_amine_second4, mols)
  match5 = rcdk::matches(query_amine_second5, mols)
  match6 = rcdk::matches(query_amine_second6, mols)
  match7 = rcdk::matches(query_amine_second7, mols)
  match8 = rcdk::matches(query_amine_second8, mols)
  if (match | match1 | match2 | match3 | match4 | match5 | match6 | match7 | match8) {
    match = TRUE
  }
  return(match)
}
function_is_amine_second = Vectorize(function_is_amine_second)




function_is_amine_tert = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_amine_tert <- '[C;H2][N]([C;H2])[C;H2]'
  query_amine_tert1 <- '[C;H2][N]([C;H2])[C;C&H1]'
  query_amine_tert2 <- '[C;H2][N]([C;H2])[C;C&C]'
  query_amine_tert3 <- '[C;C&H1][N]([C;H2])[C;H2]'
  query_amine_tert4 <- '[C;C&H1][N]([C;H2])[C;C&H1]'
  query_amine_tert5 <- '[C;C&H1][N]([C;H2])[C;C&C]'
  query_amine_tert6 <- '[C;C&C][N]([C;H2])[C;H2]'
  query_amine_tert7 <- '[C;C&C][N]([C;H2])[C;C&H1]'
  query_amine_tert8 <- '[C;C&C][N]([C;H2])[C;C&C]'
  query_amine_tert9 <- '[C;H2]~[N]([C;C&H1])~[C;H2]'
  query_amine_tert10 <- '[C;H2]~[N]([C;C&H1])~[C;C&H1]'
  query_amine_tert11 <- '[C;H2][N]([C;C&H1])[C;C&C]'
  query_amine_tert12 <- '[C;C&H1][N]([C;C&H1])[C;H2]'
  query_amine_tert13 <- '[C;C&H1][N]([C;C&H1])[C;C&H1]'
  query_amine_tert14 <- '[C;C&H1][N]([C;C&H1])[C;C&C]'
  query_amine_tert15 <- '[C;C&C][N]([C;C&H1])[C;H2]'
  query_amine_tert16 <- '[C;C&C][N]([C;C&H1])[C;C&H1]'
  query_amine_tert17 <- '[C;C&C][N]([C;C&H1])[C;C&C]'
  query_amine_tert18 <- '[C;H2]~[N]([C;C&C])~[C;H2]'
  query_amine_tert19 <- '[C;H2]~[N]([C;C&C])~[C;C&H1]'
  query_amine_tert20 <- '[C;H2][N]([C;C&C])[C;C&C]'
  query_amine_tert21 <- '[C;C&H1][N]([C;C&C])[C;H2]'
  query_amine_tert22 <- '[C;C&H1][N]([C;C&C])[C;C&H1]'
  query_amine_tert23 <- '[C;C&H1][N]([C;C&C])[C;C&C]'
  query_amine_tert24 <- '[C;C&C][N]([C;C&C])[C;H2]'
  query_amine_tert25 <- '[C;C&C][N]([C;C&C])[C;C&H1]'
  query_amine_tert26 <- 'CC(C)N(C(C)C)C(C)C '
  match = rcdk::matches(query_amine_tert, mols)
  match1 = rcdk::matches(query_amine_tert1, mols)
  match2 = rcdk::matches(query_amine_tert2, mols)
  match3 = rcdk::matches(query_amine_tert3, mols)
  match4 = rcdk::matches(query_amine_tert4, mols)
  match5 = rcdk::matches(query_amine_tert5, mols)
  match6 = rcdk::matches(query_amine_tert6, mols)
  match7 = rcdk::matches(query_amine_tert7, mols)
  match8 = rcdk::matches(query_amine_tert8, mols)
  match9 = rcdk::matches(query_amine_tert9, mols)
  match10 = rcdk::matches(query_amine_tert10, mols)
  match11 = rcdk::matches(query_amine_tert11, mols)
  match12 = rcdk::matches(query_amine_tert12, mols)
  match13 = rcdk::matches(query_amine_tert13, mols)
  match14 = rcdk::matches(query_amine_tert14, mols)
  match15 = rcdk::matches(query_amine_tert15, mols)
  match16 = rcdk::matches(query_amine_tert16, mols)
  match17 = rcdk::matches(query_amine_tert17, mols)
  match18 = rcdk::matches(query_amine_tert18, mols)
  match19 = rcdk::matches(query_amine_tert19, mols)
  match20 = rcdk::matches(query_amine_tert20, mols)
  match21 = rcdk::matches(query_amine_tert21, mols)
  match22 = rcdk::matches(query_amine_tert22, mols)
  match23 = rcdk::matches(query_amine_tert23, mols)
  match24 = rcdk::matches(query_amine_tert24, mols)
  match25 = rcdk::matches(query_amine_tert25, mols)
  match26 = rcdk::matches(query_amine_tert26, mols)
  if (match | match1 | match2 | match3 | match4 | match5 | match6 | match7 | match8 | match9 | match10 | match11 | match12 | match13 | match14 | match15 | match16 | match17 | match18 | match19 | match20 | match21 | match22 | match23 | match24 | match25 | match26) {
    match = TRUE
  }
  return(match)
}
function_is_amine_tert = Vectorize(function_is_amine_tert)


#checking

# function_is_amine_tert('CN=C(O)O/N=C(/C(=O)N(C)C)SC')
# 
# mols_x <- parse.smiles(c('CC(C)N(C(C)C)C(C)C'))
# query_x <- '[C;C&C][N]([C;C&C])[C;C&C]'
# rcdk::matches(query_x, mols_x)
# 


function_is_aniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_aniline <- 'c1ccc(cc1)N'
  match = rcdk::matches(query_aniline, mols)
  return(match)
}
function_is_aniline = Vectorize(function_is_aniline)



function_is_quinoline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_quinoline <- 'c1ccc2c(c1)cccn2'
  match = rcdk::matches(query_quinoline, mols)
  return(match)
}
function_is_quinoline= Vectorize(function_is_quinoline)


#aromatic amines



function_is_methyl_aniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_methyl_aniline <- 'Cc1ccc(N)cc1'
  match = rcdk::matches(query_methyl_aniline, mols)
  return(match)
}
function_is_methyl_aniline= Vectorize(function_is_methyl_aniline)


function_is_chloroaniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_chloroaniline <- 'Nc1ccc(Cl)cc1'
  match = rcdk::matches(query_chloroaniline, mols)
  return(match)
}
function_is_chloroaniline= Vectorize(function_is_chloroaniline)



function_is_4_nitroaniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_4_nitroaniline <- 'c1cc(ccc1N)N(=O)=O'
  match = rcdk::matches(query_4_nitroaniline, mols)
  return(match)
}
function_is_4_nitroaniline= Vectorize(function_is_4_nitroaniline)

#function_is_4_nitroaniline('c1cc(c(cc1Cl)N(=O)=O)N')


function_is_3_nitroaniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_3_nitroaniline <- 'c1cc(cc(c1)N(=O)=O)N'
  match = rcdk::matches(query_3_nitroaniline, mols)
  return(match)
}
function_is_3_nitroaniline= Vectorize(function_is_3_nitroaniline)


function_is_2_nitroaniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_2_nitroaniline <- 'Nc1ccccc1N(=O)=O'
  match = rcdk::matches(query_2_nitroaniline, mols)
  return(match)
}
function_is_2_nitroaniline= Vectorize(function_is_2_nitroaniline)


	

function_is_imidazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_imidazole <- 'c1c[nH]cn1'
  match = rcdk::matches(query_imidazole, mols)
  return(match)
}
function_is_imidazole= Vectorize(function_is_imidazole)


function_is_4_nitroimidazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_4_nitroimidazole <- 'c1c([nH]cn1)N(=O)=O'
  match = rcdk::matches(query_4_nitroimidazole, mols)
  return(match)
}
function_is_4_nitroimidazole= Vectorize(function_is_4_nitroimidazole)



function_is_2_nitroimidazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_2_nitroimidazole <- 'c1c[nH]c(n1)N(=O)=O'
  match = rcdk::matches(query_2_nitroimidazole, mols)
  return(match)
}
function_is_2_nitroimidazole= Vectorize(function_is_2_nitroimidazole)



function_is_benzimidazole= function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_benzimidazole <- 'c1ccc2[nH]cnc2c1'
  match = rcdk::matches(query_benzimidazole, mols)
  return(match)
}
function_is_benzimidazole = Vectorize(function_is_benzimidazole)


function_is_piperidine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_piperidine <- 'C1CCNCC1'
  match = rcdk::matches(query_piperidine, mols)
  return(match)
}
function_is_piperidine= Vectorize(function_is_piperidine)


function_is_triazole123 = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_triazole123 <- 'c1c[nH]nn1'
  match = rcdk::matches(query_triazole123, mols)
  return(match)
}
function_is_triazole123 = Vectorize(function_is_triazole123)


function_is_triazole124 = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_triazole124 <- 'c1nc[nH]n1'
  match = rcdk::matches(query_triazole124, mols)
  return(match)
}
function_is_triazole124 = Vectorize(function_is_triazole124)


function_is_tetrazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_tetrazole <- 'c1nnn[nH]1'
  match = rcdk::matches(query_tetrazole, mols)
  return(match)
}
function_is_tetrazole = Vectorize(function_is_tetrazole)




function_is_purine= function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_purine <- 'c1ncc2nc[nH]c2n1'
  match = rcdk::matches(query_purine, mols)
  return(match)
}
function_is_purine = Vectorize(function_is_purine)


function_is_morpholine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_morpholine <- 'C1COCCN1'
  match = rcdk::matches(query_morpholine, mols)
  return(match)
}
function_is_morpholine = Vectorize(function_is_morpholine)


function_is_pyrazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrazole <- 'c1cn[nH]c1'
  match = rcdk::matches(query_pyrazole, mols)
  return(match)
}
function_is_pyrazole = Vectorize(function_is_pyrazole)


function_is_pyrimidine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrimidine <- 'c1cncnc1'
  match = rcdk::matches(query_pyrimidine, mols)
  return(match)
}
function_is_pyrimidine = Vectorize(function_is_pyrimidine)


function_is_pyrazine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrazine <- 'c1cnccn1'
  match = rcdk::matches(query_pyrazine, mols)
  return(match)
}
function_is_pyrazine = Vectorize(function_is_pyrazine)


function_is_pyridazine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyridazine <- 'c1ccnnc1'
  match = rcdk::matches(query_pyridazine, mols)
  return(match)
}
function_is_pyridazine = Vectorize(function_is_pyridazine)




function_is_pyrrole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrrole <- 'c1cc[nH]c1'
  match = rcdk::matches(query_pyrrole, mols)
  return(match)
}
function_is_pyrrole = Vectorize(function_is_pyrrole)


function_is_pyrrolidine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrrolidine <- 'C1CCNC1'
  match = rcdk::matches(query_pyrrolidine, mols)
  return(match)
}
function_is_pyrrolidine = Vectorize(function_is_pyrrolidine)



function_is_carbazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_carbazole <- 'c1ccc2c(c1)[nH]c3ccccc23'
  match = rcdk::matches(query_carbazole, mols)
  return(match)
}
function_is_carbazole = Vectorize(function_is_carbazole)


function_is_indoline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_indoline <- 'C1Cc2ccccc2N1'
  match = rcdk::matches(query_indoline, mols)
  return(match)
}
function_is_indoline = Vectorize(function_is_indoline)



function_is_indole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_indole <- 'c1ccc2[nH]ccc2c1'
  match = rcdk::matches(query_indole, mols)
  return(match)
}
function_is_indole = Vectorize(function_is_indole)


function_is_thiophene = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_thiophene <- 'c1ccsc1'
  match = rcdk::matches(query_thiophene, mols)
  return(match)
}
function_is_thiophene = Vectorize(function_is_thiophene)



function_is_oxazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_oxazole <- 'c1cocn1'
  match = rcdk::matches(query_oxazole, mols)
  return(match)
}
function_is_oxazole = Vectorize(function_is_oxazole)


function_is_thiazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_thiazole <- 'c1cscn1'
  match = rcdk::matches(query_thiazole, mols)
  return(match)
}
function_is_thiazole = Vectorize(function_is_thiazole)


function_is_sulfadiazine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_sulfadiazine <- 'Nc1ccc(cc1)S(=O)(=O)Nc2ncccn2'
  match = rcdk::matches(query_sulfadiazine, mols)
  return(match)
}
function_is_sulfadiazine = Vectorize(function_is_sulfadiazine)


function_is_sulfonamide = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_sulfonamide <- 'S(N)(=O)=O'
  match = rcdk::matches(query_sulfonamide, mols)
  return(match)
}
function_is_sulfonamide = Vectorize(function_is_sulfonamide)



function_is_acridine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_acridine <- 'c1ccc2nc3ccccc3cc2c1'
  match = rcdk::matches(query_acridine, mols)
  return(match)
}
function_is_acridine = Vectorize(function_is_acridine)




function_is_amide = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_amide <- 'C(=N)O'
  match = rcdk::matches(query_amide, mols)
  return(match)
}
function_is_amide = Vectorize(function_is_amide)



function_is_amide2 = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_amide2 <- 'C(N)=O'
  match = rcdk::matches(query_amide2, mols)
  return(match)
}
function_is_amide2 = Vectorize(function_is_amide2)


function_is_carbamate = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_carbamate <- 'N=C(O)O'
  match = rcdk::matches(query_carbamate, mols)
  return(match)
}
function_is_carbamate = Vectorize(function_is_carbamate)

function_is_carbamate2 = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_carbamate2 <- 'NC(=O)O'
  match = rcdk::matches(query_carbamate2, mols)
  return(match)
}
function_is_carbamate2 = Vectorize(function_is_carbamate2)


function_is_3_hydroxypyridine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_3_hydroxypyridine <- 'Oc1cccnc1'
  match = rcdk::matches(query_3_hydroxypyridine, mols)
  return(match)
}
function_is_3_hydroxypyridine = Vectorize(function_is_3_hydroxypyridine)


function_is_2_chloropyridine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_2_chloropyridine <- 'Clc1ccccn1'
  match = rcdk::matches(query_2_chloropyridine, mols)
  return(match)
}
function_is_2_chloropyridine = Vectorize(function_is_2_chloropyridine)


function_is_25_chloroaniline = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_25_chloroaniline <- 'Nc1cc(Cl)ccc1Cl'
  match = rcdk::matches(query_25_chloroaniline, mols)
  return(match)
}
function_is_25_chloroaniline = Vectorize(function_is_25_chloroaniline)



#--------------------

dataset = dataset %>%
  mutate(is_pyridine = function_is_pyridine(SMILES),
         is_2_chloropyridine = function_is_2_chloropyridine(SMILES),
         is_amine_primary = function_is_amine_primary(SMILES),
         is_quinoline = function_is_quinoline(SMILES),
         is_amine_tert = function_is_amine_tert(SMILES),
         is_amine_second = function_is_amine_second(SMILES),
         is_aniline = function_is_aniline(SMILES),
         is_methylaniline = function_is_methyl_aniline(SMILES),
         is_chloroaniline = function_is_chloroaniline(SMILES),
         is_4_nitroaniline = function_is_4_nitroaniline(SMILES),
         is_3_nitroaniline = function_is_3_nitroaniline(SMILES),
         is_imidazole = function_is_imidazole(SMILES),
         is_piperidine = function_is_piperidine(SMILES),
         is_triazole124 = function_is_triazole124(SMILES),
         is_morpholine = function_is_morpholine(SMILES),
         is_pyrazine = function_is_pyrazine(SMILES),
         is_pyrrolidine = function_is_pyrrolidine(SMILES),
         is_indoline = function_is_indoline(SMILES),
         is_oxazole = function_is_oxazole(SMILES),
         is_thiazole = function_is_thiazole(SMILES),
         is_purine = function_is_purine(SMILES),
         is_sulfonamide = function_is_sulfonamide(SMILES),
         is_pyrimidine = function_is_pyrimidine(SMILES),
         is_pyridazine = function_is_pyridazine(SMILES),
         is_amide = function_is_amide(SMILES),
         is_amide2 = function_is_amide2(SMILES),
         is_2_nitroaniline = function_is_2_nitroaniline(SMILES),
         is_carbamate = function_is_carbamate(SMILES),
         is_carbamate2 = function_is_carbamate2(SMILES),
         is_4_nitroimidazole = function_is_4_nitroimidazole(SMILES),
         is_2_nitroimidazole = function_is_2_nitroimidazole(SMILES),
         is_benzimidazole = function_is_benzimidazole(SMILES),
         is_3_hydroxypyridine = function_is_3_hydroxypyridine(SMILES),
         is_25_chloroaniline =function_is_25_chloroaniline(SMILES),
         is_acridine = function_is_acridine(SMILES)) %>%
  mutate(pKa_base_pred = case_when(
    is_pyrrolidine ~ 11.3,
    is_piperidine ~ 11,
    is_amine_second ~ 10.73,
    is_amine_primary ~ 10.64,
    is_amine_tert ~ 9.81,
    is_morpholine ~ 8.4,
    is_imidazole ~ 6.95,
    is_benzimidazole ~ 5.79,
    is_acridine ~ 5.6,
    is_pyridine ~ 5.2,
    is_methylaniline ~ 5.08,
    is_quinoline ~ 4.92,
    is_indoline ~ 4.9,
    is_3_hydroxypyridine ~ 4.79,
    is_aniline ~ 4.63, 
    is_chloroaniline ~ 4.15,
    is_thiazole ~ 2.5,
    is_purine ~ 2.5,
    is_3_nitroaniline ~ 2.47,
    is_pyridazine ~ 2.3,
    is_triazole124 ~ 2.2,
    is_25_chloroaniline ~ 1.57,
    is_pyrimidine ~ 1.3,
    is_4_nitroaniline ~ 1,
    is_oxazole ~ 0.8,
    is_pyrazine ~ 0.6,
    is_2_chloropyridine ~ 0.49,
    is_2_nitroaniline ~ -0.28,
    is_4_nitroimidazole ~ -0.3,
    is_2_nitroimidazole ~ -0.5,
    is_carbamate ~ -2.11,
    is_carbamate2 ~ -2.11,
    is_sulfonamide ~ -2.7,
    is_amide ~ -5,
    is_amide2 ~ -5,
    TRUE ~ NA_real_
  ))
for(i in 1:nrow(dataset)) {
  if (dataset[i, 38]==TRUE) {
    dataset[i, 41]<- 4.79
  } else if (dataset[i, 16]==TRUE) {
    dataset[i, 41]<- 2.47
  } else  if (dataset[i, 15]==TRUE) {
    dataset[i, 41]<- 1
  } else if (dataset[i, 7]==TRUE) {
    dataset[i, 41]<- 0.49
  } else if (dataset[i, 32]==TRUE) {
    dataset[i, 41]<- -0.28
  } else if (dataset[i, 35]==TRUE) {
    dataset[i, 41]<- -0.3
  } else if (dataset[i, 36]==TRUE) {
    dataset[i, 41]<- -0.5
  }
}


graph = ggplot(data = dataset) +
  geom_point(mapping = aes(x = pKa_base,
                           y = pKa_base_pred,
                           color = name)) +
  theme(legend.position="none")

ggplotly(graph)


source('my_theme_Anneli.R')
ggplot(data = dataset) +
  geom_point(mapping = aes(x =pKa_base, y =pKa_base_pred)) +
  theme(aspect.ratio = 1)+
  my_theme


sum <- lm(dataset$pKa_base_pred ~ dataset$pKa_base)
summary(sum)


# R^2 0.5921

dataset$pKa_base_pred[is.na(dataset$pKa_base_pred)] <- -5
dataset$pKa_base[is.na(dataset$pKa_base)] <- -5

RMSE(dataset$pKa_base_pred, dataset$pKa_base)
#5.153146

write_delim(dataset,
            "data_with_pka_base2.csv",
            delim = ",")



#checking




#mols_y <- parse.smiles(c('c1cc(ccc1N)N(=O)=O'))
#query_y <- 'c1ccc(cc1)N'
#rcdk::matches(query_y, mols_y)

#mols_x <- parse.smiles(c('C(CC(=N)O)[C@@H](C(=O)O)N'))
#query_x <- '[N;H2][CC][C]'
#query_x <- '[C;H1C1][N;H2]'
#query_x <- '[C;C2][N;H2]'


#query_x <- 'C(=O)N'
#query_x <- 'C(=N)O'
#query_x <- 'C(N)=O'
#rcdk::matches(query_x, mols_x)


#mols_y <- parse.smiles(c('c1cc(c(CC(=NC(=N)N)O)c(c1)Cl)Cl'))
#query_y <- 'Clc1ccccc1'
#rcdk::matches(query_y, mols_y)

#function_is_nitro = function(SMILES) {
#mols <- parse.smiles(SMILES)
#query_nitro <- 'N(=O)=O'
#match = rcdk::matches(query_nitro, mols)
#return(match)
#}
#function_is_nitro = Vectorize(function_is_nitro)

#function_is_nitro2 = function(SMILES) {
 # mols <- parse.smiles(SMILES)
  #query_nitro2 <- '[N+][O-])=O'
  #match = rcdk::matches(query_nitro2, mols)
  #return(match)
#}
#function_is_nitro2 = Vectorize(function_is_nitro2)


# function_is_oxamyl1 = function(SMILES) {
#   mols <- parse.smiles(SMILES)
#   query_oxamyl1 <- 'CNC(=O)ON=C(/SC)C(=O)N(C)C'
#   match = rcdk::matches(query_oxamyl1, mols)
#   return(match)
# }
# function_is_oxamyl1 = Vectorize(function_is_oxamyl1)
# 
# 
# function_is_oxamyl2 = function(SMILES) {
#   mols <- parse.smiles(SMILES)
#   query_oxamyl2 <- 'CN=C(O)O/N=C(/C(=O)N(C)C)SC'
#   match = rcdk::matches(query_oxamyl2, mols)
#   return(match)
# }
# function_is_oxamyl2 = Vectorize(function_is_oxamyl2)
