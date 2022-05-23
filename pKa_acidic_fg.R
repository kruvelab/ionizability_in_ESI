library(rcdk)
library(ggplot2)
library(plotly)

setwd("C:/Users/chem5199/OneDrive - Kruvelab/Chimnaz/Data/ionization efficiencies/WP 2 IE predictions")

dataset = read_delim('pKa_logP_data_for_IE_data.csv',
                     delim = ",",
                     col_names = TRUE,
                     trim_ws = TRUE)



dataset = dataset[-c(397, 398, 399, 400),]

list <- as.vector(dataset$SMILES)

#acidic functional groups

function_is_benzoic_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_benzoic_acid <- 'c1ccc(cc1)C(=O)[O;H]'
  match = rcdk::matches(query_benzoic_acid, mols)
  return(match)
}
function_is_benzoic_acid = Vectorize(function_is_benzoic_acid)


function_is_phenol = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_phenol <- 'c1ccc(cc1)[O;H]'
  match = rcdk::matches(query_phenol, mols)
  return(match)
}
function_is_phenol = Vectorize(function_is_phenol)


function_is_methanol = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_methanol <- 'C[O;H]'
  match = rcdk::matches(query_methanol, mols)
  return(match)
}
function_is_methanol = Vectorize(function_is_methanol)





function_is_acetic_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_acetic_acid <- 'CC(=O)[O;H]'
  match = rcdk::matches(query_acetic_acid, mols)
  return(match)
}
function_is_acetic_acid = Vectorize(function_is_acetic_acid)


function_is_hydroxy_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_hydroxy_acid <- '[O;H]C(=O)c1ccccc1O'
  match = rcdk::matches(query_hydroxy_acid, mols)
  return(match)
}
function_is_hydroxy_acid = Vectorize(function_is_hydroxy_acid)


function_is_trichloroacetic_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_trichloroacetic_acid <- 'C(=O)(C(Cl)(Cl)Cl)[O;H]'
  match = rcdk::matches(query_trichloroacetic_acid, mols)
  return(match)
}
function_is_trichloroacetic_acid = Vectorize(function_is_trichloroacetic_acid)


function_is_sulfonic_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_sulfonic_acid <- 'S(=O)(=O)[O;H]'
  match = rcdk::matches(query_sulfonic_acid, mols)
  return(match)
}
function_is_sulfonic_acid = Vectorize(function_is_sulfonic_acid)


function_is_phosphate = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_phosphate <- 'P[O;H](O)=O'
  match = rcdk::matches(query_phosphate, mols)
  return(match)
}
function_is_phosphate = Vectorize(function_is_phosphate)


function_is_salicylic_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_salicylic_acid <- 'c1ccc(c(c1)C(=O)[O;H])O'
  match = rcdk::matches(query_salicylic_acid, mols)
  return(match)
}
function_is_salicylic_acid = Vectorize(function_is_salicylic_acid)


function_is_carboxylic_acid = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_carboxylic_acid <- '[O;H]C(=O)'
  match = rcdk::matches(query_carboxylic_acid, mols)
  return(match)
}
function_is_carboxylic_acid = Vectorize(function_is_carboxylic_acid)



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



function_is_benzimidazole= function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_benzimidazole <- 'c1ccc2[nH]cnc2c1'
  match = rcdk::matches(query_benzimidazole, mols)
  return(match)
}
function_is_benzimidazole = Vectorize(function_is_benzimidazole)


function_is_purine= function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_purine <- 'c1ncc2nc[nH]c2n1'
  match = rcdk::matches(query_purine, mols)
  return(match)
}
function_is_purine = Vectorize(function_is_purine)




function_is_pyrazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrazole <- 'c1cn[nH]c1'
  match = rcdk::matches(query_pyrazole, mols)
  return(match)
}
function_is_pyrazole = Vectorize(function_is_pyrazole)




function_is_pyrrole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_pyrrole <- 'c1cc[nH]c1'
  match = rcdk::matches(query_pyrrole, mols)
  return(match)
}
function_is_pyrrole = Vectorize(function_is_pyrrole)



function_is_carbazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_carbazole <- 'c1ccc2c(c1)[nH]c3ccccc23'
  match = rcdk::matches(query_carbazole, mols)
  return(match)
}
function_is_carbazole = Vectorize(function_is_carbazole)



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



function_is_thiazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_thiazole <- 'c1cscn1'
  match = rcdk::matches(query_thiazole, mols)
  return(match)
}
function_is_thiazole = Vectorize(function_is_thiazole)



function_is_sulfonamide = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_sulfonamide <- 'S(N)(=O)=O'
  match = rcdk::matches(query_sulfonamide, mols)
  return(match)
}
function_is_sulfonamide = Vectorize(function_is_sulfonamide)

function_is_saccharine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_saccharine <- 'c1ccc2c(c1)C(=NS2(=O)=O)O'
  match = rcdk::matches(query_saccharine, mols)
  return(match)
}
function_is_saccharine = Vectorize(function_is_saccharine)



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


function_is_3_hydroxypyridine = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_3_hydroxypyridine <- 'Oc1cccnc1'
  match = rcdk::matches(query_3_hydroxypyridine, mols)
  return(match)
}
function_is_3_hydroxypyridine = Vectorize(function_is_3_hydroxypyridine)



function_is_indazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_indazole <- 'c1ccc2c(c1)cn[nH]2'
  match = rcdk::matches(query_indazole, mols)
  return(match)
}
function_is_indazole = Vectorize(function_is_indazole)


function_is_6_nitroindazole = function(SMILES) {
  mols <- parse.smiles(SMILES)
  query_6_nitroindazole <- 'c1cc(cc2c1c[nH]n2)N(=O)=O'
  match = rcdk::matches(query_6_nitroindazole, mols)
  return(match)
}
function_is_6_nitroindazole = Vectorize(function_is_6_nitroindazole)



dataset = dataset %>%
  mutate(is_phenol = function_is_phenol(SMILES),
         is_benzoic_acid = function_is_benzoic_acid(SMILES),
         is_methanol = function_is_methanol(SMILES),
         is_acetic_acid = function_is_acetic_acid(SMILES),
         is_hydroxyacid = function_is_hydroxy_acid(SMILES),
         is_trichloroacetic_acid = function_is_trichloroacetic_acid(SMILES),
         is_sulfonic_acid = function_is_sulfonic_acid(SMILES),
         is_phosphate = function_is_phosphate(SMILES),
         is_salicylic_acid = function_is_salicylic_acid(SMILES),
         is_carboxylic_acid = function_is_carboxylic_acid(SMILES),
         is_sulfonamide = function_is_sulfonamide(SMILES),
         is_pyrrole = function_is_pyrrole(SMILES),
         is_carbazole = function_is_carbazole(SMILES),
         is_indole = function_is_indole(SMILES),
         is_thiophen = function_is_thiophene(SMILES),
         is_thiazole = function_is_thiazole(SMILES),
         is_imidazole = function_is_imidazole(SMILES),
         is_triazole123 = function_is_triazole123(SMILES),
         is_triazole124 = function_is_triazole124(SMILES),
         is_purine = function_is_purine(SMILES),
         is_pyrazole = function_is_pyrazole(SMILES),
         is_tetrazole = function_is_tetrazole(SMILES),
         is_saccharine = function_is_saccharine(SMILES),
         is_4_nitroimidazole = function_is_4_nitroimidazole(SMILES),
         is_2_nitroimidazole = function_is_2_nitroimidazole(SMILES),
         is_3_hydroxypyridine = function_is_3_hydroxypyridine(SMILES),
         is_indazole = function_is_indazole(SMILES),
         is_6_nitroindazole = function_is_6_nitroindazole(SMILES),
         is_benzimidazole = function_is_benzimidazole(SMILES)) %>%
  mutate(pKa_acid_pred = case_when(
    is_sulfonic_acid ~ -2.5,
    is_trichloroacetic_acid ~ 0.66,
    is_saccharine ~ 1.31,
    is_phosphate ~ 2.34,
    is_salicylic_acid ~ 2.79,
    is_hydroxyacid ~ 2.97,
    is_benzoic_acid ~ 4.2,
    is_acetic_acid ~ 4.76,
    is_tetrazole ~ 4.9,
    is_carboxylic_acid ~ 5,
    is_4_nitroimidazole ~ 8.31,
    is_2_nitroimidazole ~ 8.72,
    is_3_hydroxypyridine ~ 8.75,
    is_purine ~ 8.9,
    is_sulfonamide ~ 9.09,
    is_triazole123 ~ 9.3,
    is_phenol ~ 10,
    is_triazole124 ~ 10.3,
    is_benzimidazole ~ 12.3,
    is_6_nitroindazole ~ 13.61,
    is_indazole ~ 14.01,
    is_imidazole ~ 14.4,
    is_methanol ~ 15.5,
    is_pyrazole ~ 19.8,
    is_carbazole ~ 19.9,
    is_indole ~ 21,
    is_pyrrole ~ 23,
    is_thiazole ~ 29.4,
    is_thiophen ~ 33,
    TRUE ~ NA_real_
  ))


graph = ggplot(data = dataset) +
  geom_point(mapping = aes(x = pKa_acid,
                           y = pKa_acid_pred,
                           color = name)) +
  theme(legend.position="none") +
  ylim(-1,20)

ggplotly(graph)

source('my_theme_Anneli.R')

ggplot(data = dataset) +
  geom_point(mapping = aes(x =pKa_acid, y =pKa_acid_pred)) +
  theme(aspect.ratio = 1)+
  ylim(-1,20) +
  my_theme


sum <- lm(dataset$pKa_acid_pred ~ dataset$pKa_acid)
summary(sum)


# R^2 0.8872

dataset$pKa_acid_pred[is.na(dataset$pKa_acid_pred)] <- 20
dataset$pKa_acid[is.na(dataset$pKa_acid)] <- 20

RMSE(dataset$pKa_acid_pred, dataset$pKa_acid)
#3.491903


write_delim(dataset,
            "data_with_pka_acid.csv",
            delim = ",")
