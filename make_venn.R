setwd("/run/media/will/phd_data1/ForDistribution_hypersusceptibility_MAC109")

library('eulerr')
#library('nVennR') #Alternative package for displaying - but less pretty

## Plot Venn of Hypersusceptible mutants
clar_df = read.csv('output/hypersus_analysis/TableS1_CLAR.csv')
clar_sus_set = as.vector(clar_df[clar_df['Prediction'] == 'Hypersusceptible','Genomic.Feature'])

moxi_df = read.csv('output/hypersus_analysis/TableS2_MOXI.csv')
moxi_sus_set = as.vector(moxi_df[moxi_df['Prediction'] == 'Hypersusceptible','Genomic.Feature'])

rfb_df = read.csv('output/hypersus_analysis/TableS3_RFB.csv')
rfb_sus_set = rfb_df[rfb_df['Prediction'] == 'Hypersusceptible','Genomic.Feature']

emb_df = read.csv('output/hypersus_analysis/TableS4_EMB.csv')
emb_sus_set = emb_df[emb_df['Prediction'] == 'Hypersusceptible','Genomic.Feature']

emb_sus_set1 = setdiff(intersect(emb_sus_set, rfb_sus_set), union(moxi_sus_set, clar_sus_set))
emb_sus_set2 = setdiff(emb_sus_set, emb_sus_set1)

plot(euler(list(Clarithromycin = clar_sus_set, Moxifloxacin = moxi_sus_set,
                 Rifabutin = rfb_sus_set, "Ethambutol (1)" = emb_sus_set1, "Ethambutol (2)" = emb_sus_set2), 
           shape='ellipse'), fills= c('#808080','cyan','red','yellow','yellow'),
quantities = list(type = c("counts")))

# myV1 = plotVenn(list(Clarithromycin = clar_sus_set, Moxifloxacin = moxi_sus_set, 
#                      Rifabutin = rfb_sus_set, Ethambutol = emb_sus_set), 
#                 setColors=c('#808080','cyan','red','yellow'), labelRegions = F, 
#                 nCycles=5000, opacity = 0.1, borderWidth = 3, outFile='output/hypersus_analysis/nVennR.svg')

## Plot Venn of hypertolerant mutants
clar_df = read.csv('output/hypersus_analysis/TableS1_CLAR.csv')
clar_tol_set = as.vector(clar_df[clar_df['Prediction'] == 'Hypertolerant','Genomic.Feature'])

moxi_df = read.csv('output/hypersus_analysis/TableS2_MOXI.csv')
moxi_tol_set = as.vector(moxi_df[moxi_df['Prediction'] == 'Hypertolerant','Genomic.Feature'])

rfb_df = read.csv('output/hypersus_analysis/TableS3_RFB.csv')
rfb_tol_set = rfb_df[rfb_df['Prediction'] == 'Hypertolerant','Genomic.Feature']

emb_df = read.csv('output/hypersus_analysis/TableS4_EMB.csv')
emb_tol_set = emb_df[emb_df['Prediction'] == 'Hypertolerant','Genomic.Feature']

plot(euler(list(Clarithromycin = clar_tol_set, Moxifloxacin = moxi_tol_set,
                Rifabutin = rfb_tol_set, Ethambutol = emb_tol_set), shape='ellipse'),
     fills= c('#808080','cyan','red','yellow'), quantities = list(type = c("counts")))

