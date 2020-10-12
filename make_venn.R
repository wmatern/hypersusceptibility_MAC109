setwd("/run/media/will/phd_data1/ForDistribution_hypersusceptibility_MAC109")

library('eulerr')
#library('nVennR') #Alternative package for displaying - but less pretty

## Plot Venn of Hypersusceptible mutants
clar_df = read.csv('output/hypersus_analysis/TableS3_CLAR.csv', row.names='Genomic.Feature')
clar_sus_set = row.names(clar_df[clar_df['Prediction'] == 'Hypersusceptible',])

moxi_df = read.csv('output/hypersus_analysis/TableS4_MOXI.csv', row.names='Genomic.Feature')
moxi_sus_set = row.names(moxi_df[moxi_df['Prediction'] == 'Hypersusceptible',])

rfb_df = read.csv('output/hypersus_analysis/TableS5_RFB.csv', row.names='Genomic.Feature')
rfb_sus_set = row.names(rfb_df[rfb_df['Prediction'] == 'Hypersusceptible',])

emb_df = read.csv('output/hypersus_analysis/TableS6_EMB.csv', row.names='Genomic.Feature')
emb_sus_set = row.names(emb_df[emb_df['Prediction'] == 'Hypersusceptible',])

emb_sus_set1 = setdiff(intersect(emb_sus_set, rfb_sus_set), union(moxi_sus_set, clar_sus_set))
emb_sus_set2 = setdiff(emb_sus_set, emb_sus_set1)

pdf('output/hypersus_analysis/hypersusceptible_venn.pdf')
plot(euler(list(Clarithromycin = clar_sus_set, Moxifloxacin = moxi_sus_set,
                 Rifabutin = rfb_sus_set, "Ethambutol (1)" = emb_sus_set1, "Ethambutol (2)" = emb_sus_set2), 
           shape='ellipse'), fills= c('#808080','cyan','red','yellow','yellow'),
quantities = list(type = c("counts")))
dev.off()

# myV1 = plotVenn(list(Clarithromycin = clar_sus_set, Moxifloxacin = moxi_sus_set, 
#                      Rifabutin = rfb_sus_set, Ethambutol = emb_sus_set), 
#                 setColors=c('#808080','cyan','red','yellow'), labelRegions = F, 
#                 nCycles=5000, opacity = 0.1, borderWidth = 3, outFile='output/hypersus_analysis/nVennR.svg')

# Output union of hypersus mutants for each drug
sus_set = sort(union(union(clar_sus_set, moxi_sus_set), union(rfb_sus_set, emb_sus_set)))

df_sus_merge = data.frame(row.names = sus_set)
df_sus_merge[['Genbank.Annotation']] = clar_df[sus_set,'Genbank.Annotation']
df_sus_merge[['Clarithromycin_Prediction']] = clar_df[sus_set,'Prediction']
df_sus_merge[['Moxifloxacin_Prediction']] = moxi_df[sus_set,'Prediction']
df_sus_merge[['Rifabutin_Prediction']] = rfb_df[sus_set,'Prediction']
df_sus_merge[['Ethambutol_Prediction']] = emb_df[sus_set,'Prediction']

write.csv(df_sus_merge, file='output/hypersus_analysis/TableS1_hypersusceptible.csv')

## Plot Venn of hypertolerant mutants
clar_df = read.csv('output/hypersus_analysis/TableS3_CLAR.csv', row.names='Genomic.Feature')
clar_tol_set = row.names(clar_df[clar_df['Prediction'] == 'Hypertolerant',])

moxi_df = read.csv('output/hypersus_analysis/TableS4_MOXI.csv', row.names='Genomic.Feature')
moxi_tol_set = row.names(moxi_df[moxi_df['Prediction'] == 'Hypertolerant',])

rfb_df = read.csv('output/hypersus_analysis/TableS5_RFB.csv', row.names='Genomic.Feature')
rfb_tol_set = row.names(rfb_df[rfb_df['Prediction'] == 'Hypertolerant',])

emb_df = read.csv('output/hypersus_analysis/TableS6_EMB.csv', row.names='Genomic.Feature')
emb_tol_set = row.names(emb_df[emb_df['Prediction'] == 'Hypertolerant',])

pdf('output/hypersus_analysis/hypertolerant_venn.pdf')
#Warning message caused by next line due to empty set (EMB):
plot(euler(list(Clarithromycin = clar_tol_set, Moxifloxacin = moxi_tol_set,
                Rifabutin = rfb_tol_set, Ethambutol = emb_tol_set), shape='ellipse'),
     fills= c('#808080','cyan','red','yellow'), quantities = list(type = c("counts")))
dev.off()

# Output union of hypertol mutants for each drug
tol_set = sort(union(union(clar_tol_set, moxi_tol_set), union(rfb_tol_set, emb_tol_set)))

df_tol_merge = data.frame(row.names = tol_set)
df_tol_merge[['Genbank.Annotation']] = clar_df[tol_set,'Genbank.Annotation']
df_tol_merge[['Clarithromycin_Prediction']] = clar_df[tol_set,'Prediction']
df_tol_merge[['Moxifloxacin_Prediction']] = moxi_df[tol_set,'Prediction']
df_tol_merge[['Rifabutin_Prediction']] = rfb_df[tol_set,'Prediction']
df_tol_merge[['Ethambutol_Prediction']] = emb_df[tol_set,'Prediction']

write.csv(df_tol_merge, file='output/hypersus_analysis/TableS2_hypertolerant.csv')
