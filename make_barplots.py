#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

in_fold = 'output/hypersus_analysis/'
out_fold = 'output/hypersus_analysis/'

#CLAR
fname = in_fold+'TableS3_CLAR.csv'

df = pd.read_csv(fname)[['Genomic Feature','Clarithromycin (5.4ug/mL) - LFC_48h', 'Prediction']]
df = df[df['Prediction'] != 'Not Significant'] #Do not plot effect sizes for not significant features
df = df.rename({'Clarithromycin (5.4ug/mL) - LFC_48h':'LFC'},axis='columns')

plt.subplots()
df = df.sort_values(by='LFC',ascending=True)
print('Number of significant hypersusceptible mutants (CLAR):',sum((df['Prediction'] == 'Hypersusceptible')))
print('Number of significant hypertolerant mutants (CLAR):',sum((df['Prediction'] == 'Hypertolerant')))
print(df.iloc[0:10])
plt.bar(df['Genomic Feature'],df['LFC'],1,edgecolor=None,color='#808080')
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change) (48 hours, 5.4ug/mL Clarithromycin)')
plt.axhline(y=0,color='k')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.savefig(out_fold+'CLAR_barplot.svg',bbox_inches='tight')

#RFB
fname = in_fold+'TableS5_RFB.csv'

df = pd.read_csv(fname)[['Genomic Feature','Rifabutin (0.63ug/mL) - LFC_48h', 'Prediction']]
df = df[df['Prediction'] != 'Not Significant'] #Do not plot effect sizes for not significant features
df = df.rename({'Rifabutin (0.63ug/mL) - LFC_48h':'LFC'},axis='columns')

plt.subplots()
df = df.sort_values(by='LFC',ascending=True)
print('Number of significant hypersusceptible mutants (RFB):',sum((df['Prediction'] == 'Hypersusceptible')))
print('Number of significant hypertolerant mutants (RFB):',sum((df['Prediction'] == 'Hypertolerant')))
print(df.iloc[0:10])
plt.bar(df['Genomic Feature'],df['LFC'],1,edgecolor=None,color='r')
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change) (48 hours, 0.63ug/mL Rifabutin)')
plt.axhline(y=0,color='k')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.savefig(out_fold+'RFB_barplot.svg',bbox_inches='tight')

#MOXI
fname = in_fold+'TableS4_MOXI.csv'

df = pd.read_csv(fname)[['Genomic Feature','Moxifloxacin (1.0ug/mL) - LFC_48h', 'Prediction']]
df = df[df['Prediction'] != 'Not Significant'] #Do not plot effect sizes for not significant features
df = df.rename({'Moxifloxacin (1.0ug/mL) - LFC_48h':'LFC'},axis='columns')

plt.subplots()
df = df.sort_values(by='LFC',ascending=True)
print('Number of significant hypersusceptible mutants (MOXI):',sum((df['Prediction'] == 'Hypersusceptible')))
print('Number of significant hypertolerant mutants (MOXI):',sum((df['Prediction'] == 'Hypertolerant')))
print(df.iloc[0:10])
plt.bar(df['Genomic Feature'],df['LFC'],1,edgecolor=None,color='c')
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change) (48 hours, 1.0ug/mL Moxifloxacin)')
plt.axhline(y=0,color='k')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.savefig(out_fold+'MOXI_barplot.svg',bbox_inches='tight')

#EMB
fname = in_fold+'TableS6_EMB.csv'

df = pd.read_csv(fname)[['Genomic Feature','Ethambutol (2.1ug/mL) - LFC_48h', 'Prediction']]
df = df[df['Prediction'] != 'Not Significant'] #Do not plot effect sizes for not significant features
df = df.rename({'Ethambutol (2.1ug/mL) - LFC_48h':'LFC'},axis='columns')

plt.subplots()
df = df.sort_values(by='LFC',ascending=True)
print('Number of significant hypersusceptible mutants (EMB):',sum((df['Prediction'] == 'Hypersusceptible')))
print('Number of significant hypertolerant mutants (EMB):',sum((df['Prediction'] == 'Hypertolerant')))
print(df.iloc[0:10])
plt.bar(df['Genomic Feature'],df['LFC'],1,edgecolor='k',color='y')
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change) (48 hours, 2.1ug/mL Ethambutol)')
plt.axhline(y=0,color='k')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.savefig(out_fold+'EMB_barplot.svg',bbox_inches='tight')
