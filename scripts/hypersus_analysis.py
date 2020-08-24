#!/usr/bin/env python3
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from JT_test import JT_test_pd
from hypersus_helper import norm_data_for_JT, estimate_lfc, summary_data

#parameters
pseudocount = 4
plots_on = True
run_calcs = True
out_fold = 'output/hypersus_analysis/'

#Seed the random number generator - ensure consistent results each time
np.random.seed(525600)

def main():
    #Initializations
    n_f = 0 #Figure index (adds one after each figure plot)

    #Read in and organized data, annot_df contains annotations of each site.
    full_df,full_annot_df = read_in_avium_data()

    #Drop the Tn samples (not necessary for final - you've already removed these)
    #full_df = full_df.drop(labels=['Tn'],level='abx',axis=1)

    ##Isolate single gene names in uid_df (for merging back after processing)
    uid_df = full_annot_df[['uid','product']].drop_duplicates(subset='uid') #Gene by gene annotations
    uid_df = uid_df[~uid_df['uid'].isna()]
    func = lambda x: x.find(';') == -1
    uid_df = uid_df[uid_df['uid'].map(func)]
    uid_df = uid_df.set_index('uid')

    lowrm_df = full_df[(full_df['0h']['In'].sum(axis=1) > 10)] #Remove sites with low read counts in the input
    annot_df = full_annot_df[(full_df['0h']['In'].sum(axis=1) > 10)]

    #Normalize to sum of the reads in each sample. #TODO: For samples with loss of diversity you should use a different normalizer - some estimate of the wildtype. 
    norm_df = full_df / full_df.sum(axis=0)

    lfc_12h_df = estimate_lfc(lowrm_df['12h'],annot_df,pseudocount,ctrl_grp = 'DMSO')
    norm_12h_df = norm_data_for_JT(lowrm_df['12h'],annot_df)

    lfc_48h_df = estimate_lfc(lowrm_df['48h'],annot_df,pseudocount,ctrl_grp = 'DMSO')
    norm_48h_df = norm_data_for_JT(lowrm_df['48h'],annot_df)

    if plots_on:
        for namex,namey in [('R++','R+++'),('M+','M++'),('C++','C+++'),('E+','E++')]:
            #12h
            plt.figure(n_f)
            plt.scatter(lfc_12h_df[namex],lfc_12h_df[namey])
            plt.xlim(xmin=-5*1.1,xmax=5*1.1)
            plt.ylim(ymin=-5*1.1,ymax=5*1.1)
            plt.xlabel('Log2-RelativeFitness_'+namex)
            plt.ylabel('Log2-RelativeFitness_'+namey)
            plt.savefig(out_fold+'12h_'+namex+'_vs_'+namey+'.png')
            n_f+=1

            #48h
            plt.figure(n_f)
            plt.scatter(lfc_48h_df[namex],lfc_48h_df[namey])
            plt.xlim(xmin=-5*1.1,xmax=5*1.1)
            plt.ylim(ymin=-5*1.1,ymax=5*1.1)
            plt.xlabel('Log2-RelativeFitness_'+namex)
            plt.ylabel('Log2-RelativeFitness_'+namey)
            plt.savefig(out_fold+'48h_'+namex+'_vs_'+namey+'.png')
            n_f+=1

    #Plot histogram of log fold change
    #plt.hist(lfc_48h_df['C++'].dropna(),bins=50, color='r')
    #plt.show()

    if run_calcs:
        # Computations using 12h timepoint only
        treatment = ['DMSO','C++','C+++']
        #merged_df = summary_data_R(norm_12h_df,lfc_12h_df,treatment,N_perm)
        merged_df = summary_data(norm_12h_df,lfc_12h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_12h_CLAR.csv')

        treatment = ['DMSO','M+','M++']
        #merged_df = summary_data_R(norm_12h_df,lfc_12h_df,treatment,N_perm)
        merged_df = summary_data(norm_12h_df,lfc_12h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_12h_MOXI.csv')
    
        treatment = ['DMSO','R++','R+++']
        #merged_df = summary_data_R(norm_12h_df,lfc_12h_df,treatment,N_perm)
        merged_df = summary_data(norm_12h_df,lfc_12h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_12h_RFB.csv')
    
        treatment = ['DMSO','E+','E++']
        #merged_df = summary_data_R(norm_12h_df,lfc_12h_df,treatment,N_perm)
        merged_df = summary_data(norm_12h_df,lfc_12h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_12h_EMB.csv')

        #Using 48h time point
        treatment = ['DMSO','C++','C+++']
        #merged_df = summary_data_R(norm_48h_df,lfc_48h_df,treatment,N_perm)
        merged_df = summary_data(norm_48h_df,lfc_48h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_48h_CLAR.csv')
    
        treatment = ['DMSO','M+','M++']
        #merged_df = summary_data_R(norm_48h_df,lfc_48h_df,treatment,N_perm)
        merged_df = summary_data(norm_48h_df,lfc_48h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_48h_MOXI.csv')
    
        treatment = ['DMSO','R++','R+++']
        #merged_df = summary_data_R(norm_48h_df,lfc_48h_df,treatment,N_perm)
        merged_df = summary_data(norm_48h_df,lfc_48h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_48h_RFB.csv')
    
        treatment = ['DMSO','E+','E++']
        #merged_df = summary_data_R(norm_48h_df,lfc_48h_df,treatment,N_perm)
        merged_df = summary_data(norm_48h_df,lfc_48h_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df.to_csv(out_fold+'summary_48h_EMB.csv')

def read_in_avium_data():
    #Read in and organized data. Return a dataframe containing organized raw read count data.

    #Read in csv files of counts
    csv1 = 'output/input.csv'
    
    full_df = pd.concat([pd.read_csv(csv1,dtype={'regulatory_class':str,'bound_moiety':str})],sort=False,axis=1)
    annot_df = full_df.iloc[:,0:6].rename(columns={'unique_identifier (locus_tag or record_id_start_end_strand)':'uid'})

    mindex = pd.MultiIndex.from_arrays([annot_df['contig'],annot_df['insertion_site']])
    full_df = full_df.drop(labels=['contig','insertion_site','unique_identifier (locus_tag or record_id_start_end_strand)','product','regulatory_class','bound_moiety'],axis=1)
    full_df.index = mindex
    annot_df=annot_df.iloc[:,2:]
    annot_df.index = mindex
    #full_df.to_csv(out_fold+'hypersus_testing.csv')

    labels_df = pd.read_csv("input/Sample_Labels.csv")
    
    #Read in data and rename for ease of reference
    samp_column_names = list()
    tuples = list()
    for i in range(len(labels_df["Index"])):
        old_name = 'read_count (' + labels_df["SRA Accession"][i] + ')'
        if labels_df["ID"][i][0:2] == 'In':
            new_name = labels_df["ID"][i]
            tuples.append(('0h','In'))
        elif labels_df["ID"][i][0:2] == 'Tn':
            #Remove these columns
            new_name = labels_df["ID"][i]
            tuples.append(('0h','Tn'))
        elif (labels_df["ID"][i][0:2] == '12') or (labels_df["ID"][i][0:2] == '48'):
            new_name = labels_df["Label2"][i] + '_' + labels_df["ID"][i][0:3]
            tuples.append((labels_df["ID"][i][0:3],labels_df["Label2"][i][:-2]))
        else:
            print(labels_df["ID"][i])
            exit("Unexpected ID")
        samp_column_names.append(old_name)
    
    #Build a multi-Index for the columns - You want to be able to refer to each group: (12h,48h) and (M+,C++...)
    full_df = full_df[samp_column_names]
    mindex = pd.MultiIndex.from_tuples(tuples, names=['time','abx'])
    full_df.columns=mindex

    #Sort the column names for convenience
    full_df = full_df.sort_index(axis=1,level=['time','abx'])
    #pd.concat([annot_df,full_df],join='inner',axis=1,levels=['time','abx']).to_csv(out_fold+'Organized_raw_data.csv')

    return(full_df, annot_df)

if __name__ == '__main__':
    main()
