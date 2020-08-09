#!/usr/bin/env python3
import pandas as pd 
import matplotlib.pyplot as plt
from JT_test import JT_test_pd
from hypersus_helper import norm_data_for_JT, estimate_lfc, summary_data

#parameters
pseudocount = 4
plots_on = False
run_calcs = True
out_fold = 'output/'

def main():
    #Initializations
    n_f = 0 #Figure index (adds one after each figure plot)

    #Read in and organized data, annot_df contains annotations of each site.
    full_df,full_annot_df = read_in_avium_data()

    #Drop the Tn samples
    full_df = full_df.drop(labels=['Tn'],level='abx',axis=1)

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
            plt.xlim(xmin=-10*1.1,xmax=10*1.1)
            plt.ylim(ymin=-10*1.1,ymax=10*1.1)
            plt.xlabel('Log2-RelativeFitness_'+namex)
            plt.ylabel('Log2-RelativeFitness_'+namey)
            plt.savefig(out_fold+'12h_'+namex+'_vs_'+namey+'.png')
            n_f+=1

            #48h
            plt.figure(n_f)
            plt.scatter(lfc_48h_df[namex],lfc_48h_df[namey])
            plt.xlim(xmin=-10*1.1,xmax=10*1.1)
            plt.ylim(ymin=-10*1.1,ymax=10*1.1)
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
    csv1 = '/home/will/baderlab/Projects/Tnseq_Avium/Results/Avium109_TnSeq_analysis/2018-06-11_TnSeq_analysis/Testing_Processing/2018-08-24_TPP_out1/LaneMerge_Pool_1.csv'
    csv2 = '/home/will/baderlab/Projects/Tnseq_Avium/Results/Avium109_TnSeq_analysis/2018-06-11_TnSeq_analysis/Testing_Processing/2018-08-24_TPP_out2/LaneMerge_Pool_2.csv'
    
    full_df = pd.concat([pd.read_csv(csv1,dtype={'regulatory_class':str,'bound_moiety':str}), pd.read_csv(csv2,dtype={'regulatory_class':str,'bound_moiety':str}).iloc[:,6:]],sort=False,axis=1)
    annot_df = full_df.iloc[:,0:6].rename(columns={'unique_identifier (locus_tag or record_id_start_end_strand)':'uid'})

    mindex = pd.MultiIndex.from_arrays([annot_df['contig'],annot_df['insertion_site']])
    full_df = full_df.drop(labels=['contig','insertion_site','unique_identifier (locus_tag or record_id_start_end_strand)','product','regulatory_class','bound_moiety'],axis=1)
    full_df.index = mindex
    annot_df=annot_df.iloc[:,2:]
    annot_df.index = mindex
    #full_df.to_csv(out_fold+'hypersus_testing.csv')

    labels_df = pd.read_excel("/home/will/baderlab/Projects/Tnseq_Avium/Results/Avium109_TnSeq_analysis/2018-06-11_TnSeq_analysis/Master_Labels.xls")

    #Read in data and rename for ease of reference
    samp_column_names = list()
    tuples = list()
    for i in range(len(labels_df["Index"])):
        old_name = 'read_count (' + labels_df["Run_UID"][i] + '_' + labels_df["Index"][i] + ')'
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


################################################################################
#The functions below are no longer used. Reference only.

#def reindex_by_uid(norm_df,annot_df):
#    #Changes the indices for the dataframe to use the uid (after removing sites with zero or multiple uids)
#
#    #Remove sites which don't contain any features:
#    norm_df = norm_df[~annot_df['uid'].isna()]
#    annot_df = annot_df[~annot_df['uid'].isna()]
#    
#    #Remove sites which contain more than one feature.
#    func = lambda x: x.find(';') == -1
#    norm_df = norm_df[annot_df['uid'].map(func)]
#    annot_df = annot_df[annot_df['uid'].map(func)]
#
#    #Compile into genes. #TODO: Make this into a function (and possibly combine with removal of multifeature sites and zero feature sites)
#    new_index = annot_df.index.to_frame()
#    new_index['uid'] = annot_df['uid']
#    new_index = new_index.drop('insertion_site',axis=1)
#    norm_df.index = new_index['uid']
#    annot_df.index = new_index['uid']
#
#    return norm_df,annot_df
#def lg2plus(x):
#    if x < math.pow(2.0,MINLOG):
#        return MINLOG
#    elif x > math.pow(2.0,MAXLOG):
#        return MAXLOG
#    else:
#        return math.log(x,2)
#
#def summary_data(norm_12h_rm0_df,mmfc_12h_df,treatment):
#    t0 = time.time()
#    df = norm_12h_rm0_df[treatment]
#    df_index = df.index
#    df = df.reset_index(drop=True)
#
#    df['pval'] = JT_test_pd(df,treatment,direction='lessthan',random_smooth=True)
#    df.index = df_index
#
#    p_ls = []
#    name_ls = []
#    for name,g in df.groupby(['uid'],sort=False):
#        name_ls.append(name)
#        Zscore,_ = stats.combine_pvalues(g['pval'],method='stouffer')
#        p_pooled = stats.norm.sf(abs(Zscore))*2
#        p_ls.append(p_pooled)
#    _,p_adj_ls,_,_ = smm.multipletests(p_ls, alpha=0.05, method='fdr_bh')
#    pval = pd.DataFrame({'uid':name_ls,'pval':p_ls,'pval-adj (BH)':p_adj_ls})
#
#    merged_df = mmfc_12h_df[treatment].merge(pval,how='left',left_on='uid',right_on='uid')
#    merged_df.index = merged_df['uid']
#    merged_df = merged_df.drop('uid',axis=1)
#    print(time.time()-t0)
#    return merged_df
#
#def summary_data_R(norm_12h_rm0_df,mmfc_12h_df,treatment,N_perm):
#    t0 = time.time()
#    df = norm_12h_rm0_df[treatment]
#    df_index = df.index
#    df = df.reset_index(drop=True)
#
#    #Define the groups used by JT calculation
#    g = []
#    treat_to_int = dict(zip(treatment,range(len(treatment)))) #Map treatment to integers for R code
#    for t in df.columns:
#        g.append(treat_to_int[t])
#
#    dfr = pandas2ri.py2ri(df.T)
#    out_dfr = gmwt.gmw(dfr,g,test='jt',nper=1,type='external',alternative='two.sided') #This calls to an R package
#    out_ls = list(out_dfr[0])
#    df['pval'] = out_ls
#    df.index = df_index
#
#    p_ls = []
#    name_ls = []
#    for name,g in df.groupby(['uid'],sort=False):
#        name_ls.append(name)
#        _,p_pooled = stats.combine_pvalues(g['pval'],method='fisher')
#        p_ls.append(p_pooled)
#    _,p_adj_ls,_,_ = smm.multipletests(p_ls, alpha=0.05, method='fdr_bh')
#    pval = pd.DataFrame({'uid':name_ls,'pval':p_ls,'pval-adj (BH)':p_adj_ls})
#
#    merged_df = mmfc_12h_df[treatment].merge(pval,how='left',left_on='uid',right_on='uid')
#    merged_df.index = merged_df['uid']
#    merged_df = merged_df.drop('uid',axis=1)
#    print(time.time()-t0)
#    return merged_df

#def summary_data(norm_12h_rm0_df,mmfc_12h_df,treatment,N_perm):
#    #This code is slow and only here for reference in case a faster implementation is attempted later.
#    C_12h_df = norm_12h_rm0_df[treatment]
#    p_ls = []
#    name_ls = []
#    for name,g in C_12h_df.groupby(['uid'],sort=False):
#        name_ls.append(name)
#        p_ls.append(pooled_JT_test(g,N_perm,direction='greaterthan'))
#    _,p_adj_ls,_,_ = smm.multipletests(p_ls, alpha=0.05, method='fdr_bh')
#    pval = pd.DataFrame({'uid':name_ls,'pval':p_ls,'pval-adj (BH)':p_adj_ls})
#    merged_df = mmfc_12h_df[treatment].merge(pval,how='left',left_on='uid',right_on='uid')
#    return merged_df
#
#def pooled_JT_test(df,N_perm,direction='lessthan'):
#    #Input Note: df is assumed to be already ordered so that the suspected smallest group is on the left and the suspected largest group is on the right.
#    i = 0
#    p_ls = []
#    group_ls = list(df.columns.values)
#    for i in range(len(df.index)):
#        data_ls = list(df.iloc[i,:])
#        p_ls.append(JT_test(data_ls,group_ls,N_perm,direction))
#
#    if not p_ls:
#        print(df)
#        exit('This case is not correctly handled yet, you need to test the behavior of smm.multipletests when NaNs are present')
#        p_pooled = np.nan
#    else:
#        _,p_pooled = stats.combine_pvalues(p_ls,method='fisher')
#
#    return p_pooled
#
#def JT_test(data_ls,group_ls,N_perm,direction='lessthan'):
#    #Input: data_ls is a list of values, group_ls is the group that each value belongs to, N_perm is the number of permutations that should be done.
#    if direction=='lessthan':
#        pass
#    elif direction=='greaterthan':
#        data_ls = data_ls[::-1]
#        group_ls = group_ls[::-1]
#    else:
#        exit('Direction: '+direction+' not recognized.')
#
#    #Check that group_ls is sorted
#    #Calculate JT statistic (B)
#    B = JT_statistic(data_ls,group_ls)
#    #Permute samples, build list of JT statistic
#    I_perm_ls = [None]*N_perm
#    for i in range(N_perm):
#        I_perm_ls[i] = (JT_statistic(np.random.permutation(data_ls),group_ls) >= B)
#
#    return sum(I_perm_ls)/N_perm
#
#def JT_statistic(data_ls,group_ls):
#    #Input: data_ls is a list of values, group_ls is the group that each value belongs to 
#    #Note: This technically implements Tapestra's statistic - though there is a linear transform to get Jonckheere's statistics. This won't affect the results of the permutation test however.
#    #Build a list-of-lists datastructure to make processing easier
#    t0 = time.time()
#    prev_g = group_ls[0]
#    curr_ls = []
#    lls = []
#    for i,g in enumerate(group_ls):
#        if g == prev_g:
#            curr_ls.append(data_ls[i])
#        else:
#            lls.append(curr_ls)
#            curr_ls = []
#            curr_ls.append(data_ls[i])
#            prev_g = g
#    lls.append(curr_ls)
#
#    k = len(set(group_ls))
#
#    #Calculate JT statistic
#    B = 0
#    for i in range(k-1):
#        for j in range(i+1,k):
#            for r in range(len(lls[i])):
#                for s in range(len(lls[j])):
#                    if lls[i][r] < lls[j][s]:
#                        D = 1
#                    elif lls[i][r] == lls[j][s]:
#                        D = 0.5
#                    else:
#                        D = 0
#                    B = B + D
#    t1 = time.time()
#    #print(t1-t0)
#    return B

if __name__ == '__main__':
    main()
