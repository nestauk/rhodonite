#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 17:07:42 2019

@author: jdjumalieva
"""

import pandas as pd
import os
import pickle
#import igraph
import numpy as np
import collections
import gensim
from sklearn.feature_extraction.text import CountVectorizer
import json
from ast import literal_eval
import re
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster.hierarchy import ward, dendrogram
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster 
import prep_for_clustering as prep
import scipy
from datasketch import MinHashLSHEnsemble, MinHash, MinHashLSH
import datetime
import networkx as nx

output_dir = '/Users/jdjumalieva/ESCoE/outputs'


def find_infrequent(df_col, threshold = 1):
    '''
    Identify skills that were mentioned fewer than threshold.
    Input is a column in a dataframe.
    In current application is used on a corpus of one year of adverts.
    '''
    count_model = CountVectorizer(tokenizer=prep.tokenize_asis, 
                        lowercase = False,
                        ngram_range=(1,1))
    X = count_model.fit_transform(df_col)
    names = count_model.get_feature_names()
    Xd = X.todense()
    sums = np.sum(Xd, axis = 0)
    not_common = np.where(sums <= threshold)[1].tolist()
    not_common_skills = [names[elem] for elem in not_common]
    return not_common_skills
    

def get_unique_sets(some_df, lower_thresh, upper_thresh, count_thresh):
    '''
    Get unique combinations of skills across all job adverts in a dataframe.
    Lower and upper thresholds stand for min and max allowable lenght of a skill set.
    Count threshold is not used at the moment, but can be used to filter out
    skill sets that occur fewer than specified number of times.
    '''
    stringed = some_df['clean_skills'].apply(lambda x: ', '.join(x))
    stringed = stringed.sort_values()
    unique_sets = pd.DataFrame(stringed.groupby(stringed).count())
    unique_sets.columns = ['count']
    unique_sets = unique_sets[unique_sets['count'] > count_thresh]
    unique_sets = unique_sets.reset_index()
    unique_sets['skill_list'] = unique_sets['clean_skills'].apply(lambda x:\
               x.split(', '))
    unique_sets['size'] = unique_sets['skill_list'].apply(lambda x: len(x))
    unique_sets = unique_sets[(unique_sets['size'] >lower_thresh) & \
                              (unique_sets['size'] <= upper_thresh)]
    
    return unique_sets


def get_skill_sets(some_set, skill_sets):
    '''
    Convert groupings of skill set IDs produced by LSH back into groupings of 
    skill sets.
    '''
    ix = some_set.split(' ')[1]
    res = skill_sets[int(ix)]
    return res


def union_find(data):
    '''Create Disjoint Data Structure from list of lists'''
    parents = {}
    def find(i):
        j = parents.get(i, i)
        if j == i:
            return i
        k = find(j)
        if k != j:
            parents[i] = k
        return k
    for l in filter(None, data):
        parents.update(dict.fromkeys(map(find, l), find(l[0])))
    merged = {}
    for k, v in parents.items():
        merged.setdefault(find(v), []).append(k)
    return list(merged.values())


def group_skill_sets_lsh(skill_sets, lsh_threshold = 0.8):
    '''
    Identify candidates within skill sets that likely have a Jaccard similarity
    above the specified threshold.
    '''
    #Hashing
    hash_objects = []
    for i in range(len(skill_sets)):
        m = MinHash(num_perm=200)
        hash_objects.append(m)
        
    for ix, skill_set in enumerate(skill_sets):
        for t in skill_set:
            hash_objects[ix].update(t.encode('utf8'))
    
    #Create LSH index
    #Jaccard threshold has to be specified at initiation
    lsh = MinHashLSH(threshold= lsh_threshold, num_perm=200)
    
    for ix, (skill_set, hash_object) in enumerate(zip(skill_sets,hash_objects)):
        skill_set_name = 'skill_set ' + str(ix)
        lsh.insert(skill_set_name, hash_object) # s[ix]
    
    #Query LSH for each unique skill set
    #t1 = datetime.datetime.now()
    candidates = []
    for ix, skill_set in enumerate(skill_sets):
        result = lsh.query(hash_objects[ix])
        candidates.append(result)
    #        print(result)    
    #    print('***************')    
    return candidates


#Read in 66 highly transversal skills
transversal = pd.read_csv(os.path.join(output_dir, 
                                       'top_transversal_skills.csv'),
                          encoding = 'utf-8')

transversal_skills = list(transversal['Skill'])
most_frequent = ['software engineering',
                 'software development']
    
#Identify communities of skills from raw online job adverts
rhodonite_sample = {}
filelist = os.listdir(os.path.join(output_dir, 'BG_by_year'))
filelist.pop(filelist.index('2012_2016_jobs_long.csv'))
filelist.pop(filelist.index('.DS_Store'))

sorted_filelist = sorted(filelist, key = lambda x: x[:4])

#Specify skill clusters in software engineering branch
soft_eng_clusters = ['software development', 'web development', 'data engineering',
               'servers and middleware', 'app development']


#rhodonite_ads_soft_eng = {}
rhodonite_communities_soft_eng = {}


for file in sorted_filelist:
    year = file[:4]
    print('Starting to process ', year)
    df = pd.read_csv(os.path.join(output_dir, 'BG_by_year', file), 
                     index_col = 0)
    sftwr_df = df[df['category'].isin(soft_eng_clusters)] #646,739
    sftwr_df['clean_skills'] = sftwr_df['clean_skills'].apply(lambda x: eval(x))
    
    #Filter out infrequent and transversal skills
    infrequent_skills = find_infrequent(sftwr_df['clean_skills'],
                                        threshold = 3)
    sftwr_df['clean_skills'] = sftwr_df['clean_skills'].apply(lambda x:\
    [elem.rstrip() for elem in x if elem not in transversal_skills +\
     infrequent_skills +\
     most_frequent])
    sftwr_df['clean_skills'] = sftwr_df['clean_skills'].apply(lambda x:\
        sorted(x))
#    rhodonite_ads_soft_eng[year] = [elem for elem in list(sftwr_df['clean_skills'])\
#                    if len(elem) >0]

    #Find disjoint sets of skills (using LSH with Jaccard similarity at 0.8)
    unique_sets = get_unique_sets(sftwr_df, 4, 20, 0) #175,462
    skill_sets = [set(elem) for elem in unique_sets['skill_list'].values]
    #Simplied iteration from 0.9 to 0.55 Jaccard similarity
    jaccard_range = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55]
    set_asides_processed = []
    for jr in jaccard_range:
        print('Jaccard Similarity ', jr)
    #    Get initial groupings of skill sets at a given jaccard similarity threshold
        candidates = group_skill_sets_lsh(skill_sets, lsh_threshold = jr)
    #    Identify sets that occur in many groupings of skill sets
    #    These are filtered out to prevent creating extremely large groupings
    #    at a disjoint set stage
        flat_candidates = (item for sublist in candidates for item in sublist)
        skill_set_count = collections.Counter(flat_candidates)
        sorted_set_count = sorted(skill_set_count.items(), key = lambda x: x[1],
                              reverse = True)
    #    Calculate thresholds for removing highly common skill sets
        outliers = np.percentile(list(skill_set_count.values()), 95) + 1
        potential_outliers = [k for k,v in sorted_set_count if v >= outliers]
    #    Separate groupings of highly common skill sets from others
    #    These are dealt with separately
        to_join = []
        set_aside = []
        outlier_set = set(potential_outliers)
        for ix, c in enumerate(candidates):
    #        print(ix)
            set_c = set(c)
            if len(set_c.intersection(outlier_set)) > 0:
                set_aside.append(c)
            else:
                to_join.append(c)    
    
    #   Find disjoint sets among groupings in to_join list
        disjoint_candidates = union_find(to_join)
    
        disjoint_skill_sets = [[get_skill_sets(skill_set, skill_sets) for skill_set\
                            in disjoint_candidate] for disjoint_candidate in\
                            disjoint_candidates if len(disjoint_candidate) <100]    
            
        communities = [set.union(*disjoint_sets) for disjoint_sets in \
                       disjoint_skill_sets if len(disjoint_sets) < 50]
        
    #    Deal with skill sets that were set aside
        
        set_aside_tups = map(tuple, set_aside)
        set_aside_tup_counts = collections.Counter(set_aside_tups)
    
        set_aside_skill_sets = []
        for k in set_aside_tup_counts.keys():
            set_aside_skill_sets.append([get_skill_sets(skill_set, skill_sets) for \
                                         skill_set in k])
    
        set_aside_communities = [set.union(*set_aside_sets) for set_aside_sets in \
                       set_aside_skill_sets]
        set_aside_candidates = group_skill_sets_lsh(set_aside_communities, \
                                                    lsh_threshold = 0.55)
        
        disjoint_set_asides = union_find(set_aside_candidates)
        disjoint_set_asides_skills = [[get_skill_sets(skill_set, set_aside_communities) for skill_set\
                            in disjoint_candidate] for disjoint_candidate in\
                            disjoint_set_asides if len(disjoint_candidate) <100]    
            
        set_aside_communities2 =[set.union(*set_aside_sets) for set_aside_sets in \
                       disjoint_set_asides_skills]
        
        set_asides_processed.append(set_aside_communities2)
        print(len(communities), len(disjoint_set_asides), len(set_aside_communities2))
        skill_sets = communities
    
    for set_aside in set_asides_processed:
        communities.append(set_aside)
    
    all_communities = [elem for elem in communities if len(elem) < 50]
   
    rhodonite_communities_soft_eng[year] = all_communities
    
    print(year, len(sftwr_df), len(communities))


#with open(os.path.join(output_dir, 'rhodonite_6y_ads_soft_eng.pkl'), 'wb') as f:
#    pickle.dump(rhodonite_ads_soft_eng, f)

with open(os.path.join(output_dir, 'rhodonite_6y_comms_soft_eng.pkl'), 'wb') as f:
    pickle.dump(rhodonite_communities_soft_eng, f)
    
#
#2012 387144
#2013 416886
#2014 321041
#2015 539573
#2016 609833
#2017 646739


##Test with one year
#test_file = sorted_filelist[-1:][0]
#year = test_file[:4]
#df = pd.read_csv(os.path.join(output_dir, 'BG_by_year', test_file), 
#                     index_col = 0) #8,907,937
#sftwr_df = df[df['category'].isin(soft_eng_clusters)] #646,739
#sftwr_df['clean_skills'] = sftwr_df['clean_skills'].apply(lambda x: eval(x))
#infrequent_skills = find_infrequent(sftwr_df['clean_skills'], threshold = 3)
#sftwr_df['clean_skills'] = sftwr_df['clean_skills'].apply(lambda x:\
#[elem.rstrip() for elem in x if elem not in transversal_skills + \
# infrequent_skills + \
# most_frequent
# ])
#sftwr_df['clean_skills'] = sftwr_df['clean_skills'].apply(lambda x:\
#    sorted(x))
#sftwr_df['size'] = sftwr_df['clean_skills'].apply(lambda x: len(x))
#sftwr_df = sftwr_df[sftwr_df['size'] >1]
##sftwr_df['num_skills'] = sftwr_df['clean_skills'].apply(lambda x: len(x))
#unique_sets = get_unique_sets(sftwr_df, 4, 20, 0) #175,462
#skill_sets = [set(elem) for elem in unique_sets['skill_list'].values]
#
#
##Simplied iteration from 0.9 to 0.55 Jaccard similarity
#jaccard_range = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55]
#set_asides_processed = []
#for jr in jaccard_range:
#    print('Jaccard Similarity ', jr)
##    Get initial groupings of skill sets at a given jaccard similarity threshold
#    candidates = group_skill_sets_lsh(skill_sets, lsh_threshold = jr)
##    Identify sets that occur in many groupings of skill sets
##    These are filtered out to prevent creating extremely large groupings
##    at a disjoint set stage
#    flat_candidates = (item for sublist in candidates for item in sublist)
#    skill_set_count = collections.Counter(flat_candidates)
#    sorted_set_count = sorted(skill_set_count.items(), key = lambda x: x[1],
#                          reverse = True)
##    Calculate thresholds for removing highly common skill sets
#    outliers = np.percentile(list(skill_set_count.values()), 95) + 1
#    potential_outliers = [k for k,v in sorted_set_count if v >= outliers]
##    Separate groupings of highly common skill sets from others
##    These are dealt with separately
#    to_join = []
#    set_aside = []
#    outlier_set = set(potential_outliers)
#    for ix, c in enumerate(candidates):
##        print(ix)
#        set_c = set(c)
#        if len(set_c.intersection(outlier_set)) > 0:
#            set_aside.append(c)
#        else:
#            to_join.append(c)    
#
##   Find disjoint sets among groupings in to_join list
#    disjoint_candidates = union_find(to_join)
#
#    disjoint_skill_sets = [[get_skill_sets(skill_set, skill_sets) for skill_set\
#                        in disjoint_candidate] for disjoint_candidate in\
#                        disjoint_candidates if len(disjoint_candidate) <100]    
#        
#    communities = [set.union(*disjoint_sets) for disjoint_sets in \
#                   disjoint_skill_sets if len(disjoint_sets) < 50]
#    
##    Deal with skill sets that were set aside
#    
#    set_aside_tups = map(tuple, set_aside)
#    set_aside_tup_counts = collections.Counter(set_aside_tups)
#
#    set_aside_skill_sets = []
#    for k in set_aside_tup_counts.keys():
#        set_aside_skill_sets.append([get_skill_sets(skill_set, skill_sets) for \
#                                     skill_set in k])
#
#    set_aside_communities = [set.union(*set_aside_sets) for set_aside_sets in \
#                   set_aside_skill_sets]
#    set_aside_candidates = group_skill_sets_lsh(set_aside_communities, \
#                                                lsh_threshold = 0.55)
#    
#    disjoint_set_asides = union_find(set_aside_candidates)
#    disjoint_set_asides_skills = [[get_skill_sets(skill_set, set_aside_communities) for skill_set\
#                        in disjoint_candidate] for disjoint_candidate in\
#                        disjoint_set_asides if len(disjoint_candidate) <100]    
#        
#    set_aside_communities2 =[set.union(*set_aside_sets) for set_aside_sets in \
#                   disjoint_set_asides_skills]
#    
#    set_asides_processed.append(set_aside_communities2)
#    print(len(communities), len(disjoint_set_asides), len(set_aside_communities2))
#    skill_sets = communities
#
#for set_aside in set_asides_processed:
#    communities.append(set_aside)
#
#all_communities = [elem for elem in communities if len(elem) < 50]
    
