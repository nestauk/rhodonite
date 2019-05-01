import numpy as np
from collections import Counter
import seaborn as sns
import scipy
from datasketch import MinHashLSHEnsemble, MinHash, MinHashLSH
from rhodonite.utilities import flatten, filter_subsets

def get_unique_groups(groups):
    '''
    Get unique combinations of skills across all job adverts in a dataframe.
    Lower and upper thresholds stand for min and max allowable lenght of a skill set.
    Count threshold is not used at the moment, but can be used to filter out
    skill sets that occur fewer than specified number of times.
    '''
    set_groups = [set(g) for g in groups]
    unique_groups = set(set_groups)
    return unique_groups

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

def groups_lsh(groups, lsh_threshold=0.8):
    '''
    Identify candidates within skill sets that likely have a Jaccard similarity
    above the specified threshold.
    '''
    #Hashing
    hash_objects = []
    for i in range(len(groups)):
        m = MinHash(num_perm=200)
        hash_objects.append(m)
        
    for ix, group in enumerate(groups):
        for t in group:
            hash_objects[ix].update(t.encode('utf8'))
    
    #Create LSH index
    lsh = MinHashLSH(threshold=lsh_threshold, num_perm=200)
    
    for ix, (group, hash_object) in enumerate(zip(groups, hash_objects)):
        group_name = 'group ' + str(ix)
        lsh.insert(group_name, hash_object) # s[ix]
    
    #Query LSH for each unique skill set
    #t1 = datetime.datetime.now()
    candidates = []
    for ix, group in enumerate(groups):
        result = lsh.query(hash_objects[ix])
        candidates.append(result)
    #        print(result)    
    #    print('***************')    
    return candidates

def get_group_sets(x, groups):
    '''
    Convert groupings of skill set IDs produced by LSH back into groupings of 
    skill sets.
    '''
    i = x.split(' ')[1]
    return groups[int(i)]

def coalesce_groups(groups, j_max=0.99, j_min=0.5, steps=10, commonality_threshold=95, max_size=100,
        remove_subsets=True):

    j_range = np.arange(j_min, j_max, (j_max - j_min)/steps)[::-1]
    set_asides_processed = []
    for j in j_range:
        candidates = groups_lsh(groups, lsh_threshold=j)
        group_set_counts = Counter(flatten(candidates))
        
        commonality_percentile = np.percentile(list(group_set_counts.values()), commonality_threshold) + 1
        commons = [k for k, v in group_set_counts.items() if v >= commonality_percentile]

        to_join = []
        set_aside = []

        common_set = set(commons)
        for i, c in enumerate(candidates):
            c = set(c)
            if len(c.intersection(common_set)) > 0:
                set_aside.append(tuple(c))
            else:
                to_join.append(tuple(c))

        disjoint_candidates = union_find(to_join)

        disjoint_group_sets = [[get_group_sets(group, groups) for group in
            disjoint_candidate] for disjoint_candidate in 
            disjoint_candidates if len(disjoint_candidate) <= max_size]

        communities = [set.union(*[set(d) for d in disjoint_sets]) for disjoint_sets in
                disjoint_group_sets if len(disjoint_sets) <= (max_size / 2)]

        set_aside_process = map(tuple, set_aside)
        set_aside_counts = Counter(set_aside_process)

        set_aside_groups = []

        for k in set_aside_counts.keys():
            set_aside_groups.append([get_group_sets(group, groups) for group in k])

        set_aside_communities = [set.union(*(set(s) for s in set_aside_group)) for 
                set_aside_group in set_aside_groups]
        set_aside_candidates = groups_lsh(set_aside_communities, lsh_threshold=j_min)

        disjoint_set_asides = union_find(set_aside_candidates)

        disjoint_set_aside_groups = [[get_group_sets(group, set_aside_communities) for group in
            disjoint_candidate] for disjoint_candidate in
            disjoint_set_asides if len(disjoint_candidates) <= max_size]

        set_aside_communities_end = [set.union(*set_aside_sets) for set_aside_sets in
                disjoint_set_aside_groups]

        set_asides_processed.append(set_aside_communities_end)
        groups = communities

    for sap in set_asides_processed:
        communities.append(sap)

    all_communities = [e for e in communities if len(e) <= (max_size / 2)]

    if remove_subsets:
        all_communities = filter_subsets((list(l) for l in all_communities))

    return all_communities

