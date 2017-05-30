from __future__ import division
import math
from collections import defaultdict
import sys

def convert_to_str(list_val):
    return '_'.join(list_val)
        
def em(paired_words, max_iter = 100):
    
    #build prob_table    
    corsp_seq = defaultdict(lambda: defaultdict(set))
    for i, a in enumerate(alignments):
        for idx in xrange(len(a[0])):
            ep, jp = a[0][idx], a[1][idx]            
            corsp_seq[ep][jp].add(i)

    #initialization  
    prob_table = defaultdict(lambda: defaultdict(float))  
    for ep in corsp_seq:
        for jp in corsp_seq[ep]:            
            prob_table[ep][jp] = 1/len(corsp_seq[ep])    

    for iteration in xrange(max_iter):
        p_x_z = defaultdict(lambda:defaultdict(float))
        p_x = defaultdict(float)
        data_prob = 0.0

        #E-step        
        for aligned_pair in alignments:
            x =  convert_to_str(aligned_pair[0])
            z = convert_to_str(aligned_pair[1])
            (w1, w2) = aligned_pair
            
            p_x_z[x][z] = 1.0
            for i in xrange(len(w1)):
                    p_x_z[x][z] *= prob_table[w1[i]][w2[i]]
            p_x[x] += p_x_z[x][z]
            
        #renormalized to get fraccount
        for x in p_x_z:
            data_prob += math.log(p_x[x])
            for z in p_x_z[x]:
                p_x_z[x][z] = p_x_z[x][z]/p_x[x]
                                    

        #M-step
        p_ep_jp = defaultdict(lambda:defaultdict(float))
        p_ep = defaultdict(float)        
        for ep in prob_table:
            p_ep[ep] = 0.0
            for jp in prob_table[ep]:
                p_ep_jp[ep][jp] = 0.0
                for pair_idx in corsp_seq[ep][jp]:                    
                    aligned_pair = alignments[pair_idx]
                    x =convert_to_str(aligned_pair[0])
                    z =convert_to_str(aligned_pair[1])
                    p_ep_jp[ep][jp] += p_x_z[x][z]
                p_ep[ep] += p_ep_jp[ep][jp]
        for ep in prob_table:            
            for jp in prob_table[ep]:
                prob_table[ep][jp] = p_ep_jp[ep][jp]/p_ep[ep]        


        print 'iteration {0} ----- corpus prob = {1}'.format(iteration,data_prob)
        non_zeros = 0
        for ep in prob_table:
            m = ep+'|->    '
            for jp in prob_table[ep]:
                if prob_table[ep][jp] > 0.01:
                    non_zeros += 1
                    m += '{0}: {1:.2f}    '.format(jp, prob_table[ep][jp])
            print m
        print 'nonzeros = {0}'.format(non_zeros)
        if non_zeros == len(prob_table):
            break
    return prob_table

def extract_epron_jpron_data(lines):
    pairs = []
    is_paired = True
    ln_num = 0
    while ln_num < len(lines):        
        l1 = lines[ln_num].strip().split()
        l2 = lines[ln_num+1].strip().split()
        ln_num += 3        
        pairs += [(l1,l2)]        
    return pairs

if __name__ == '__main__':
    arg = sys.argv
    print arg
    if len(arg) > 1:
        max_iter = int(arg[1])
    
    fin = sys.stdin
    lines = fin.readlines()

    pairs_list = extract_epron_jpron_data(lines)    
    prob_table = em(pairs_list, max_iter)
    with open('epron_jpron.prob', 'w') as fp:
        for ep in prob_table:
            for jp in prob_table[ep]:
                if prob_table[ep][jp] > 0.01:                    
                    fp.write('{0} : {1} # {2}\n'.format(ep, jp, prob_table[ep][jp]))