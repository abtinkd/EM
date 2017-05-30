from __future__ import division
import math
from collections import defaultdict
import sys


def align(rem_short, rem_long, limit):
    if rem_short == 1:
        return [rem_long]

    min_map_cnt = max([0, int(math.ceil(rem_long/(rem_short-1)) - limit)])
    max_map_cnt = min([limit, rem_long])

    out = []
    for m in xrange(min_map_cnt, max_map_cnt+1):
        a = align(rem_short-1, rem_long-m, limit)
        if isinstance(a[0], list):
            out += [[m] + x for x in a]
        else:
            out += [[m] + [x] for x in a]
    return out



def gen_legal_alignments(seq1, seq2, max_sub_size = 3):
    if (len(seq2) - len(seq1)) / len(seq1) > max_sub_size-1:
        return None

    if len(seq1) == 1:                
        return [(seq1, [' '.join(seq2)])]

    ex =  align(len(seq1), len(seq2) - len(seq1), max_sub_size-1)            
    aligned_indices = [[x+1 for x in (ex[i])] for i in xrange(len(ex))]        
    # mapping_loc = defaultdict(lambda: defaultdict(set))
    alignments = []
    for al_ix_list in aligned_indices:
        idx = 0
        aligned_seq = []
        for i in xrange(len(seq1)):
            sub_seq2 = ' '.join(seq2[idx:idx+al_ix_list[i]])
            # mapping_loc[seq1[i]][sub_seq2].add(len(alignments))
            aligned_seq += [sub_seq2]
            idx += al_ix_list[i]
        alignments += [(seq1,aligned_seq)]
    
    return alignments

def convert_to_str(list_val):
    return '_'.join(list_val)
        
def em(paired_words, max_iter = 100):
    alignments = []    
    for p_w in paired_words:
        english, katakana = p_w[0], p_w[1]
        alignments += gen_legal_alignments(english, katakana)

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