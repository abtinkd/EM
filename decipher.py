from __future__ import division
from collections import defaultdict
from convert_wfsa_prob import convert_wfsa_prob
import sys

ALPHABET = '_abcdefghijklmnopqrstuvwxyz'

def build_prob_table():    
    prob_table = defaultdict(lambda: defaultdict(float))
    for alp1 in ALPHABET:
        for alp2 in ALPHABET:
            prob_table[alp1][alp2] = 1/len(ALPHABET)
    prob_table['</s>']['</s>'] = 1.0    
    return prob_table

def build_icecream_prob_table():
    prob_table = defaultdict(lambda: defaultdict(float))
    prob_table['1']['C'] = 0.7
    prob_table['1']['H'] = 0.1
    prob_table['2']['C'] = 0.2
    prob_table['2']['H'] = 0.2
    prob_table['3']['C'] = 0.1
    prob_table['3']['H'] = 0.7
    prob_table['</s>']['</s>'] = 1.0
    return prob_table



def decipher(ciphered_list, LM, iteration = 100):
    # prob_table = build_prob_table()
    prob_table = build_icecream_prob_table()
    
    cipher = prob_table.keys()
    cipher.remove('</s>')
    states = prob_table[cipher[0]].keys()

    m_prev_post = LM[0]
    m_post_prev = LM[1]    

    for i in xrange(iteration):        
        for ciph_line in ciphered_list:
            
            alpha, beta = defaultdict(lambda: defaultdict(float)), defaultdict(lambda: defaultdict(float))                
            start_index = 0
            
            if '<s><s>' in m_prev_post.keys():
                ciph_line = ['<s><s>'] + [c for c in ciph_line.strip()] + ['</s></s>']
                start_index = 1
            elif '<s>' in m_prev_post:                
                ciph_line = ['<s>'] + [c for c in ciph_line.strip()] + ['</s>'] 
                start_index = 1
                alpha[0]['<s>'] = 1.0
                beta[len(ciph_line)-1]['</s>'] = 1.0
                                                                
            #calculate alpha & beta
            end_index = len(ciph_line)-start_index
            for index in xrange(start_index, end_index):                
                for prev in alpha[index - start_index]:
                    for post in m_prev_post[prev]:                                                    
                        alpha[index][post] += alpha[index - start_index][prev]* \
                                        m_prev_post[prev][post]*prob_table[ciph_line[index]][post]                                                                        
                for post in beta[end_index-index+1]:
                    for prev in m_post_prev[post]:                                                                           
                        beta[end_index-index][prev] += beta[end_index-index+1][post]* \
                                        m_post_prev[post][prev]*prob_table[ciph_line[end_index-index+1]][post]

            
            #calculate alpha*beta
            prod_ab = defaultdict(lambda: defaultdict(float))
            prob_state = defaultdict(lambda: defaultdict(float))
            prob_input_state = defaultdict(lambda: defaultdict(float))
            tot_prob = defaultdict(float)                    
            for index in xrange(start_index, end_index):
                for st in states:  
                    prod_ab[index][st] = alpha[index][st] * beta[index][st]                                                            
                    tot_prob[index] += prod_ab[index][st]

                for st in states:
                    prob_state[index][st] = prod_ab[index][st] / tot_prob[index]
                    prob_input_state[ciph_line[index]][st] += prod_ab[index][st] / tot_prob[index]

            for inp in prob_table:
                for st in states:                    
                    prob_table[inp][st] = prob_input_state[inp][st]/ sum([prob_state[i][st] for i in prob_state])

            
            for i in prob_table:
                for j in prob_table[i]:
                    print 'P({0}|{1}) = {2}'.format(i,j,prob_table[i][j])
            raw_input()
        




def build_model(filename):
    model_0_1 = defaultdict(lambda: defaultdict(float))
    model_1_0 = defaultdict(lambda: defaultdict(float))
    with open(filename, 'r') as fp:
        for line in fp:
            line  = line.replace(':','').replace('#','').strip().split()
            if len(line) != 3:
                return None
            model_0_1[line[0]][line[1]] = float(line[2])
            model_1_0[line[1]][line[0]] = float(line[2])

    return model_0_1, model_1_0

def build_icecream_model():
    model_0_1 = defaultdict(lambda: defaultdict(float))
    model_1_0 = defaultdict(lambda: defaultdict(float))

    model_0_1['C']['C'] = 0.8
    model_0_1['C']['H'] = 0.1
    model_0_1['C']['</s>'] = 0.1
    model_0_1['H']['C'] = 0.1
    model_0_1['H']['H'] = 0.8
    model_0_1['H']['</s>'] = 0.1
    model_0_1['<s>']['C'] = 0.5
    model_0_1['<s>']['H'] = 0.5
    model_0_1['<s>']['</s>'] = 0.0

    for i in model_0_1:
        for j in model_0_1[i]:
            model_1_0[j][i] = model_0_1[i][j]    

    return (model_0_1, model_1_0)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        ciphered_file = sys.argv[1]
        if len(sys.argv) > 2:
            model_prob_file = sys.argv[2]
        else:
            model_prob_file = 'bigram.wfsa.prob'
            max_iterations = 100
            if len(sys.argv) > 3:
                max_iterations = sys.argv[3]
            else:
                max_iterations = 100
    else:
        ciphered_file = 'cipher.txt'
        model_prob_file = 'bigram.wfsa.prob'
        max_iterations = 100

    with open(ciphered_file, 'r') as fp:
        ciphered_list = fp.readlines()

    LM = build_model(model_prob_file)

    LM = build_icecream_model()
    ciphered_list = [''.join(map(str,[2,3,3,2,3,2,3,2,2,3,1,3,3,1,1,1,2,1,1,1,3,1,2,1,1,1,2,3,3,2,3,2,2]))]

    print decipher(ciphered_list, LM, max_iterations)