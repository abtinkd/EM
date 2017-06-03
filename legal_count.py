from __future__ import division
from math import factorial as fact
from em import gen_legal_alignments
import sys

def comb(n,k):    
    return fact(n)/(fact(n-k)*fact(k))

def legal_alignment(n,m,k=3):
    if m/n > k:
        return 0
    
    count1 = len(gen_legal_alignments('a'*n, 'b'*m, k))

    d = (m-n)//k    
    r = (m-n)%k    
    if m-n>k-1:
        print 'a' #comb(m-1, m-n) = fact(m-1)/(fact(n-1)*fact(m-n))
        count2 = comb(m-1, m-n)  - n*comb(m-k-1, m-n-k) + (((m-n-k)//k)>0)*comb(n,d)*comb(n-1+r, r)
    else:
        print 'b'
        count2  = comb(m-1, m-n)

    return count1, count2#, (fact(m-1)/(fact(n-1)*fact(m-n)) - count1),(fact(m-k-1)/(fact(n-1)*fact(m-n-k)))


if __name__=='__main__':
    n,m = int(sys.argv[1]), int(sys.argv[2])
    if len(sys.argv) == 4:        
        print legal_alignment(n,m,int(sys.argv[3]))
    else:
        print legal_alignment(n,m)

