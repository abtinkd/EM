import sys

def convert_wfsa_prob(filename):
    with open(filename, 'r') as fpr:
        with open(filename+'.prob', 'w') as fpw:
            for line in fpr:
                line = line.replace('(','').replace(')', '').strip().split()
                if len(line) < 3:
                    continue
                if len(line) == 3:
                    fpw.write('{0} : {1} # {2}\n'.format(line[0], line[2], 1.0))
                else:
                    fpw.write('{0} : {1} # {2}\n'.format(line[0], line[2], line[3]))

    

if __name__=='__main__':
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename ='bigram.wfsa'

    convert_wfsa_prob(filename)

