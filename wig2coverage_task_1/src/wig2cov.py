import wiggelen
import os
import argparse
import gzip

# global counts
n0 = 0
n1 = 0
pp1 = 0
pp2 = 0
p1 = 0
p2 = 0
wp = 0
N0 = 0
N1 = 0
LT = 0

def parseOptions():
    epilog = """produce interval file annotated with counts of covered bases in this interval + previous non-interval region"""
    desc = "sums of covered bases in this interval ."
    parser = argparse.ArgumentParser(description=desc, epilog=epilog)

    parser.add_argument('-w','--wig_coverage_file', type=str, help='mutect 1  coverage wig file')
    parser.add_argument("-i","--interval_list_bed", type=str, help="interval_list_bed file")
    parser.add_argument("-o","--output_area", type=str, help="output_area",default=".")
    args = parser.parse_args()
    return args

def is_int(value):
  try:
    int(value)
    return True
  except:
    return False

def genome_pos(chr,pos):
    if is_int(chr):
        ichr=int(chr)
    elif chr == 'X':
        ichr = 23
    elif chr == 'Y':
        ichr=24
    elif chr == 'M':
        ichr=25
    else:
        ichr=100

    return (ichr-1)*1000000000+int(pos)

if __name__ == '__main__':

    args = parseOptions()

    wig_coverage_file = args.wig_coverage_file
    interval_list_bed = args.interval_list_bed
    output_area=args.output_area

    #f=open('/Users/stewart/Downloads/gaf_20111020+broad_wex_1.1_hg19.bed','rt')
    #of=open('/Users/stewart/Downloads/gaf_20111020+broad_wex_1.1_hg19.cov.bed','wt')
    base = os.path.basename(wig_coverage_file)

    if wig_coverage_file.endswith('.gz'):
        fw=gzip.open(wig_coverage_file)
        base, file_extension = os.path.splitext(base)
    else:
        fw=open(wig_coverage_file,'rt')

    if interval_list_bed.endswith('.gz'):
        fbed=gzip.open(interval_list_bed)
    else:
        fbed=open(interval_list_bed,'rt')

    base, file_extension = os.path.splitext(base)
    ofile=output_area+'/'+base+'.covered_bases.tsv'

    of=open(ofile,'wt')
    L=fbed.readline().strip().split()
    # assume order of bed and wig both in chromosome position order
    n0=0
    n1=0
    pp1=0
    pp2=0
    c0='A'
    p1=genome_pos(L[0],L[1])
    p2=genome_pos(L[0],L[2])
    n=0
    LT=LT+p2-p1

    for x in wiggelen.walk(fw):
        n=n+1
        wp=genome_pos(x[0],x[1])
        if wp > p2:   # write out interval line, read next interval (pp1,pp2,p1,p2), reset n1,n2
            of.write('%s\t%s\t%d\t%d\t%d\n' % (L[0],L[1],int(L[2])-int(L[1]),n0,n1))
            if (L[0] != c0):
                print('%s\t%d' % (L[0],n))
                c0=L[0]

            N0 = N0 + n0
            N1 = N1 + n1
            n0=0
            n1=x[2]
            pp1=p1
            pp2=p2
            L = fbed.readline().strip().split()
            if len(L)<3:
                break

            p1 = genome_pos(L[0], L[1])
            p2 = genome_pos(L[0], L[2])
            LT = LT + p2 - p1
            if p1<pp1:
                print('interval list out of order\n')
                print('suggest something like:')
                print("cat test1.bed | sed 's/^X/23/g' | sed 's/^Y/24/g' | sed 's/^M/25/g'| grep -v '^G' | sort -t$'\t'  -k1,1n -k2,1n | sed 's/^23/X/g' | sed 's/^24/Y/g'| sed 's/^25/M/g'")
                print('%s %s %s' % (L[0],L[1],L[2]))
                of.close()
                exit(2)

        elif wp<p1:
            n0=n0+x[2]
        elif (wp>=p1)&(wp<=p2):
            n1=n1+x[2]

    of.flush()
    os.fsync(of)
    of.close()

    file0=output_area+'/'+base+'.covered_bases_summary.txt'
    file1=output_area+'/'+base+'.total_targeted_bases.txt'
    file2=output_area+'/'+base+'.total_non_targeted_covered_bases.txt'
    file3=output_area+'/'+base+'.total_targeted_covered_bases.txt'

    of=open(file0,'wt')
    of.write('%d\n%d\n%d\n' % (LT,N0,N1))
    of.close()
    of=open(file1,'wt')
    of.write('%d\n' % LT)
    of.close()
    of=open(file2,'wt')
    of.write('%d\n' % N0)
    of.close()
    of=open(file3,'wt')
    of.write('%d\n' % N1)
    of.close()
