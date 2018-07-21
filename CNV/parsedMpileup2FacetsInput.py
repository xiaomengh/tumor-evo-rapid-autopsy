from sys import stdin, argv

normalSampleColumn=int(argv[1])
samples = []
for line in stdin:
    line=line.strip()
    a = line.split('\t')
    if a[0]=='#chrom':
        for i in range(int((len(a)-4)/3)):
           parsed_name = a[3*i+4].split('_')
           if len(parsed_name)==2:
               s=parsed_name[0]
           else:
               s='_'.join(parsed_name[0:(len(parsed_name)-1)])
           samples.append(s)
        outFiles = [open(t+'_Germ1.facets.sort.txt', mode='wt') for t in samples]
    else:
        for j in range(len(outFiles)):
            outFiles[j].write('\n'+'\t'.join(a[0:2])+'\t'+'\t'.join(a[normalSampleColumn:normalSampleColumn+2])+'\t'+'\t'.join(a[3*j+4:3*j+6]))

[f.close() for f in outFiles]
