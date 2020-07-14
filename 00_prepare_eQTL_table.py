import pandas as pd
import numpy as np

import re
def remove_chr(ss):
    return [ re.sub('chr', '', s) for s in ss ]

gene = 'ENSG00000008128'

df = pd.read_parquet('~/scratch/mixQTL-GTExV8/mixqtl/Liver/mixqtl.Liver_GTEx_eGene.cis_qtl_pairs.mixQTL.chr1.parquet')
kk = df[ df.phenotype_id == 'ENSG00000008128' ]
# print(kk.columns)
kk = kk[ ['variant_id', 'tstat_trc', 'pval_trc' ]]
kk['variant_id'] = remove_chr(kk['variant_id'].values)

# kk.to_csv('eQTL_score.tsv', sep='\t', index=False, header=False)

h1 = np.loadtxt('/home/t.cri.yliang/test_plasma/hap1.txt.gz')
h2 = np.loadtxt('/home/t.cri.yliang/test_plasma/hap2.txt.gz')
geno = h1 + h2
geno = geno.astype(int)
row_sd = np.std(geno, axis=1)
# geno = geno[row_sd != 0, :]
geno_mm = h1 - h2
geno_mm = geno_mm.astype(int)
row_sd_mm = np.std(geno_mm, axis=1)
good_ind = np.logical_and(row_sd != 0, row_sd_mm != 0)
geno = geno[ good_ind, :]

kk = kk[ good_ind ]
kk.to_csv('eQTL_score.tsv', sep='\t', index=False, header=False)

ld = np.corrcoef(geno)
# print(geno[:3, :3])
# geno = np.concatenate((kk.variant_id.values[:, np.newaxis], geno), axis=1)
np.savetxt('tmp_genotype.tsv', geno, fmt="%d", delimiter='\t')
np.savetxt('ld.tsv', ld, delimiter='\t')
kk[['variant_id']].to_csv('temp', index=False, header=False)
header = [ 'Id' ] + [ 'indiv_' + str(i) for i in range(h1.shape[1]) ]
with open('tt2', 'w') as f:
    f.write('\t'.join(header) + '\n')
import os
os.system('paste temp tmp_genotype.tsv > ttgenotype.tsv')
os.system('cat tt2 ttgenotype.tsv > genotype.tsv')
os.system('rm tt2')
os.system('rm ttgenotype.tsv')
os.system('rm tmp_genotype.tsv')
os.system('rm temp')

a1 = np.loadtxt('/home/t.cri.yliang/test_plasma/ref.txt.gz')
a2 = np.loadtxt('/home/t.cri.yliang/test_plasma/alt.txt.gz')
a1 = a1.astype(int)
a2 = a2.astype(int)
# a1[a1 == 0] = 0
# a2[a2 == 0] = 0
df = pd.DataFrame({'indiv': [ 'indiv_' + str(i) for i in range(h1.shape[1]) ], 'a1': a1, 'a2': a2})
df.columns = ['SAMPLE_ID', 'H1_COUNT', 'H2_COUNT']
df.to_csv('ASE.tsv', sep='\t', index=False, header=True)



