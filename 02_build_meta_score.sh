out1=ase_het.txt
ase=test_mapase
python ~/labshare/softwares/aim/sample_pipeline/computeLDMatrix.py $ase genotype.tsv $out1

out2=ase_het_psd.txt
Rscript ~/labshare/softwares/aim/sample_pipeline/nearPD.R $out1 $out2


eqtl=eQTL_score.tsv
score1=meta_score1.txt
score2=meta_score2.txt
Rscript ~/labshare/softwares/aim/sample_pipeline/computeMetaScores.R $ase $eqtl $score1 $score2

ld=ld.tsv
metald=meta_ld.tsv
Rscript ~/labshare/softwares/aim/sample_pipeline/computeMetaLD.R $ld $ase_het.txt $metald


