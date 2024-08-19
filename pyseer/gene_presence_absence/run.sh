pyseer --phenotypes ../invasive.pheno --phenotype-column "Invasive" \
    --pres ../../pangenome/bactopia-runs/SHL_igas-20240616-024658/panaroo/gene_presence_absence.Rtab \
    --similarity ../phylogeny_K.tsv --lmm --uncompressed \
    --output-patterns gene_presence_absence_patterns.txt --cpu 8 \
    > gene_presence_absence_results.txt 2>> gene_presence_absence_results.log

sleep 5 

cat <(head -1 gene_presence_absence_results.txt) <(LC_ALL=C awk -v pval=$(python ../scripts/count_patterns_microgwas.py --threshold gene_presence_absence_patterns.txt) '$4<pval {print $0}' gene_presence_absence_results.txt) > gene_presence_absence_patterns_filtered.txt

python ../scripts/gpa_summary.py


pyseer --phenotypes ../invasive.pheno --phenotype-column "Invasive" \
    --pres ../../pangenome/bactopia-runs/SHL_igas-20240616-024658/panaroo/struct_presence_absence.Rtab \
    --similarity ../phylogeny_K.tsv --lmm --uncompressed \
    --output-patterns struct_presence_absence_patterns.txt --cpu 8 \
    > struct_presence_absence_results.txt 2>> struct_presence_absence_results.log

sleep 5 

cat <(head -1 struct_presence_absence_results.txt) <(LC_ALL=C awk -v pval=$(python ../scripts/count_patterns_microgwas.py --threshold struct_presence_absence_patterns.txt) '$4<pval {print $0}' struct_presence_absence_results.txt) > struct_presence_absence_patterns_filtered.txt
