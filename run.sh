#!/bin/zsh
set -euo pipefail
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#bactopia --samples bactopia_sample_sheet.csv --outdir bactopia --run_name SHL_igas -profile docker 
#bactopia --wf emmtyper --bactopia bactopia --oudtir emmtyper --run_name SHL_igas -profile docker 
#bactopia --wf pangenome --bactopia bactopia --outdir pangenome --run_name SHL_igas -profile docker --use_panaroo

#conda activate virulence
## Then perform assembly typing to identify M1UK SNP sets
#cd assembly_snptyper/assembly_snptyper
#ls -d -1 ../../bactopia/*/main/annotator/prokka/*.fna.gz > fastas.txt
#assembly_snptyper --vcf data/M1UK.vcf --reference data/MGAS5005.fa --list_input fastas.txt -p 4 --verbose > output.txt
#cp output.txt ..
#cd ../../

## Analysis virulence and adehsion factors
#cd virulence_adhesion_factors
#mkdir -p vf_hits genes
#for ffn in ../bactopia/*/main/annotator/prokka/*.ffn.gz; do zcat $ffn > genes/$(basename $ffn | sed 's/.gz//'); done
#
#makeblastdb -in virulence_factors.fas -dbtype nucl
#for gene_set in genes/*.ffn; 
#    do 
#        blastn -query $gene_set -db virulence_factors.fas -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus' > vf_hits/$(basename $gene_set | sed 's/.ffn//').out6
#        echo $gene_set
#done
## CovR and CovS alignments
#python get_covR_covS.py
#cd covR_covS_alignment
#cat covR_reference_seq.fas covR_alleles.ffn | minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 -t 4 covR_reference_seq.fas - -o covR_mapped.sam
#gofasta sam toMultiAlign -s covR_mapped.sam -t 4 --reference covR_reference_seq.fas > covR_aligned.fas
#cat covS_reference_seq.fas covS_alleles.ffn | minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 -t 4 covS_reference_seq.fas - -o covS_mapped.sam
#gofasta sam toMultiAlign -s covS_mapped.sam -t 4 --reference covS_reference_seq.fas > covS_aligned.fas
#cd ../..

## Run poppunk
#cd poppunk
#mkdir -p genomes
#for genome in ../bactopia/*/main/annotator/prokka/*.fna.gz; 
#do 
#    zcat $genome >> genomes/$(basename $genome | sed 's/.gz//')
#done
#cat /dev/null > genomes.txt
#for i in $(find genomes -name '*.fna' ); do 
#    echo -e "$(basename $i | sed 's/\.fna//')\t$i" >> genomes.txt
#done 

#poppunk_assign --db gas_database --query genomes.txt --output output --threads 10
#poppunk_visualise --ref-db gas_database --query-db output --tree mst --output poppunk_mst --microreact 
#cd ..

conda activate pyseer
cd pyseer
#ls -d -1 $PWD/../bactopia/*/main/annotator/prokka/*.fna.gz > refs.txt
#
#mkdir -p gffs
#for gff in $(ls -d -1 $PWD/../bactopia/*/main/annotator/prokka/*.gff.gz);
#do 
#    zcat $gff > gffs/$(basename $gff | sed 's/\.gz//');
#done
ls -d -1 $PWD/gffs/*.gff > gffs.txt
paste -d '\t' refs.txt gffs.txt <(yes "ref" | head -n $(cat refs.txt | wc -l)) > ref_and_gff.txt
rm gffs.txt 

#unitig-caller --call --refs refs.txt --pyseer --threads 8 --write-graph
#python phylogeny_distance.py --lmm ../pangenome/bactopia-runs/SHL_igas-20240616-024658/iqtree/core-genome.treefile > phylogeny_K.tsv
#python phylogeny_distance.py ../pangenome/bactopia-runs/SHL_igas-20240616-024658/iqtree/core-genome.treefile > phylogeny_distance.tsv
#
## fix names in unitig data
#cat unitig_caller.pyseer | sed 's/\.fna//g' | gzip > untig_caller_genomes.fixed_names.pyseer.gz
#rm unitig_caller.pyseer

#cut -f1,3 ../gas_metadata.csv | sed 's/GAS Isolate/samples/' | sed 's/Non-Invasive/0/' | sed 's/Invasive/1/' | sed 's/Location/Invasive/' > invasive.pheno
#cat ../poppunk/output/output_clusters.csv | sed 's/,/\t/' | sed '1d' > lineages.txt
#pyseer --kmers untig_caller_genomes.fixed_names.pyseer.gz --lmm --similarity phylogeny_K.tsv --distance phylogeny_distance.tsv --phenotypes invasive.pheno --lineage --lineage-file lineage_effects.txt --lineage-clusters lineages.txt > pyseer_genome_unitigs_lineages_results.txt
echo "P-value threshold used: $(python scripts/count_patterns.py unitig_patterns.txt | grep "Threshold"  | cut -f2)"
cat <(head -1 pyseer_genome_unitigs_lineages_results.txt) <(awk '$4<6.16E-07 {print $0}' pyseer_genome_unitigs_lineages_results.txt) > significant_lineage_kmers.txt


annotate_hits_pyseer significant_lineage_kmers.txt ref_and_gff.txt significant_lineage_kmers_annotated.txt
python scripts/summarise_annotations.py --nearby --unadj-p significant_lineage_kmers_annotated.txt > significant_lineage_kmers_annotated_summary.txt
cd ..

#annotate_hits_pyseer significant_kmers.txt references.txt annotated_kmers.txt
#mkdir -p gffs; 
#for i in /workspace/fin/SHL/bacteria/analyses/invasive_group_a_streptococcus_2024_clean/pyseer/../20240122_gas/*/main/annotator/prokka/*.gff.gz; do zcat $i > gffs/$(basename $i | sed 's/\.gz//'); done
#paste  <(cat refs.txt)  <(cut -d '/' -f 15 refs.txt | sed 's/\.fna\.gz/\.gff/' | sed 's/^/gffs\//') | sed 's/$/\tdraft/' > ref_genomes_gffs_annot.txt
