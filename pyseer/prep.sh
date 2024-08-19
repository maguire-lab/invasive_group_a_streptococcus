python scripts/phylogeny_distance.py --lmm ../../pangenome/bactopia-runs/SHL_igas-20240616-024658/iqtree/core-genome.treefile > phylogeny_K.tsv
python scripts/phylogeny_distance.py ../../pangenome/bactopia-runs/SHL_igas-20240616-024658/iqtree/core-genome.treefile > phylogeny_distances.tsv
cat ../poppunk/output/output_clusters.csv | sed 's/,/\t/' | sed '1d' > lineages.txt

ls -d -1 $PWD/../bactopia/*/main/qc/*.fastq.gz > reads.txt
ls -d -1 $PWD/../bactopia/*/main/annotator/prokka/*.fna.gz > refs.txt
unitig-caller --call --refs refs.txt --reads reads.txt --pyseer --threads 8 --write-graph
# fix names in unitig data
cat untig_caller_genomes.pyseer | sed 's/\.fna//g' | gzip > untig_caller_genomes.fixed_names.pyseer.gz



mkdir -p gffs
for gff in $(ls -d -1 $PWD/../bactopia/*/main/annotator/prokka/*.gff.gz);
do 
    zcat $gff > gffs/$(basename $gff | sed 's/\.gz//');
done
ls -d -1 $PWD/gffs/*.gff > gffs.txt
paste -d '\t' refs.txt gffs.txt > references.txt


