from Bio import SeqIO
from pathlib import Path

if __name__ == "__main__":
    all_genes = []

    for ffn in Path('genes').glob('*.ffn'):
        for record in SeqIO.parse(ffn, 'fasta'):
            all_genes.append(record)

    covr_accs = []
    covs_accs = []
    for vf_hits in Path("vf_hits").glob('*.out6'):
        with open(vf_hits) as f:
            for line in f:
                line = line.split('\t')

                if line[1] == 'covR':
                    covr_accs.append(line[0])

                if line[1] == 'covS':
                    covs_accs.append(line[0])

    covr_alleles = []
    covs_alleles = []
    for record in all_genes:
        if record.id in covr_accs:
            covr_alleles.append(record)
        elif record.id in covs_accs:
            covs_alleles.append(record)

    SeqIO.write(covr_alleles, 'covR_covS_alignment/covR_alleles.ffn', 'fasta')
    SeqIO.write(covs_alleles, 'covR_covS_alignment/covS_alleles.ffn', 'fasta')
