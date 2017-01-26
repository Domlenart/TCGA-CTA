import timeit
import fbget


with open('common-acrosome-cta') as f:
    list_of_genes = f.readlines()
list_of_genes = [x.strip() for x in list_of_genes]


def chunk(seq, size):

    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def get_genes_data(patient_barcode, cohort):

    for part in chunk(list_of_genes, 50):

        gene_data = fbget.mrnaseq(format="tsv", gene=part, cohort=cohort, barcode=patient_barcode, sample_type='TP')

print timeit.timeit(get_genes)
