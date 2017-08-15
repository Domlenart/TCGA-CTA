import os
from firebrowse import fbget

CD = os.path.dirname(os.path.abspath(__file__))
SEQUESTERED_GENES = [i.rstrip() for i in open('pure-sequestered-cta')]
LIST_OF_CANCERS = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'COADREAD', 'DLBC', 'ESCA', 'GBM', 'GBMLGG', 'HNSC',
                   'KICH', 'KIPAN', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG',
                   'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'STES', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

print 'List of cancers that will be analyzed: ' + str(LIST_OF_CANCERS)
print 'List of genes to be analyzed: ' + str(SEQUESTERED_GENES)
print 'Current script working directory: ' + str(CD)




def get_clinical_data(cohort):
    if not os.path.exists(cohort):
        os.makedirs(cohort)
        working_directory = CD + '/' + cohort
        os.chdir(working_directory)
        filename = str(cohort) + "_clinical.txt"
        f = open(filename, "w")
        clinical_data = str(fbget.clinical(format='tsv', cohort=cohort, sort_by='tcga_participant_barcode'))
        f.write(clinical_data)
        f.close()
    else:
        print "Clinical data for " + str(cohort) + " exists. Skipping download."


def get_mRNASeq_data(cohort):
    working_directory = CD + '/' + cohort
    os.chdir(working_directory)
    working_directory = working_directory + "/Gene_data_RAW"
    if not os.path.exists(working_directory):
        try:
            os.makedirs("Gene_data_RAW")
        except OSError:
            pass
        os.chdir(working_directory)
        for gene_name in SEQUESTERED_GENES:
            filename = str(cohort) + '_' + str(gene_name).upper() + '_mRNASeq_RAW.txt'
            f = open(filename, 'w')
            mRNASeq_data = str(fbget.mrnaseq(format='tsv', gene=gene_name, cohort=cohort, sample_type='TP',
                                             sort_by='tcga_participant_barcode'))
            f.write(mRNASeq_data)
            f.close()
    else:
        print "mRNA_RAW data for " + str(cohort) + " exists. Skipping download."


def check_mRNASeq_data(cohort):

    ''' Remove nonexistent clinical patients from mRNASeq files'''

    working_dir = CD + '/' + cohort
    os.chdir(working_dir)

    try:
        filename = str(cohort) + "_clinical.txt"
        clinical_data = open(filename)
    except OSError:
        print "Can't open clinical file - TRIM STAGE"

    header = clinical_data.readline().rstrip().split('\t')

    clinical_data = [line.rstrip().split('\t') for line in clinical_data]
    clinical_barcodes = [sample[0] for sample in clinical_data]

    working_dir = CD + '/' + cohort + '/Gene_data_RAW'
    os.chdir(working_dir)

    for gene in SEQUESTERED_GENES:

        try:
            filename = str(cohort) + '_' + str(gene).upper() + '_mRNASeq_RAW.txt'
            gene_data = open(filename)
        except OSError:
            print "Can't open gene data for " + str(gene) + ' - TRIM STAGE'

        header = gene_data.readline().rstrip().split('\t')

        gene_data = [line.rstrip().split('\t') for line in gene_data]
        gene_barcodes = [sample[0] for sample in clinical_data]

        if set(gene_barcodes) == set(clinical_barcodes):
            pass
        else:
            print 'Barcodes DO NOT overlap for gene ' + str(gene)
            break
    print 'No problems with mRNASeq in ' + str(cohort)

def process_data(cohort):
    get_clinical_data(cohort)
    get_mRNASeq_data(cohort)

process_data("BRCA")

key = 'x'
dict = {key: ['a', 'b']}
print dict

check_mRNASeq_data('BRCA')