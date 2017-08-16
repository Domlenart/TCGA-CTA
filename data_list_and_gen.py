import os
from firebrowse import fbget
import pickle

from pandas import DataFrame
import pandas as pd
import lifelines
from matplotlib import pyplot as plt

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

def merge_clinical_with_mrna(cohort):

    ''' Order of data:
    1. Barcode
    2. Death status (1 = dead, 0 = alive)
    3. Days to death (NA if alive)
    4. Days to last followup (NA if dead) '''

    os.chdir(str(CD + '/' + cohort))
    if not os.path.isfile(str(CD + '/' + cohort + '/' + 'merged_data_' + cohort)):

        print 'Merged file for ' + str(cohort) + 'does not exist, processing...'

        clinical_data = open(str(cohort) + '_clinical.txt')
        clinical_header = clinical_data.readline().rstrip().split('\t')
        clinical_header_barcode_pos = clinical_header.index('tcga_participant_barcode')
        clinical_vital_pos = clinical_header.index('vital_status')
        clinical_to_death_pos = clinical_header.index('days_to_death')
        clinical_last_followup = clinical_header.index('days_to_last_followup')

        clinical_data = [line.rstrip().split('\t') for line in clinical_data]
        relevant_data = []
        clinical_relevant_data_header = ['tcga_participant_barcode', 'death_status', 'days_to_death',
                                         'days_to_last_followup', 'time_alive']
        relevant_data.append(clinical_relevant_data_header)

        # Extracting clinical data
        for sample in clinical_data:
            barcode, vital, to_death, last_followup = sample[clinical_header_barcode_pos], sample[clinical_vital_pos], \
                                                      sample[clinical_to_death_pos], sample[clinical_last_followup]
            if vital == 'alive':
                time_alive = last_followup
                vital = 0

            else:
                time_alive = to_death
                vital = 1

            lst = [barcode, vital, to_death, last_followup, time_alive]

            relevant_data.append(lst)

        os.chdir(str(CD + '/' + cohort + '/Gene_data_RAW'))

        # Extracting gene data
        for gene in SEQUESTERED_GENES:
            relevant_data[0].append(gene)

            try:
                filename = str(cohort) + '_' + str(gene).upper() + '_mRNASeq_RAW.txt'
                gene_data = open(filename)
            except OSError:
                print "Can't open gene data for " + str(gene) + ' - MERGE STAGE'

            gene_data_header = gene_data.readline().rstrip().split('\t')
            gene_data_z_score_pos = gene_data_header.index('z-score')

            gene_data = [line.rstrip().split('\t') for line in gene_data]
            id = 0

            for line in gene_data:
                id += 1
                if line[0] == relevant_data[id][0]:
                    z_score = line[gene_data_z_score_pos]
                    relevant_data[id].append(z_score)

                else:
                    relevant_data[id].append('Not found')
                    z_score = line[gene_data_z_score_pos]
                    id += 1
                    relevant_data[id].append(z_score)

        os.chdir(str(CD + '/' + cohort))

        # Removing patients that have wrong survival data from file
        for i in relevant_data:
            if i == relevant_data[0]:
                pass
            else:
                if i[4] == 'NA' or int(i[4]) < 0:
                    relevant_data.pop(relevant_data.index(i))
                    print "Removed entry"

        # Dumping a list with relevant data to file
        with open(str('merged_data_' + cohort), 'wb') as fp:
            pickle.dump(relevant_data, fp)

        return 0
    else:
        with open(str('merged_data_' + cohort), 'rb') as fp:
            DATA = pickle.load(fp)

        return DATA


def kaplan_meier_plot_and_stat(data):

    for i in data:
        print i

    header = data[0]
    data = data[1:]
    pandas_data = DataFrame(data=data, columns=header)

    # Need to convert pd data series to numeric, because it is incorrectly recognized as object.
    kaplan_meier_time = pd.to_numeric(pandas_data['time_alive'])
    kaplan_meier_event = pandas_data['death_status']
    fitter = lifelines.KaplanMeierFitter()
    fitter.fit(kaplan_meier_time, event_observed=kaplan_meier_event)
    fitter.plot()
    plt.show()


def process_data(cohort):
    get_clinical_data(cohort)
    get_mRNASeq_data(cohort)
    check_mRNASeq_data(cohort)



x = merge_clinical_with_mrna('BRCA')

kaplan_meier_plot_and_stat(x)