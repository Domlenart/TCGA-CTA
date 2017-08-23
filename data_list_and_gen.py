import os
from firebrowse import fbget
import pickle
from pandas import DataFrame
import pandas as pd
import lifelines
from lifelines.statistics import multivariate_logrank_test
from matplotlib import pyplot as plt
import numpy
from math import sqrt, ceil
import csv

CD = os.path.dirname(os.path.abspath(__file__))
SEQUESTERED_GENES = [i.rstrip() for i in open('overlapping-acrosome-CD')]
# LAML has no rnaseq data for solid tumor (no suprise)
# 'UCEC', 'UCS', 'UVM'are female only so no interested in them for now
LIST_OF_CANCERS = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'COADREAD', 'DLBC', 'ESCA', 'GBM', 'GBMLGG', 'HNSC',
                   'KICH', 'KIPAN', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG',
                   'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'STES', 'TGCT', 'THCA', 'THYM']


#print 'List of cancers that will be analyzed: ' + str(LIST_OF_CANCERS)
#print 'List of genes to be analyzed: ' + str(SEQUESTERED_GENES)
print 'Script directory: ' + str(CD)


def find_two_closest_integers(X):
    N = ceil(sqrt(X))
    M = 0
    while True:
        if (X % N) == 0:
            M = X/N
            break
        else:
            N+=1

    result = []
    result.append(N)
    result.append(M)
    return result


def get_clinical_data(cohort):
    os.chdir(CD)
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


def merge_clinical_with_mrna_new(cohort):
    os.chdir(str(CD + '/' + cohort))

    # Load clinical data and convert it to pandas dataframe with first row as column names
    with open(str(cohort + '_clinical.txt')) as clinical_file:
        clinical_header = clinical_file.readline().rstrip().split('\t')
        clinical_data = [line.rstrip().split('\t') for line in clinical_file]

    pandas_data = pd.DataFrame.from_records(clinical_data, columns=clinical_header)

    os.chdir(str(CD + '/' + cohort + '/Gene_data_RAW'))

    # Get z-score data for all genes for each patient and append it to clinical file
    for i in SEQUESTERED_GENES:
        with open(str(cohort + '_' + i.upper() + '_mRNASeq_RAW.txt')) as gene_file:
            gene_header = gene_file.readline().rstrip().split('\t')
            gene_data = [line.rstrip().split('\t') for line in gene_file]
            gene_data = pd.DataFrame.from_records(gene_data, columns=gene_header)
            gene_data = gene_data[['tcga_participant_barcode', 'z-score']]
            gene_data.rename(columns={'z-score': i}, inplace=True)
        pandas_data = pandas_data.merge(gene_data, left_on='tcga_participant_barcode', right_on='tcga_participant_barcode')

    # Merging "days to death" with "..to last followup" and removing any rows where it still shows up as NA or 0
    pandas_data['time'] = pandas_data.apply(lambda row: row["days_to_last_followup"] if row["days_to_death"] == 'NA' else row["days_to_death"], axis=1)
    pandas_data['death_status'] = pandas_data.apply(lambda row: 0 if row['vital_status'] == 'alive' else 1, axis=1)
    pandas_data = pandas_data[pandas_data.time != 'NA']
    pandas_data = pandas_data[pandas_data.time >= 1]


    # Saving output as human readable and pickled data
    os.chdir(str(CD + '/' + cohort))
    pandas_data.to_csv('merged_data.csv')
    pandas_data.to_pickle('merged_data')


def label_genes_with_expression(cohort):
    os.chdir(str(CD + '/' + cohort))
    data = pd.read_pickle('merged_data')
    # Drop column where everything is NaN
    data = data.replace('nan', numpy.NaN)
    data = data.dropna(axis=1, how='all')

    def apply_expression_description(row, column):
        if row[column] < 0:
            return 'Underexpressed'
        elif row[column] > 0.5:
            return 'Overexpressed'
        else:
            return 'Normal_expression'

    for i in SEQUESTERED_GENES:
        try:
            new_col_name = str(i) + '_expression'
            gene = []
            gene.append(i)
            data[i] = data[i].apply(pd.to_numeric)
            data[new_col_name] = data.apply(apply_expression_description, axis=1, args=gene)
            # print data[new_col_name].value_counts()

        except KeyError:
            pass

    data.to_csv('merged_and_labeled.csv')
    data.to_pickle('merged_and_labeled')


def kaplan_meier_plot_and_stat_new(cohort):

    # TODO: Add p value onto plots.
    # TODO: Add another version of logrank test to compare between each group instead of ANOVA

    os.chdir(str(CD + '/' + cohort))
    # Load pandas data
    data = pd.read_pickle('merged_and_labeled')

    def run_survival(data, gene_name):
        ay = plt.subplot(111)
        ay.set_title(gene_name)

        gene = gene_name + '_expression'
        gene = ''.join(gene)
        genders = ['male', 'female']
        group_by = ['gender']
        group_by.append(gene)
        # print group_by
        gene_groups = ['Underexpressed', 'Overexpressed', 'Normal_expression']

        kmf = lifelines.KaplanMeierFitter()

        grouped_data = data.groupby(group_by)

        for gene_group in gene_groups:
                for gender in genders:
                    try:
                        pre_tuple_list = [gender, gene_group]
                        group = tuple(pre_tuple_list)

                        # print 'tuple: ' + str(group)
                        d = grouped_data.get_group(group)
                        kaplan_meier_time = pd.to_numeric(d['time'])
                        kaplan_meier_event = d['death_status']

                        # TODO: Change label to display N
                        n_patients = [len(d)]
                        pre_tuple_list.append(n_patients)
                        label = str(pre_tuple_list)

                        kmf.fit(kaplan_meier_time, kaplan_meier_event, label=label)

                        kmf.plot(ax=ay, show_censors=True, ci_show=False)
                    except KeyError:
                        # print "No " + str(gender) + ' in gene' + str(gene_group)
                        pass

        event_durations = data.as_matrix(columns=['time'])

        data['stat_col'] = data[gene] + data['gender']
        group_labels = data.as_matrix(columns=['stat_col'])

        event = numpy.array(data.as_matrix(columns=['death_status']))

        result = multivariate_logrank_test(event_durations, group_labels, event, 0.85)

        os.chdir(str(CD + '/' + cohort))
        if not os.path.exists(CD + '/' + cohort + '/results'):
            os.makedirs('results')
        os.chdir(str(CD + '/' + cohort + '/results'))

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.savefig(str(gene_name + '_' + str(result.is_significant) + '.png'), bbox_inches='tight')

        # f.close()
        plt.close()

    genes = []
    for i in SEQUESTERED_GENES:
        if i in data:
            genes.append(i)

    for gene in genes:
        gene = str(gene)
        run_survival(data, gene)


def process_data_new(cohort):
    print "Processing data for: " + cohort
    get_clinical_data(cohort)
    print "Clinical data downloaded"
    get_mRNASeq_data(cohort)
    print "mRNASeq data downloaded"
    check_mRNASeq_data(cohort)
    print "Merging data"
    merge_clinical_with_mrna_new(cohort)
    print "Files merged. Labeling"
    label_genes_with_expression(cohort)
    print "Plotting."
    kaplan_meier_plot_and_stat_new(cohort)
    print "-" * 45


for i in LIST_OF_CANCERS:
    process_data_new(i)
