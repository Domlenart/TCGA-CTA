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
SEQUESTERED_GENES = [i.rstrip() for i in open('pure-sequestered-cta')]
LIST_OF_CANCERS = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'COADREAD', 'DLBC', 'ESCA', 'GBM', 'GBMLGG', 'HNSC',
                   'KICH', 'KIPAN', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG',
                   'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'STES', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

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


def merge_clinical_with_mrna(cohort):

    ''' Order of data:
    0. Barcode
    1. Gender
    2. Death status (1 = dead, 0 = alive) - to use in censoring
    3. Days to death (NA if alive)
    4. Days to last followup (NA if dead)
    5. Time to death OR last followup
    6 - end. Data for each gene
    '''

    os.chdir(str(CD + '/' + cohort))
    if not os.path.isfile(str(CD + '/' + cohort + '/' + 'merged_data_' + cohort)):

        print 'Merged file for ' + str(cohort) + 'does not exist, processing...'

        clinical_data = open(str(cohort) + '_clinical.txt')
        clinical_header = clinical_data.readline().rstrip().split('\t')
        clinical_header_barcode_pos = clinical_header.index('tcga_participant_barcode')
        clinical_vital_pos = clinical_header.index('vital_status')
        clinical_to_death_pos = clinical_header.index('days_to_death')
        clinical_last_followup = clinical_header.index('days_to_last_followup')
        clinical_gender = clinical_header.index('gender')

        clinical_data = [line.rstrip().split('\t') for line in clinical_data]
        relevant_data = []
        clinical_relevant_data_header = ['tcga_participant_barcode', 'gender', 'death_status', 'days_to_death',
                                         'days_to_last_followup', 'time_alive']
        relevant_data.append(clinical_relevant_data_header)

        # Extracting clinical data
        for sample in clinical_data:
            barcode, gender,  vital, to_death, last_followup = sample[clinical_header_barcode_pos], sample[clinical_gender], \
                                                               sample[clinical_vital_pos], sample[clinical_to_death_pos], \
                                                               sample[clinical_last_followup]
            if vital == 'alive':
                time_alive = last_followup
                vital = 0

            else:
                time_alive = to_death
                vital = 1

            lst = [barcode, gender, vital, to_death, last_followup, time_alive]

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
                # If Barcodes match then append Z score
                try:
                    if line[0] == relevant_data[id][0]:
                        z_score = line[gene_data_z_score_pos]
                        relevant_data[id].append(z_score)

                    else:
                        relevant_data[id].append('Not found')
                        z_score = line[gene_data_z_score_pos]
                        id += 1
                        relevant_data[id].append(z_score)
                except IndexError:
                    pass
        os.chdir(str(CD + '/' + cohort))

        # Removing patients that have wrong survival data or no rnaseq from file
        for line in relevant_data:

            # Skip 1st line
            if line == relevant_data[0]:
                pass
            else:
                # Remove entries that have wrong survival time data
                if line[5] == 'NA' or int(line[5]) < 0 or int(line[5]) == 0:
                    relevant_data.pop(relevant_data.index(line))
                    print "Removed entry due to wrong lifespan data"

                else:
                    # For each column with gene data check if they are valid, if not - remove.
                    try:
                        for i in range(len(relevant_data[0][6:])):
                            if line[i] == 'Not found':
                                relevant_data.pop(relevant_data.index(line))
                                print "Removed entry due to wrong RNASeq data"
                                break
                            else:
                                pass
                    except IndexError:
                        print "Weird indexing problem occured, removing sample"
                        relevant_data.pop(relevant_data.index(line))

        # Temporary solution for LUAD data bug:

        if cohort == "LUAD":
            print "running fix for LUAD DATA"
            relevant_data.pop(1)

        # Dumping a list with relevant data to file
        with open('merged_data', 'wb') as ff:
            writer = csv.writer(ff)
            writer.writerows(relevant_data)
        with open(str('merged_data_' + cohort), 'wb') as fp:
            pickle.dump(relevant_data, fp)

    else:
        pass


def merge_clinical_with_mrna_new(cohort):
    os.chdir(str(CD + '/' + cohort))

    with open(str(cohort + '_clinical.txt')) as clinical_file:
        clinical_header = clinical_file.readline().rstrip().split('\t')
        clinical_data = [line.rstrip().split('\t') for line in clinical_file]

    pandas_data = pd.DataFrame.from_records(clinical_data, columns=clinical_header)

    os.chdir(str(CD + '/' + cohort + '/Gene_data_RAW'))

    for i in SEQUESTERED_GENES:
        with open(str(cohort + '_' + i.upper() + '_mRNASeq_RAW.txt')) as gene_file:
            gene_header = gene_file.readline().rstrip().split('\t')
            gene_data = [line.rstrip().split('\t') for line in gene_file]
            gene_data = pd.DataFrame.from_records(gene_data, columns=gene_header)
            gene_data = gene_data[['tcga_participant_barcode', 'z-score']]
            gene_data.rename(columns={'z-score': i}, inplace=True)

        pandas_data = pandas_data.merge(gene_data, left_on='tcga_participant_barcode', right_on='tcga_participant_barcode')

    # Merging "days to death" with "..to last followup" and removing any rows where it still shows up as NA or 0
    pandas_data["time"] = pandas_data.apply(lambda row: row["days_to_last_followup"] if row["days_to_death"] == 'NA' else row["days_to_death"], axis=1)

    pandas_data = pandas_data[pandas_data.time != 'NA']
    pandas_data = pandas_data[pandas_data.time > 0]

    os.chdir(str(CD + '/' + cohort))
    pandas_data.to_csv('test')


def label_genes_with_quartiles(cohort):
    os.chdir(str(CD + '/' + cohort))

    with open(str('merged_data_' + cohort), 'rb') as fp:
        data = pickle.load(fp)

    # Make labels for each quartile that will be added as last columns at the rightmost side
    clinical_labels = [i for i in data[0]]
    gene_labels = [i + ' Quartile' for i in data[0][6:]]

    # Strip the line with labels
    data = data[1:]
    # Determine quartiles for each gene expression. Apply quartiles to main data structure

    genes_with_varying_expression = []
    genes_with_no_changes_in_expression = []

    NEW_DATA_LIST = [clinical_labels + gene_labels] + data
    to_append_each_line = ['' for i in range(len(SEQUESTERED_GENES))]
    for i in NEW_DATA_LIST:
        i += to_append_each_line

    for gene in range(len(SEQUESTERED_GENES)):
        list_with_gene_expression_score = []
        for line in data:
            list_with_gene_expression_score.append(line[gene+6])

        list_with_gene_expression_score = [float(i) for i in list_with_gene_expression_score]
        dict_key = SEQUESTERED_GENES[gene]

        numpy_array = numpy.array(list_with_gene_expression_score)
        percentile_50 = numpy.percentile(numpy_array, 50)
        # print str("Percentile 50 for " + str(dict_key) + "is equal to: " + str(percentile_50))
        percentile_75 = numpy.percentile(numpy_array, 75)
        # print str("Percentile 75 for " + str(dict_key) + "is equal to: " + str(percentile_75))

        list_of_quartile_determination = []
        index_of_column_to_edit = NEW_DATA_LIST[0].index(str(dict_key) + ' Quartile')

        # Make labels for each patient
        for i in list_with_gene_expression_score:
            if i < percentile_50:
                list_of_quartile_determination += ['Median']
            elif i > percentile_75:
                list_of_quartile_determination += ['Third_quartile']
            else:
                list_of_quartile_determination += ['Between_50_and_75']

        for i, determination in enumerate(list_of_quartile_determination):
            NEW_DATA_LIST[i+1][index_of_column_to_edit] = determination

        if percentile_50 == percentile_75:
            genes_with_no_changes_in_expression.append(SEQUESTERED_GENES[gene])
        else:
            genes_with_varying_expression.append(SEQUESTERED_GENES[gene])

    # Need to remove empty elements from the header.
    length_of_first_row = len(NEW_DATA_LIST[1])
    NEW_DATA_LIST[0] = NEW_DATA_LIST[0][:length_of_first_row]

    print 'Changes in expression in: ' + str(len(genes_with_varying_expression)) + ' :' + str(genes_with_varying_expression)

    with open('Genes_with_changed_expression', 'wb') as f:
        pickle.dump(genes_with_varying_expression, f)

    print 'No change in expression in: ' + str(len(genes_with_no_changes_in_expression)) + ' :' + str(genes_with_no_changes_in_expression)

    # Time to remove data for genes that do not change their expression
    gene_index = []
    gene_quartile_index = []

    for i in genes_with_no_changes_in_expression:
        gene_index.append(NEW_DATA_LIST[0].index(i))
        gene_quartile_index.append(NEW_DATA_LIST[0].index(str(i + ' Quartile')))

    indexes_to_remove = gene_index + gene_quartile_index

    # Sorting the list in reverse order, so removing by index removes correct entries
    indexes_to_remove = sorted(indexes_to_remove, reverse=True)
    for i in indexes_to_remove:
        for line in NEW_DATA_LIST:
            line.pop(i)

    # Exporting new data file
    with open(str('merged_data_with_quartiles_' + cohort), 'wb') as fp:
        pickle.dump(NEW_DATA_LIST, fp)


def kaplan_meier_plot_and_stat(cohort):

    os.chdir(str(CD + '/' + cohort))
    with open(str('merged_data_with_quartiles_' + cohort), 'rb') as fp:
        data = pickle.load(fp)
    with open('Genes_with_changed_expression', 'rb') as f:
        changed_expression = pickle.load(f)

    # Loading patient/gene data as pandas DataFrame, 1st row will be header
    header = data[0]
    data = data[1:]
    pandas_data = DataFrame(data=data, columns=header, dtype=float)

    try:
        pandas_data = pandas_data.drop(pandas_data[pandas_data.time_alive == 'NA'].index)
    except TypeError:
        pass

    # This is IT! That's how to easily group by 2 columns in pandas!
    '''
    grouped_data = pandas_data.groupby(['gender', 'Wdr65 Quartile'])
    group = grouped_data.get_group(('female', 'Median'))
    test_data = pd.to_numeric(group['time_alive'])
    '''

    def run_survival(data, gene_name):
        ay = plt.subplot(111)
        ay.set_title(gene_name)

        gene = gene_name + ' Quartile'
        gene = ''.join(gene)
        genders = ['male', 'female']
        group_by = ['gender']
        group_by.append(gene)
        # print group_by
        gene_groups = ['Median', 'Third_quartile', 'Between_50_and_75']

        kmf = lifelines.KaplanMeierFitter()

        grouped_data = data.groupby(group_by)

        for gene_group in gene_groups:
            for gender in genders:
                try:
                    pre_tuple_list = [gender, gene_group]
                    group = tuple(pre_tuple_list)

                    # print 'tuple: ' + str(group)
                    d = grouped_data.get_group(group)
                    kaplan_meier_time = pd.to_numeric(d['time_alive'])
                    kaplan_meier_event = d['death_status']

                    kmf.fit(kaplan_meier_time, kaplan_meier_event, label=group)

                    kmf.plot(ax=ay, show_censors=False, ci_show=False)
                except KeyError:
                    # print "No " + str(gender) + ' in gene' + str(gene_group)
                    pass

        event_durations = pandas_data.as_matrix(columns=['time_alive'])

        pandas_data['stat_col'] = pandas_data[gene] + pandas_data['gender']
        group_labels = pandas_data.as_matrix(columns=['stat_col'])

        event = numpy.array(pandas_data.as_matrix(columns=['death_status']))

        result = multivariate_logrank_test(event_durations, group_labels, event, 0.85)

        os.chdir(str(CD + '/' + cohort))
        if not os.path.exists(CD + '/' + cohort + '/results'):
            os.makedirs('results')
        os.chdir(str(CD + '/' + cohort + '/results'))

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.savefig(str(gene_name + '_' + str(result.is_significant) + '.png'), bbox_inches='tight')

        f.close()
        plt.close()

    for gene in changed_expression:
        gene = str(gene)
        run_survival(pandas_data, gene)


def process_data(cohort):
    print "Processing data for: " + cohort
    get_clinical_data(cohort)
    print "Clinical data downloaded"
    get_mRNASeq_data(cohort)
    print "mRNASeq data downloaded"
    check_mRNASeq_data(cohort)
    merge_clinical_with_mrna(cohort)
    print "Files merged"
    label_genes_with_quartiles(cohort)
    print "Quartiles determined. Plotting."
    kaplan_meier_plot_and_stat(cohort)
    print "-" * 45


merge_clinical_with_mrna_new('ACC')


