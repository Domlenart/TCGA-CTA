import fbget
import csv
from multiprocessing.dummy import Pool as ThreadPool


# opens txt file with genes, loads them up into list and splits newline symbol from each one
with open('pure-sequestered-cta') as f:
    list_of_genes = f.readlines()

list_of_genes = [x.strip() for x in list_of_genes]
print list_of_genes


# splits the gene file into parts to feed into fbget (because otherwise it timeouts on too many genes requested from DB)
def chunk(seq, size):

    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


# gets gene log2 data in float format for the current patient for all genes
def get_genes_data(patient_barcode):

    high_exp = 0    # stores how many genes are expressed more that in control

    for part in chunk(list_of_genes, 15):

        gene_data = fbget.mrnaseq(format="tsv", gene=part, barcode=patient_barcode, sample_type='TP')    # getting genetic data

        gene_file = open(patient_barcode, 'w')                                              # making and opening a file
        gene_file.write(gene_data)                                                          # writing data to a file
        gene_file.close()                                                                   # closing file

        reader = csv.DictReader(open(patient_barcode), delimiter='\t')            # opens the file and reads header

        for row in reader:
            if row['expression_log2'] != "None" and float(row['expression_log2']) > 1.5:
                high_exp += 1

    print patient_barcode + '\t' + str(high_exp)
    return high_exp


def clinical_checker(cohort=""):                             # gets patient data from each cohort separately

    y = fbget.clinical(format="tsv", cohort=cohort)          # gets clinical data
    text_file = open(cohort, "w")  # creates the file for the data
    text_file.write(y)  # writes the data into the file
    text_file.close()  # closes the file to prevent corruption
    reader = csv.DictReader(open(cohort), delimiter='\t')  # opens the file in read mode
    list_of_patients = []
    list_of_deaths = []
    list_of_ongoing = []

    for row in reader:
        if row['days_to_death'] != "NA" and row['gender'] == 'female':
            patient_id = row['tcga_participant_barcode']
            death = row['days_to_death']
            print patient_id + '\t' + death
            list_of_patients.append(patient_id)
            list_of_deaths.append(death)
        elif row['days_to_last_followup'] != "NA" and row['gender'] == 'female':
            patient_id = row['tcga_participant_barcode']
            last_followup = row['days_to_last_followup']
            print patient_id + '\t' + last_followup
            list_of_patients.append(patient_id)

    print "_______________________________"

    pool = ThreadPool(8)
    results = pool.map(get_genes_data, list_of_patients)
    pool.close()
    pool.join()


clinical_checker("BRCA")
