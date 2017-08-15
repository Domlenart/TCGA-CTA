import csv
import transposer

# transposer.transpose(i='KIRCclin.txt', o='KIRCclinTRANSPOSED.txt', d='\t')
# transposer.transpose(i="KIRCrna.txt", o= 'KIRCrnaTRANSPOSED.txt', d='\t')

clinical_reader = csv.DictReader(open("KIRCclinTRANSPOSED.txt"), delimiter='\t')
rna_reader = csv.DictReader(open("KIRCrnaTRANSPOSED.txt"), delimiter='\t')
next(rna_reader)
males = 0
females = 0

for row in clinical_reader:
    if row['patient.gender'] == 'male':
        males += 1
        barcode = str(row['patient.bcr_patient_barcode'])
        barcode = barcode.upper()
        print barcode
        rna_reader = csv.DictReader(open("KIRCrnaTRANSPOSED.txt"), delimiter='\t')
        for rna_row in rna_reader:
            if rna_row["Hybridization REF"][:12] == barcode:
                print rna_row['ARMC5|79798']

    else:
        females += 1
        barcode = str(row['patient.bcr_patient_barcode'])
        barcode = str(row['patient.bcr_patient_barcode'])
        barcode = barcode.upper()
        print barcode
        rna_reader = csv.DictReader(open("KIRCrnaTRANSPOSED.txt"), delimiter='\t')
        for rna_row in rna_reader:
            if rna_row["Hybridization REF"][:12] == barcode:
                print rna_row['ARMC5|79798']


print 'males:' + str(males) + '\t' + 'females:' + str(females)