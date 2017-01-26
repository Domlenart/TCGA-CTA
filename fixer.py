with open('common-acrosome-cta') as f:
    list_of_seq_genes = f.readlines()
list_of_seq_genes = [x.strip() for x in list_of_seq_genes]
print list_of_seq_genes

with open('common-RB-cta') as f2:
    list_of_non_seq_genes = f2.readlines()
list_of_non_seq_genes = [x.strip() for x in list_of_non_seq_genes]
print list_of_non_seq_genes

overlap = []
count = 0
for x in list_of_seq_genes:
    for y in list_of_non_seq_genes:
        if x == y:
            overlap.append(x)
            count += 1
            list_of_seq_genes.remove(x)

print overlap
print count
print list_of_seq_genes

with open('pure-sequestered-cta', 'w') as f3:
    for i in list_of_seq_genes:
        f3.write("%s \n" % i)

with open('overlapping-acrosome-CD', 'w') as f4:
    for i in overlap:
        f4.write("%s \n" % i)
