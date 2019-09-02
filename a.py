from Bio import SeqIO
sequencias=[]
count=0
for fasta in SeqIO.parse("Sequencia.fasta", "fasta"):
    sequencias.append(fasta.seq)

a=str(sequencias)
b=[a[i:i+3] for i in range(0, len(a), 3)]
for i in b:
    if i=='TAG' or i=='CGA' or i=='CGG' or i=='AGG' or i=='TAG' or i=='GGA' or i=='ATA' or i=='CTA':
        count+=1
print(count)
