from Bio.Seq import *
from Bio.Restriction import *
from Bio import SeqIO
import numpy
import json

def optim():
    codbias=[
    {'GCT':0.19,'GCC':0.25,'GCA':0.22,'GCG':0.34},
    {'CGT':0.42,'CGC':0.37,'CGA':0.05,'CGG':0.08,'AGA':0.04,'AGG':0.03},
    {'AAT':0.39,'AAC':0.61},
    {'GAT':0.59,'GAC':0.41},
    {'TGT':0.43,'TGC':0.57},
    {'CAA':0.31,'CAG':0.69},
    {'GAA':0.7,'GAG':0.3},
    {'GGT':0.38,'GGC':0.4,'GGA':0.09,'GGG':0.13},
    {'CAT':0.52,'CAC':0.48},
    {'ATT':0.47,'ATC':0.46,'ATA':0.07},
    {'TTA':0.11,'TTG':0.11,'CTT':0.10,'CTC':0.10,'CTA':0.03,'CTG':0.55},
    {'AAA':0.76,'AAG':0.24},
    {'ATG':1.0},
    {'TTT':0.51,'TTC':0.49},
    {'CCC':0.1,'CCA':0.2,'CCT':0.16,'CCG':0.55},
    {'TCT':0.19,'TCC':0.17,'TCA':0.12,'TCG':0.13,'AGT':0.13,'AGC':0.27},
    {'ACT':0.21,'ACC':0.43,'ACA':0.3,'ACG':0.23},
    {'TGG':1.},
    {'TAT':0.53,'TAC':0.47},
    {'GTT':0.29,'GTC':0.2,'GTA':0.17,'GTG':0.34},
    {'TGA':0.3,'TAA':0.62,'TAG':0.09}]
    stdnacodtable={
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCC':'P','CCA':'P','CCT':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N',
    'AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C',
    'TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R',
    'AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G','TGA':'*','TAA':'*','TAG':'*',
    'AUG':'M'
    }
    diccod={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,
    'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'*':20,}
    def rarecod(org):
        rare=[]
        for dic in org:
            for cod,bias in dic.items():
                if bias<0.1:
                    rare.append(cod)
        return rare

    def abcod(aa, codbias):
        aaindex=diccod[aa]
        aacodons=codbias[aaindex]
        inverse = [(value, key) for key, value in aacodons.items()]
        return max(inverse)[1]
    seq=''.join(json.load(open('data.json')))
    rare=rarecod(codbias)
    list=[seq[i:i+3] for i in range(0, len(seq), 3)]
    for c in rare:
        if c in list:
            loc=list.index(c)
            prot=stdnacodtable[c]
            list[loc]=abcod(prot,codbias)
    list=''.join(list)
    json.dump(list, open('data.json','w'))

def optimcod(cod):
    diccod={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,
    'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'*':20,}
    stdnacodtable={
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCC':'P','CCA':'P','CCT':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N',
    'AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C',
    'TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R',
    'AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G','TGA':'*','TAA':'*','TAG':'*',
    'AUG':'M'
    }
    codbias=[
    {'GCT':0.19,'GCC':0.25,'GCA':0.22,'GCG':0.34},
    {'CGT':0.42,'CGC':0.37,'CGA':0.05,'CGG':0.08,'AGA':0.04,'AGG':0.03},
    {'AAT':0.39,'AAC':0.61},
    {'GAT':0.59,'GAC':0.41},
    {'TGT':0.43,'TGC':0.57},
    {'CAA':0.31,'CAG':0.69},
    {'GAA':0.7,'GAG':0.3},
    {'GGT':0.38,'GGC':0.4,'GGA':0.09,'GGG':0.13},
    {'CAT':0.52,'CAC':0.48},
    {'ATT':0.47,'ATC':0.46,'ATA':0.07},
    {'TTA':0.11,'TTG':0.11,'CTT':0.10,'CTC':0.10,'CTA':0.03,'CTG':0.55},
    {'AAA':0.76,'AAG':0.24},
    {'ATG':1.0},
    {'TTT':0.51,'TTC':0.49},
    {'CCC':0.1,'CCA':0.2,'CCT':0.16,'CCG':0.55},
    {'TCT':0.19,'TCC':0.17,'TCA':0.12,'TCG':0.13,'AGT':0.13,'AGC':0.27},
    {'ACT':0.21,'ACC':0.43,'ACA':0.3,'ACG':0.23},
    {'TGG':1.},
    {'TAT':0.53,'TAC':0.47},
    {'GTT':0.29,'GTC':0.2,'GTA':0.17,'GTG':0.34},
    {'TGA':0.3,'TAA':0.62,'TAG':0.09}]
    def abcod(aa):
        aaindex=diccod[aa]
        aacodons=codbias[aaindex]
        inverse = [(value, key) for key, value in aacodons.items()]
        if aa=='W':
            return 'TGG'
        elif aa=='M':
            return 'ATG'
        else:
            return max(inverse)[1]
    count=1
    seqcod=Seq(cod)
    index= diccod[seqcod.translate()]
    listedseq= ''.join(json.load(open('data.json')))
    list=[listedseq[i:i+3] for i in range(0, len(listedseq), 3)]
    if cod=='ATG':
        return cod
    elif cod=="TGG":
        return cod
    else:
        while count:
            if cod == abcod(seqcod.translate()):
                del codbias[index][cod]
            elif len(codbias[index])==0:
                codbias= codbias
            else:
                abundant=abcod(seqcod.translate())
                for c in list:
                    if cod in list:
                        loc=list.index(cod)
                        list[loc]=abundant
                count-=1

    joined=[''.join(list)]
    json.dump(joined, open('data.json','w'))

def ecoriswitch():
    listedseq=''.join(json.load(open('data.json')))
    seq= Seq(listedseq)
    sites=EcoRI.search(seq)
    coddist=[listedseq[i:i+3] for i in range(0, len(listedseq), 3)]
    for item in sites:
        if (item % 3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='ATC'
        elif ((item-1)%3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='AAC'
        elif ((item-2)%3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='GAG'
    joinedstring=[''.join(str(x) for x in coddist)]
    json.dump(joinedstring, open('data.json','w'))

def pstiswitch():
    listedseq=''.join(json.load(open('data.json')))
    seq= Seq(listedseq)
    sites=PstI.search(seq)
    sitesm = [x - 3 for x in sites]
    coddist=[listedseq[i:i+3] for i in range(0, len(listedseq), 3)]
    for item in sitesm:
        if ((item) % 3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='CAA'
        elif ((item-1)%3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='GCG'
        elif ((item-2)%3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='TGT'
    joinedstring=[''.join(str(x) for x in coddist)]
    json.dump(joinedstring, open('data.json','w'))

def speiswitch():
    listedseq=''.join(json.load(open('data.json')))
    seq=Seq(listedseq)
    sites=SpeI.search(seq)
    sitesm = [x - 1 for x in sites]
    coddist=[listedseq[i:i+3] for i in range(0, len(listedseq), 3)]
    for item in sitesm:
        if ((item) % 3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='CTG'
        elif ((item-1)%3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='ACC'
        elif ((item-2)%3)==0:
            ind=int(numpy.floor((item/3)))
            if coddist[ind]=='AAC':
                coddist[ind]='AAT'
            elif coddist[ind]=='TAC':
                coddist[ind]='TAT'
            elif coddist[ind]=='CAC':
                coddist[ind]='CAT'
            elif coddist[ind]=='GAC':
                coddist[ind]='GAT'
    joinedstring=[''.join(str(x) for x in coddist)]
    json.dump(joinedstring, open('data.json','w'))



def xbaiswitch():
    listedseq=''.join(json.load(open('data.json')))
    seq= Seq(listedseq)
    sites=XbaI.search(seq)
    sitesm = [x - 1 for x in sites]
    coddist=[listedseq[i:i+3] for i in range(0, len(listedseq), 3)]
    for item in sitesm:
        if ((item) % 3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='CTG'
        elif ((item-1)%3)==0:
            ind=int(numpy.floor((item/3)))
            coddist[ind]='AGC'
        elif ((item-2)%3)==0:
            ind=int(numpy.floor((item/3)))
            if coddist[ind]=='ATC':
                coddist[ind]='ATT'
            elif coddist[ind]=='TTC':
                coddist[ind]='TTT'
            elif coddist[ind]=='CTC':
                coddist[ind]='CTG'
            elif coddist[ind]=='GTC':
                coddist[ind]='GTG'

    joinedstring=[''.join(str(x) for x in coddist)]
    json.dump(joinedstring, open('data.json','w'))

def breakseq():
    seq=''.join(json.load(open('data.json')))
    return [seq[i:i+3] for i in range(0, len(seq), 3)]

def joinseq(list):
    return ''.join(list)

def criarfasta(titulo,nome):
    seq=''.join(json.load(open('data.json')))
    arquivo = open(nome+'.fasta', 'w')
    arquivo.write('>'+titulo+'\n')
    arquivo.write(seq)
    arquivo.close()

def lerfasta(arq):
    cabecalho=[]
    sequencias=[]
    for fasta in SeqIO.parse(arq, "fasta"):
        cabecalho.append(fasta.id)
        sequencias.append(fasta.seq)
    return [cabecalho, sequencias]

def hairpincheckseq():
    stringedrna=str(Seq(''.join(json.load(open('data.json')))).transcribe())[:40]
    arquivo = open('rna.fasta', 'w')
    arquivo.write('>'+'40 primeiros PBs para serem analisados pelo RNASTRUCTURE'+'\n')
    arquivo.write(stringedrna)
    arquivo.close()
