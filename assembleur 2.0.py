#Binome :
# Keravis Mael
# Fresneau Abel

""" SOMMAIRE """
# 0 - Presentation
    # lignes -> 1 a 41
# I - Touver le brin d'ADN
    # lignes -> 42 a 100
# II - Informations complementaires sur l ADN
    # lignes -> 101 a 122
# III - L ARN
    # lignes ->123 a 159
# IV - Proteine
    # lignes -> 160 a 222
# V - Informations complementaires sur les proteines
    # lignes -> 223 a 341
# VI - Nos essais/tests infructueux
    # lignes -> 342 a 487
# VII - Credits
    # lignes -> 488 a 511

""" Apres quelques heures de travail, de recherche, d echec, de reussite et de colaboration, nous avons le plaisir de vous presenter notre projet assembleur.
    Il aura comme fonction de trouver un brin d ADN complet a partir de plusieurs brins issus de PCR, la traduire en ARN puis la transcrire en proteine
    (selon le type d'organisme etudie). L'ensemble des resultats sera presente dans une interface tkinter a la fin du programme."""
    
# Importation des programmes utilises :


import timeit
from collections import Counter
# Cet import sert a simplifier notre programme en evitant l ajout d'une boucle while ou for (partie II)

from turtle import * # Pour l amour de la bio-informatique (petit bonus)

from prettytable import PrettyTable # tableau dans les resultats

""" I - Trouver le brin d ADN """

# Plusieurs idees ont ete proposees dans notre binome :
# -> La premiere consistait a prendre un certain nombre de nucleotide choisis manuellement (donc avec une interaction), puis de comparer les sequences deux a deux avec des motifs de la taille imposee.
# -> La deuxieme consistait a entrer manuellement deux sequences et les comparer nucleotide par nucleotide, neanmoins elle ne correspondait pas aux attendus.
# -> La troisieme et derniere idee etait d associer la deuxieme avec une ouverture de dossier, ce qui rend la contingence plus automatique et fluide.
# Nous avons donc naturellement previlegie la derniere methode.

# Nous avons reparti le travail en deux parties,
# Abel a imagine le gros de programme sur la gestion de l ouverture du fichier et de une comparaison minimisant les missmatch et a tenter d incorporer la gestion des gap,
# Mael a reutilise son programme (methode 2) pour simplifier le programme d Abel, majoritairement sur le nom des variables, il a aussi programmer la partie ARN/proteine.
# La partie "information de la proteine" a ete effectue ensemble.

read_liste = open("my_read.txt").read().splitlines() # Ouverture du fichier my_read.txt
read_liste = sorted(read_liste, key=len, reverse=True) # On classe les read par ordre decroissant pour commencer avec la plus longue sequence

seq1=read_liste.pop(0) # Selection de la premiere read a comparer
Final=[] 

####### Premiere fonction du programme #######

# utilite : Comparer les nucleotides de 2 sequences une a une puis supprimer les nucleotides en commun et obtenir le brin d ADN le plus proche de la sequences recherchee.

def comparaison(seq2):  # on cree une fonction qui va servir a comparer une sequence avec notre read principale (read alpha)
    seq2=read_liste[seq2]  # on recupere la read numero seq2 dans la variable eponime
    d=0 # la variable d correspond a notre possition dans la read cible 
    for b in range(len(seq1)): # cette boucle va prendre 1 par 1 les nucleotides de la read alpha 
        d=0 # on se met a 0 dans la read cible (seq2)
        for a in range(len(seq2)): # puis on va l augmenter 1 par 1 dans notre read cible
             # si la nucleotide de l alpha numero d correspond a la nucleotide numero b de la read cible. cette regle sera aussi valable pour la nucleotide suivante (permet le mismatch)
            if((b+d)!=(len(seq1)-1))and((seq2[a]==seq1[b+d])or(seq2[a+1]==seq1[b+d+1])):
                d=d+1 # on peux aller regarder la nucleotide suivante dans la cible 
            else: # si se n est pas le cas :
                break # on retourne a la ligne 69 mais avec une valeur de d = d+1 (comme on a casse cette iterration de la boucle ligne 70)
            if((b+d)==(len(seq1)-1)): # si au contraire on continue dans la boucle ligne 70 et si on arrive avec notre possition de d+b = la longeur de la read alpha ca veux dire on a un recouvrement max
                return(d+1) # on peut donc renvoyer la valeur de recouvrement + 1 soit finalement la position dans la read cible
    return(0) # dans le cas tres rare, voir impossible on renvoie 0

while read_liste:  # tant qu il y a des read dans notre notre liste on va essayer de les ajouter
    for i in range(len(read_liste)): # pour cela on va creer une liste de score avec le recouvrement de chaque liste avec la read alpha actuelle.
        Final.append(comparaison(i)) # on cree cette liste 
    Lebon=read_liste.pop(Final.index(max(Final),0,len(Final))) # on va prendre la read avec le meilleur recouvrement (max(Final)) et la definir comme la prochaine a ajouter 
    seq1=seq1+Lebon[max(Final):] # et nous pouvons l'ajouter a notre read alpha 
    Final=[] # puis on reboot notre liste de score pret pour une nouvelle comparaison (necessaire vu qu on utilise append)

# Le code inclu les missmatch, verif effecuee avec un dossier my_read_missmatch.txt
# Le code n'inclue pas les GAP, le nombre de nucleotide augmente. Plusieurs essais infructueux pour le resoudre (voir partie VI). Test effectue avec un document my_read_gap.txt

print("La sequence d ADN est : ", seq1)

ADN = seq1 # On renomme notre variable pour une meilleure comprehension

""" II - Informations complementaires sur l ADN """

# Apres avoir complete notre programme avec le minimum demande on a voulu le completer avec des informations telles que :
# Le nombre de nucleotide, la composition en A, T,C et G, les %A+T et %C+G car ce sont generalement des infomations utiles a savoir sur une sequence inconnue.

####### Deuxieme fonction du programme #######

# Cette ligne sert de compteur a nucleotide :

A,T,C,G =  [ADN.count('A'),ADN.count('T'), ADN.count('C'), ADN.count('G')]

table = PrettyTable()
table.field_names = ["nucleotides", "nombre"]

table.add_row(["A", A])
table.add_row(["T", T])
table.add_row(["C", C])
table.add_row(["G", G])
table.add_row(["%A+T", (round(((A+T)/len(ADN)*100),2))])
table.add_row(["%G+C", (round(((C+G)/len(ADN)*100),2))])

print(table)

ADN2 = [ADN[i:i+2] for i in range(0, len(ADN), 2)]
d = [ADN2.count('AA'), ADN2.count('AT'), ADN2.count('AC'),
     ADN2.count('AG'), ADN2.count('TA'), ADN2.count('TT'),
     ADN2.count('TC'), ADN2.count('TG'), ADN2.count('CA'),
     ADN2.count('CT'), ADN2.count('CC'), ADN2.count('CG'),
     ADN2.count('GA'), ADN2.count('GT'), ADN2.count('GC'),
     ADN2.count('GG')]

ADN3 = [ADN[i:i+2] for i in range(1, len(ADN), 2)]
e = [ADN3.count('AA'), ADN3.count('AT'), ADN3.count('AC'),
     ADN3.count('AG'), ADN3.count('TA'), ADN3.count('TT'),
     ADN3.count('TC'), ADN3.count('TG'), ADN3.count('CA'),
     ADN3.count('CT'), ADN3.count('CC'), ADN3.count('CG'),
     ADN3.count('GA'), ADN3.count('GT'), ADN3.count('GC'),
     ADN3.count('GG')]


AA, AT, AC, AG, TA, TT, TC, TG, CA, CT, CC, CG, GA, GT, GC, GG = [d[0]+e[0], d[1]+e[1], d[2]+e[2],
                                                                  d[3]+e[3], d[4]+e[4], d[5]+e[5],
                                                                  d[6]+e[6], d[7]+e[7], d[8]+e[8],
                                                                  d[9]+e[9], d[10]+e[10], d[11]+e[11],
                                                                  d[12]+e[12], d[13]+e[13], d[14]+e[14],
                                                                  d[15]+e[15]]

# Créer le tableau et ajouter des en-têtes de colonne
table = PrettyTable()
table.field_names = ["Dinucleotides", "nombre"]

# Ajouter des données au tableau
table.add_row(["AA", AA])
table.add_row(["AT", AT])
table.add_row(["AC", AC])
table.add_row(["AG", AG])
table.add_row(["TA", TA])
table.add_row(["TT", TT])
table.add_row(["TC", TC])
table.add_row(["TG", TG])
table.add_row(["CA", CA])
table.add_row(["CT", CT])
table.add_row(["CC", CC])
table.add_row(["CG", CG])
table.add_row(["GA", GA])
table.add_row(["GT", GT])
table.add_row(["GC", GC])
table.add_row(["GG", GG])

# Personnaliser l'apparence du tableau
table.align = "l"   # aligner le texte à gauche
table.border = True   # ajouter des bordures stylisées
table.header_style = "title"   # styliser les en-têtes de colonne

# Afficher le tableau
print(table)



""" III - L ARN """
    
# Cette partie du programme nous tenait a coeur en tant que futurs biologistes.
# On souhaitait afficher l'ARN dans le programme.
# Mael a voulu ajouter une partie interactive pour demander si la sequence d ADN se situe dans une cellule prokaryote ou eukaryote car plusieurs facteurs sont a prendre en consideration
# Dans le cas de la bacterie, ADN = pre-ARN = ARNm, donc il est simple de trouver la premiere sequence proteique.
# Dans le cas des eukaryotes, ADN = pre-ARN ne correspond pas a l ARN car il y a un epissage.
        # Trois problemes se sont poses :
            # 1 - Nous n avons pas encore les connaissances necessaire en deuxieme annee de biologie pour correctement integrer les epissages.
            # 2 - L epissage alternatif demande un niveau de maitrise trop important de la bio informatique et de la genetique pour nous.
            # 3 - La sequence etudiee n est pas assez longue pour voir si le programme fonctionne correctement (dans le cas ou on serait chez les eukaryotes et qu on integre les epissages).

# Les lignes de codes en dessous servent a rendre le programme interactif, en fonction de la reponse donnee, la reponse affichee varie.

ARN = ADN.replace("T" , "U") # Ici l ADN = ARNm donc les lignes sont plus simples
    
print("La sequence pre-ARN correspond a l ARN, donc la sequence du brin trouve par le programme est la suivante : ", ARN)

""" IV - Proteine """

# On insere un petit tableau de donnee avec les correspondances codons/acides amines

codon_table = {
"ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
"ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
"AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
"AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
"CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
"CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
"CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
"GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
"GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
"TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
"TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
"TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
"TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
}

####### Troisieme fonction du programme #######

# Cette fonction sert a separer les codons 3 par 3 a la suite du codon start, jusqu a un codon stop.

#pourrait être simplifie si on prend ADN et qu on met range(0 1 ou 3, len(ADN), 3)
ORF1 = []
ORF2 = []
ORF3 = []
for i in range(0, len(ADN), 3): # boucle à travers la séquence par groupes de 3
   ORF1 = [ADN[i:i+3] for i in range(0, len(ADN), 3)]
   ORF2 = [ADN[i:i+3] for i in range(1, len(ADN), 3)]
   ORF3 = [ADN[i:i+3] for i in range(2, len(ADN), 3)]

ORF1 = "".join(ORF1[:i])
ORF2 = "".join(ORF2[:i])
ORF3 = "".join(ORF3[:i])

proteine1= []

# Parcourir la séquence par groupes de trois nucléotides
for i in range(0, len(ORF1), 3):
    codon = ORF1[i:i+3]
    if len(codon) == 3:
        proteine1.append(codon_table[codon])

# Convertir la liste d'acides aminés en une chaîne
proteine1 = ''.join(proteine1)

proteine2 = []

# Parcourir la séquence par groupes de trois nucléotides
for i in range(0, len(ORF2), 3):
    codon = ORF2[i:i+3]
    if len(codon) == 3:
        proteine2.append(codon_table[codon])

# Convertir la liste d'acides aminés en une chaîne
proteine2 = ''.join(proteine2)

proteine3 = []

# Parcourir la séquence par groupes de trois nucléotides
for i in range(0, len(ORF3), 3):
    codon = ORF3[i:i+3]
    if len(codon) == 3:
        proteine3.append(codon_table[codon])

# Convertir la liste d'acides aminés en une chaîne
proteine3 = ''.join(proteine3)


def liste_prot(proteinex):
    seq = proteinex

    proteinex = []
    start = "M"
    stop = "_"

    while True:
        start_index = seq.find(start)
        if start_index == -1:
            break

        seq = seq[start_index:]
        stop_index = seq.find(stop)
        if stop_index == -1:
            break

        protein = seq[:stop_index+1]
        seq = seq[stop_index+1:]
    
        proteinex.append(protein)

# Filtrer les séquences de protéines qui ne commencent pas par M ou ne se terminent pas par _
    proteinesx_filtered = [p for p in proteinex if p.startswith(start) and p.endswith(stop)]

# Afficher les séquences de protéines trouvées
    return(proteinesx_filtered)

a=liste_prot(proteine1)
a.sort(key=lambda x: len(x), reverse=True)
b=liste_prot(proteine2)
b.sort(key=lambda x: len(x), reverse=True)
c=liste_prot(proteine3)
c.sort(key=lambda x: len(x), reverse=True)

table = PrettyTable()
table.field_names = ["ORFs", "sequence", "translation possible M (ATG-> stop)"]

# Ajouter des données au tableau
table.add_row(["ORF1", ("".join(ORF1[:i])), a])
table.add_row(["", "",""])
table.add_row(["OFR2", ("".join(ORF2[:i])), b])
table.add_row(["", "",""])
table.add_row(["ORF3", ("".join(ORF3[:i])), c])

table.align = "l"   # aligner le texte à gauche
table.border = True   # ajouter des bordures stylisées
table.header_style = "title"   # styliser les en-têtes de colonne

# Afficher le tableau
print(table)

print("Quelle ORF souhaitez vous étudier ? (seule la plus longue proteine sera etudiee) ")
reponse = input("1, 2 ou 3 ? ")

if reponse == "1" :
    proteine = a[0]
    print("Proteine etudiee : ", a[0])
if reponse == "2" :
    proteine = b[0]
    print("Proteine etudiee : ", b[0])
if reponse == "3" :
    proteine = c[0]
    print("Proteine etudiee : ", c[0])

    
""" V - Informations complementaires sur les proteines """


# On entre une table de donnee complete avec des listes numeriques, binaires et regroupant certaines proprietes [phi, nb C, nb O, nb H, nb N,nb S,
# masse acide amine, Hydroxylation, Methylation, Phosphorylation] avec 1 = oui et 0 = non.
# /!\ La methylation depend du phi. l acetylation est trop compliquee car elle depend des N ter. On ne les a donc pas fait.
# pistes pour pousser plus loins le programme -> Ajouter les N et O Glycosylations si presence de sucres dans l'organisme
#                                             -> Ajouter le phi de la proteine
    
# phi= [1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2]
nb_C = [3, 6, 4, 4, 3, 5, 5, 2, 6, 6, 6, 6, 5, 9, 5, 3, 4, 11, 9, 5]
nb_O = [2, 2, 3, 4, 2, 2, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2]
nb_H = [7, 15, 8, 6, 7, 8, 10, 5, 10, 13, 13, 12, 12, 11, 9, 7, 9,12, 11, 11]
nb_N = [1, 4, 3, 1, 1, 1, 2, 1, 3, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1]
nb_S = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
MM = [89.0935, 174.2017, 132.1184,133.1032, 121.159, 147.1299,146.1451, 75.0669, 155.1552,131.1736, 131.1736, 146.1882,149.2124, 165.19, 115.131,105.093, 119.1197, 204.2262,181.1894, 117.1469]
Hydroxylation = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0]
Methylation = [0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1]
Phosphorylation = [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]
nom = ["Adenine","Arginine","Asparagine",
       "Aspartate", "Cysteine", "Glutamate",
       "Glutamine", "Glycine", "Histidine",
       "Isoleucine", "Leucine", "Lysine",
       "Methionine", "Phenylalanine","Proline",
       "Serine", "Threonine", "Tryptophane",
       "Tyrosine","Valine"]
# Cette fonction donne le nombre de chaque AA de la sequence proteique (par type d acide amine)
# nombre de chaque acide amine dans la proteine trouvee

a = [proteine.count('A'),proteine.count('R'), proteine.count('N'),
     proteine.count('D'), proteine.count('C'), proteine.count('E'),
     proteine.count('Q'), proteine.count('G'), proteine.count('H'),
     proteine.count('I'), proteine.count('L'), proteine.count('K'),
     proteine.count('M'), proteine.count('F'), proteine.count('P'),
     proteine.count('S'), proteine.count('T'), proteine.count('W'),
     proteine.count('Y'), proteine.count('V')]

A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V = [proteine.count('A'),proteine.count('R'), proteine.count('N'),proteine.count('D'), proteine.count('C'), proteine.count('E'), proteine.count('Q'), proteine.count('G'), proteine.count('H'),proteine.count('I'), proteine.count('L'), proteine.count('K'),proteine.count('M'), proteine.count('F'), proteine.count('P'),proteine.count('S'), proteine.count('T'), proteine.count('W'),proteine.count('Y'), proteine.count('V')]

i =0
# phi = 0
b=[]
nb_Carbone = 0
nb_Hydrogene = 0
nb_Oxygene = 0
nb_Azote = 0
nb_Soufre = 0
Masse_Molaire = 0
Hydroxylation_ = 0
Methylation_ = 0
Phosphorylation_ = 0
table = PrettyTable()
table.field_names = ["Acide aminé", "nombre", "pourcentage"]


for i in range(19): # Ici on va calculer le nombre de chaque element en fonction du nombre d acide amine, i correspond aux acides amines dans la liste a
    nb_Carbone = nb_Carbone + a[i]*nb_C[i]
    nb_Hydrogene = nb_Hydrogene + a[i]*nb_H[i]
    nb_Oxygene = nb_Oxygene + a[i]*nb_O[i]
    nb_Azote = nb_Azote + a[i]*nb_N[i]
    nb_Soufre = nb_Soufre + a[i]*nb_S[i]
    Masse_Molaire = Masse_Molaire + a[i]*MM[i] # sans les laisons peptidiques, on les enlevera ligne 299
    Hydroxylation_ = Hydroxylation_ + a[i]*Hydroxylation[i] 
    Methylation_ = Methylation_ + a[i]*Methylation[i]
    Phosphorylation_ = Phosphorylation_ + a[i]*Phosphorylation[i]
"""
        table = PrettyTable()
        table.field_names = ["Acide aminé", "nombre", "pourcentage", "MM", "n C", "n H", "n O", "n N", "n S", "nb Methylation", "nb Methylation", "nb Phosphorylation"]
        table.add_row([nom[i], a[i], ((a[i]/len(proteine))*100), (Masse_Molaire + a[i]*MM[i]), (nb_Carbone + a[i]*nb_C[i]), (nb_Hydrogene + a[i]*nb_H[i]), (nb_Oxygene + a[i]*nb_O[i]),(nb_Azote + a[i]*nb_N[i]), (nb_Soufre + a[i]*nb_S[i]), (Hydroxylation_ + a[i]*Hydroxylation[i]), (Methylation_ + a[i]*Methylation[i]), (Phosphorylation_ + a[i]*Phosphorylation[i])])
        table.add_row(["Total", len(proteine), 100, (Masse_Molaire - 18*len(proteine)), nb_Carbone, nb_Hydrogene, nb_Oxygene, nb_Azote, nb_Soufre, Hydroxylation_, Methylation_, Phosphorylation_])
        table.align = "l"   # aligner le texte à gauche
        table.border = True   # ajouter des bordures stylisées
        table.header_style = "title"   # styliser les en-têtes de colonne
        print(table) #ValueError: Field names must be unique
        
        table.add_row(["Adenine", a[0], (a[0]/len(proteine))*100], Masse_Molaire + a[0]*MM[0], (nb_Carbone + a[0]*nb_C[0]), (nb_Hydrogene + a[0]*nb_H[0]), (nb_Oxygene + a[0]*nb_O[0]),(nb_Azote + a[0]*nb_N[0]), (nb_Soufre + a[0]*nb_S[0]), (Hydroxylation_ + a[0]*Hydroxylation[0]), (Methylation_ + a[0]*Methylation[0]), (Phosphorylation_ + a[0]*Phosphorylation[0])]))
        table.add_row(["Arginine", a[1], (a[1]/len(proteine))*100], Masse_Molaire + a[1]*MM[1], (nb_Carbone + a[1]*nb_C[1]), (nb_Hydrogene + a[1]*nb_H[1]), (nb_Oxygene + a[1]*nb_O[1]),(nb_Azote + a[1]*nb_N[1]), (nb_Soufre + a[1]*nb_S[1]), (Hydroxylation_ + a[1]*Hydroxylation[1]), (Methylation_ + a[1]*Methylation[1]), (Phosphorylation_ + a[1]*Phosphorylation[1])]))
        table.add_row(["Asparagine", a[2], (a[2]/len(proteine))*100], Masse_Molaire + a[2]*MM[2], (nb_Carbone + a[2]*nb_C[2]), (nb_Hydrogene + a[2]*nb_H[2]), (nb_Oxygene + a[2]*nb_O[2]),(nb_Azote + a[2]*nb_N[2]), (nb_Soufre + a[2]*nb_S[2]), (Hydroxylation_ + a[2]*Hydroxylation[2]), (Methylation_ + a[2]*Methylation[2]), (Phosphorylation_ + a[2]*Phosphorylation[2])]))
        table.add_row(["Aspartate", a[3], (a[3]/len(proteine))*100], Masse_Molaire + a[3]*MM[3], (nb_Carbone + a[3]*nb_C[3]), (nb_Hydrogene + a[3]*nb_H[3]), (nb_Oxygene + a[3]*nb_O[3]),(nb_Azote + a[3]*nb_N[3]), (nb_Soufre + a[3]*nb_S[3]), (Hydroxylation_ + a[3]*Hydroxylation[3]), (Methylation_ + a[3]*Methylation[3]), (Phosphorylation_ + a[3]*Phosphorylation[3])]))
        table.add_row(["Cysteine", a[4], (a[4]/len(proteine))*100], Masse_Molaire + a[4]*MM[4], (nb_Carbone + a[4]*nb_C[4]), (nb_Hydrogene + a[4]*nb_H[4]), (nb_Oxygene + a[4]*nb_O[4]),(nb_Azote + a[4]*nb_N[4]), (nb_Soufre + a[4]*nb_S[4]), (Hydroxylation_ + a[4]*Hydroxylation[4]), (Methylation_ + a[4]*Methylation[4]), (Phosphorylation_ + a[4]*Phosphorylation[4])]))
        table.add_row(["Glutamate", a[5], (a[5]/len(proteine))*100], Masse_Molaire + a[5]*MM[5], (nb_Carbone + a[5]*nb_C[5]), (nb_Hydrogene + a[5]*nb_H[5]), (nb_Oxygene + a[5]*nb_O[5]),(nb_Azote + a[5]*nb_N[5]), (nb_Soufre + a[5]*nb_S[5]), (Hydroxylation_ + a[5]*Hydroxylation[5]), (Methylation_ + a[5]*Methylation[5]), (Phosphorylation_ + a[5]*Phosphorylation[5])]))
        table.add_row(["Glutamine", a[6], (a[6]/len(proteine))*100], Masse_Molaire + a[6]*MM[6], (nb_Carbone + a[6]*nb_C[6]), (nb_Hydrogene + a[6]*nb_H[6]), (nb_Oxygene + a[6]*nb_O[6]),(nb_Azote + a[6]*nb_N[6]), (nb_Soufre + a[6]*nb_S[6]), (Hydroxylation_ + a[6]*Hydroxylation[6]), (Methylation_ + a[6]*Methylation[6]), (Phosphorylation_ + a[6]*Phosphorylation[6])]))
        table.add_row(["Glycine", a[7], (a[7]/len(proteine))*100], Masse_Molaire + a[7]*MM[7], (nb_Carbone + a[7]*nb_C[7]), (nb_Hydrogene + a[7]*nb_H[7]), (nb_Oxygene + a[7]*nb_O[7]),(nb_Azote + a[7]*nb_N[7]), (nb_Soufre + a[7]*nb_S[7]), (Hydroxylation_ + a[7]*Hydroxylation[7]), (Methylation_ + a[7]*Methylation[7]), (Phosphorylation_ + a[7]*Phosphorylation[7])]))
        table.add_row(["Histidine", a[8], (a[8]/len(proteine))*100], Masse_Molaire + a[8]*MM[8], (nb_Carbone + a[8]*nb_C[8]), (nb_Hydrogene + a[8]*nb_H[8]), (nb_Oxygene + a[8]*nb_O[8]),(nb_Azote + a[8]*nb_N[8]), (nb_Soufre + a[8]*nb_S[8]), (Hydroxylation_ + a[8]*Hydroxylation[8]), (Methylation_ + a[8]*Methylation[8]), (Phosphorylation_ + a[8]*Phosphorylation[8])]))
        table.add_row(["Isoleucine", a[9], (a[9]/len(proteine))*100], Masse_Molaire + a[9]*MM[9], (nb_Carbone + a[9]*nb_C[9]), (nb_Hydrogene + a[9]*nb_H[9]), (nb_Oxygene + a[9]*nb_O[9]),(nb_Azote + a[9]*nb_N[9]), (nb_Soufre + a[9]*nb_S[9]), (Hydroxylation_ + a[9]*Hydroxylation[9]), (Methylation_ + a[9]*Methylation[9]), (Phosphorylation_ + a[9]*Phosphorylation[9])]))
        table.add_row(["Leucine", a[10], (a[10]/len(proteine))*100], Masse_Molaire + a[10]*MM[10], (nb_Carbone + a[10]*nb_C[10]), (nb_Hydrogene + a[10]*nb_H[10]), (nb_Oxygene + a[10]*nb_O[10]),(nb_Azote + a[10]*nb_N[10]), (nb_Soufre + a[10]*nb_S[10]), (Hydroxylation_ + a[10]*Hydroxylation[10]), (Methylation_ + a[10]*Methylation[10]), (Phosphorylation_ + a[10]*Phosphorylation[10])]))
        table.add_row(["Lysine", a[11], (a[11]/len(proteine))*100], Masse_Molaire + a[11]*MM[11], (nb_Carbone + a[11]*nb_C[11]), (nb_Hydrogene + a[11]*nb_H[11]), (nb_Oxygene + a[11]*nb_O[11]),(nb_Azote + a[11]*nb_N[11]), (nb_Soufre + a[11]*nb_S[11]), (Hydroxylation_ + a[11]*Hydroxylation[11]), (Methylation_ + a[11]*Methylation[11]), (Phosphorylation_ + a[11]*Phosphorylation[11])]))
        table.add_row(["Methionine", a[12], (a[12]/len(proteine))*100], Masse_Molaire + a[12]*MM[12], (nb_Carbone + a[12]*nb_C[12]), (nb_Hydrogene + a[12]*nb_H[12]), (nb_Oxygene + a[12]*nb_O[12]),(nb_Azote + a[12]*nb_N[12]), (nb_Soufre + a[12]*nb_S[12]), (Hydroxylation_ + a[12]*Hydroxylation[12]), (Methylation_ + a[12]*Methylation[12]), (Phosphorylation_ + a[12]*Phosphorylation[12])]))
        table.add_row(["Phenylalanine", a[13], (a[13]/len(proteine))*100], Masse_Molaire + a[13]*MM[13], (nb_Carbone + a[13]*nb_C[13]), (nb_Hydrogene + a[13]*nb_H[13]), (nb_Oxygene + a[13]*nb_O[13]),(nb_Azote + a[13]*nb_N[13]), (nb_Soufre + a[13]*nb_S[13]), (Hydroxylation_ + a[13]*Hydroxylation[13]), (Methylation_ + a[13]*Methylation[13]), (Phosphorylation_ + a[13]*Phosphorylation[13])]))
        table.add_row(["Proline", a[14], (a[14]/len(proteine))*100], Masse_Molaire + a[14]*MM[14], (nb_Carbone + a[14]*nb_C[14]), (nb_Hydrogene + a[14]*nb_H[14]), (nb_Oxygene + a[14]*nb_O[14]),(nb_Azote + a[14]*nb_N[14]), (nb_Soufre + a[14]*nb_S[14]), (Hydroxylation_ + a[14]*Hydroxylation[14]), (Methylation_ + a[14]*Methylation[14]), (Phosphorylation_ + a[14]*Phosphorylation[14])]))
        table.add_row(["Serine", a[15], (a[15]/len(proteine))*100], Masse_Molaire + a[15]*MM[15], (nb_Carbone + a[15]*nb_C[15]), (nb_Hydrogene + a[15]*nb_H[15]), (nb_Oxygene + a[15]*nb_O[15]),(nb_Azote + a[15]*nb_N[15]), (nb_Soufre + a[15]*nb_S[15]), (Hydroxylation_ + a[15]*Hydroxylation[15]), (Methylation_ + a[15]*Methylation[15]), (Phosphorylation_ + a[15]*Phosphorylation[15])]))
        table.add_row(["Threonine", a[16], (a[16]/len(proteine))*100], Masse_Molaire + a[16]*MM[16], (nb_Carbone + a[16]*nb_C[16]), (nb_Hydrogene + a[16]*nb_H[16]), (nb_Oxygene + a[16]*nb_O[16]),(nb_Azote + a[16]*nb_N[16]), (nb_Soufre + a[16]*nb_S[16]), (Hydroxylation_ + a[16]*Hydroxylation[16]), (Methylation_ + a[16]*Methylation[16]), (Phosphorylation_ + a[16]*Phosphorylation[16])]))
        table.add_row(["Tryptophane", a[17], (a[17]/len(proteine))*100], Masse_Molaire + a[17]*MM[17], (nb_Carbone + a[17]*nb_C[17]), (nb_Hydrogene + a[17]*nb_H[17]), (nb_Oxygene + a[17]*nb_O[17]),(nb_Azote + a[17]*nb_N[17]), (nb_Soufre + a[17]*nb_S[17]), (Hydroxylation_ + a[17]*Hydroxylation[17]), (Methylation_ + a[17]*Methylation[17]), (Phosphorylation_ + a[17]*Phosphorylation[17])]))
        table.add_row(["Tyrosine", a[18], (a[18]/len(proteine))*100], Masse_Molaire + a[18]*MM[18], (nb_Carbone + a[18]*nb_C[18]), (nb_Hydrogene + a[18]*nb_H[18]), (nb_Oxygene + a[18]*nb_O[18]),(nb_Azote + a[18]*nb_N[18]), (nb_Soufre + a[18]*nb_S[18]), (Hydroxylation_ + a[18]*Hydroxylation[18]), (Methylation_ + a[18]*Methylation[18]), (Phosphorylation_ + a[18]*Phosphorylation[18])]))
        table.add_row(["Valine", a[19], (a[19]/len(proteine))*100], Masse_Molaire + a[19]*MM[19], (nb_Carbone + a[19]*nb_C[19]), (nb_Hydrogene + a[19]*nb_H[19]), (nb_Oxygene + a[19]*nb_O[19]),(nb_Azote + a[19]*nb_N[19]), (nb_Soufre + a[19]*nb_S[19]), (Hydroxylation_ + a[19]*Hydroxylation[19]), (Methylation_ + a[19]*Methylation[19]), (Phosphorylation_ + a[19]*Phosphorylation[19])]))
        table.add_row(["Total", len(proteine), 100, (Masse_Molaire - 18*len(proteine)), nb_Carbone, nb_Hydrogene, nb_Oxygene, nb_Azote, nb_Soufre, Hydroxylation_, Methylation_, Phosphorylation_])
"""

table.add_row(["Adenine", A, (round(((A/len(proteine))*100),2))])
table.add_row(["Arginine", R, (round(((R/len(proteine))*100),2))])
table.add_row(["Asparagine", N, (round(((N/len(proteine))*100),2))])
table.add_row(["Aspartate", D, (round(((D/len(proteine))*100),2))])
table.add_row(["Cysteine", C, (round(((C/len(proteine))*100),2))])
table.add_row(["Glutamate", Q, (round(((Q/len(proteine))*100),2))])
table.add_row(["Glutamine", G, (round(((G/len(proteine))*100),2))])
table.add_row(["Histidine", H, (round(((H/len(proteine))*100),2))])
table.add_row(["Isoleucine", I, (round(((I/len(proteine))*100),2))])
table.add_row(["Leucine", L, (round(((L/len(proteine))*100),2))])
table.add_row(["Lysine", K, (round(((K/len(proteine))*100),2))])
table.add_row(["Methionine", M, (round(((M/len(proteine))*100),2))])
table.add_row(["Phenylalanine", F, (round(((F/len(proteine))*100),2))])
table.add_row(["Proline", P, (round(((P/len(proteine))*100),2))])
table.add_row(["Serine", S, (round(((S/len(proteine))*100),2))])
table.add_row(["Threonine", T, (round(((T/len(proteine))*100),2))])
table.add_row(["Tryptophane", W, (round(((W/len(proteine))*100),2))])
table.add_row(["Tyrosine", Y, (round(((Y/len(proteine))*100),2))])
table.add_row(["Valine", V, (round(((V/len(proteine))*100),2))])
       
table.align = "l"   # aligner le texte à gauche
table.border = True   # ajouter des bordures stylisées
table.header_style = "title"   # styliser les en-têtes de colonne

print(table)

table = PrettyTable()
table.field_names = ["composition atomes", "nombre"]

table.add_row(["nb Carbon", nb_Carbone])
table.add_row(["nb Hydrogen", nb_Hydrogene])
table.add_row(["nb Oxygene", nb_Oxygene])
table.add_row(["nb Azote", nb_Azote])
table.add_row(["nb Soufre", nb_Soufre])
table.add_row(["Masse molaire en g/mol", (round((Masse_Molaire),2))])

table.align = "l"   # aligner le texte à gauche
table.border = True   # ajouter des bordures stylisées
table.header_style = "title"   # styliser les en-têtes de colonne

print(table)

table = PrettyTable()
table.field_names = ["PTM", "nombre"]
table.add_row(["Hydroxylation", Hydroxylation_])
table.add_row(["Methylation", Methylation_])
table.add_row(["Phosphorylation", Phosphorylation_])

table.align = "l"   # aligner le texte à gauche
table.border = True   # ajouter des bordures stylisées
table.header_style = "title"   # styliser les en-têtes de colonne

print(table)

e260 = round(((A * 15.4 + C * 7.5 + G * 11.7 + T * 9.2 + T * 15.4 + G * 7.5 + C * 11.7 + A  * 9.2) * 0.9 * 1000),2)
print ("coefficient d'extinction molaire de l'ADN (A260) : ", e260, "mol/L/cm")
e280 = (a[17]*5500 + a[18]*1490 + a[4]*125)
print ("coefficient d'extinction molaire de de la proteine (A280) : ", e280, "mol/L/cm")
  


print("Fin du programme")
print("Merci d'avoir participe a ce projet")

# Banque a idee
# Toutes les infos de la proteine dans une def et faire un Tkinter avec la def ? Comme ca on peut mettre toutes les proteines
# Ajouter la matrice pour les Gap et missmatchs

def grow1 (population) :
    for i in range (0,len(population), 1) :
        if population[i] > 1:
            population[i] -= 1
        else :
            population.append(9)
            population[i] = 7
    return(population)

bacto = [3,4,3,1,2]
print(bacto)
print(grow1(bacto))
print(grow1(grow1(bacto)))