# annot
Automatise les premières étapes de Annotathon, output les résultats bruts demandé sous forme de fichier txt, un tableau récapitulant toutes les informations (nombres d'orf, coordonnées, nombre de nucleotide, d'acide aminé, nombre de résultat blastp contre NR, et pour l'ORF choisi le meilleur alignement -pourcentage d'identité, longueur alignement, evalue + nombre de résultat blastp contre SW)
et un graph des résutats de blastp pour un ORF choisi afin d'aider à choisir une valeur seuil ![Screenshot](/HMP1_7827070/ORF3/blastNRplot.pdf)

#### Protocole 
getorf -minsize 180 entre deux codons stop > blastp evalue min 1e-10 contre NR et SW 

#### Déroulement : 

Enregistrer une sequence issue de Annotathon dans seq.fasta 

lancer bio.py 

Un dossier nommé d'aprés l'id de la séquence est crée. Ex: ./HMP1_7827060

Les résultats bruts de getorf et du blastp de tout les orfs sont enregistrés Ex : ./HMP1_7827060/ORF_HMP1_7827060.fasta et ./HMP1_7827060/blastp_HMP1_7827060

Un tableau récapitulant les informations est imprimé to stdout. 

On vous demande de choisir un ORF à étudié à partir de maintenant.

Un sous dossier portant le nom de l'ORF choisi est crée Ex : /HMP1_7827060/ORF3

L'orf est enregistré dans un fichier au format fasta portant son nom Ex: ./HMP1_7827060/ORF3/ORF3.fasta

Cet ORF est blastp contre swissprot avec une evalue de 1e-10 et les résultats sont enregistrés Ex :./HMP1_7827060/ORF3/blastsw

Les résultats du blastp pour l'ORF sélectionner sont mis sous forme de graphique. Ex :./HMP1_7827070/ORF3/blastNRplot.pdf


Le tableau est enregistré /HMP1_7827060/ORF3/tableau.csv



### Requis : 

GNU/Linux

suite ebi 

suite blast+


### Python : 

pandas

matplotlib

BioPython (SeqIO)
