# annot
Automatise les premières étapes de Annotathon, output les résultats bruts demandé sous forme de fichier txt, un tableau récapitulant toutes les informations (Ex : ![tableau](/HMP1_7827070/ORF3/tableau.csv)
 )
et un graph (Ex : ![Plot](/HMP1_7827070/ORF3/blastNRplot.pdf)
 ) des résutats de blastp pour un ORF choisi afin d'aider à choisir une valeur seuil 

#### Protocole 
getorf -minsize 180 entre deux codons stop | blastp evalue min 1e-10, 5000 resultats max, contre NR et swissprot 

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

Les résultats du blastp pour l'ORF sélectionné sont mis sous forme de graphique. Ex :./HMP1_7827070/ORF3/blastNRplot.pdf


Le tableau est enregistré /HMP1_7827060/ORF3/tableau.csv



### Requis : 

Unix filesystem

GNU grep (perl regex)   

ebi getorf 

blast+


#### Python : 

pandas

matplotlib

BioPython (SeqIO)
