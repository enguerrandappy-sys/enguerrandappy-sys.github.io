# CC1 ADM 
## Enguerrand Appy-Ruat
Ce script s'appuie sur le "DADA2 Pipeline Tutorial (1.16) et l'analyse des données MiSeq_SOP y étant associées.

## Chargement du package dada2
```{r}
library(dada2)
```
## Vérification de la version du package 
```{r}
packageVersion("dada2")
```
## Ouverture du fichier de données 
```{r}
path <- "~/MiSeq_SOP" 
list.files(path)
```
##[1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
##[3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
##[5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
##[7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
##[9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
##[11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
##[13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
##[15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
##[17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
##[19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
##[21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
##[23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
##[25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
##[27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
##[29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
##[31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
##[33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
##[35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
##[37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
##[39] "filtered"                      "HMP_MOCK.v35.fasta"           
##[41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
##[43] "mouse.dpw.metadata"            "mouse.time.design"            
##[45] "stability.batch"               "stability.files"
[41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
[43] "mouse.dpw.metadata"            "mouse.time.design"            
[45] "stability.batch"               "stability.files"              

La fonction path sert à localiser le fichier, puis list.files va permettre d'ouvrir le fichier et vérifier la présence de nos données.

## Création des objets fnFs et fnRs, listes des fichiers Forward et Reverse
```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```
La fonction Sort permet un tri par ordre alphabétique, puis path, pattern="_R1_001.fastq" va retenir tous les fichiers contenant "_R1_001.fastq" dans le dossier "path"

Récupération des fichiers Forward(R1) et Reverse(R2), puis extraction des noms des échantillons.

## Visualisation de qualité des lectures
### Forward
```{r}
plotQualityProfile(fnFs[1:2])
```
![](unnamed-chunk-4-1.png)<!-- -->
Ceci est un graphique permettant de constater la qualité des nucléotides à travers les lectures.

Le dégradé de noir représente la fréquence de score pour les différentes positions de lecture

La ligne verte es la moyenne du score de qualité pour chacune des positions.

La ligne orange (du haut) correspond au quartiles.

La ligne rouge représente la proportion de lecture qui atteignent bel et bien les différenets positions. 

La lecture est ici de bonne qualité jusqu'à environs 240 paire de bases.

### Reverse
```{r}
plotQualityProfile(fnRs[1:2])
```
![](unnamed-chunk-6-1.png)<!-- -->

La qualité des lectures Reverse est bien plus basse que celle des lectures forward. En effet, une chute de cette qualité peut être observée à partir d'environ 155 paire de bases. Cela n'est pas étonnant pour des données Illumina, mais dada2 présente un avantage vis à vis de cela en intégrant des informations de bonne qualité dans son modèle d'erreur, compensant donc ces séquences de "mauvaise" qualité. 

## Nommer fichiers filtrés 
```{r}
filtFs <- file.path(path, "filtered", # Crée un chemin vers un sous-dosier "filtered" à l'intérieur du dossier 'path"
                    paste0(sample.names, "_F_filt.fastq.gz")) # Crée le nom du fichier filtré pour chaque échantillon

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names # Attribution de ces noms aux fichiers
```

Cette étape permet en premier lieu de créer un nouveau dossier "filtered" et de créer des noms à ces fichiers filtrés ; avant de leur attribuer. 

Nos fichiers sont donc maintenant reconnaissables et leur nom est lié avec celui de leur échantillon.

## Filtrage des séquences
```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), # Tronque les lectures Forward à 240 bases et Reverse à 160 bases
              maxN=0, # Elimine toutes lecture contenant une base inconnue
              maxEE=c(2,2), # Maximum d'erreurs par lecture
              truncQ=2, # Tronque la lecture dès que son score de qualité est inférieur à 2
              rm.phix=TRUE, # Supprime les lectures provenant du génome PhiX
              compress=TRUE, multithread=FALSE) 
head(out)
```
##                               reads.in reads.out
## F3D0_S188_L001_R1_001.fastq       7793      7113
## F3D1_S189_L001_R1_001.fastq       5869      5299
## F3D141_S207_L001_R1_001.fastq     5958      5463
## F3D142_S208_L001_R1_001.fastq     3183      2914
## F3D143_S209_L001_R1_001.fastq     3178      2941
## F3D144_S210_L001_R1_001.fastq     4827      4312

Cette ligne de code permet de nous montrer les lectures avant (reads.in) et après filtrage (reads.out). Cela nous permet de constater que la majorité des lectures ont été conservées malgrès le filtrage.

## Modèle d'erreur 
```{r}
errF <- learnErrors(filtFs, multithread=FALSE)
```
```{r}
errR <- learnErrors(filtRs, multithread=FALSE)
```
```{r}
plotErrors(errF, nominalQ=TRUE)
```
```{r}
errR <- learnErrors(filtRs, multithread=FALSE)
```


```{r}
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
```

```{r}
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
```
```{r}
dadaFs[[1]]
```
```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```
```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```
```{r}
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
```
```{r}
sum(seqtab.nochim)/sum(seqtab)
```

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "~/dossier_1/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
```

```{r}
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```
```{r}
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

```{r}
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```
