---
title: "DADA2 pipeline tutorial"
output: github_document
---

DADA2 
```{r, eval=FALSE}
#téléchargement du manager de Bioconductor pour ensuite télécharger dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.23")
```


```{r}
#ajout du package dada2 dans la bibliothèque
library(dada2); packageVersion("dada2")
```

```{bash, eval=FALSE}
#téléchargement et extraction des données fastq
wget "https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"
unzip "miseqsopdata.zip"
```

```{r}
#définit le chemin pour accéder aux données fastq, et les affiche
path <- "~/MiSeq_SOP" 
list.files(path)
```
```{r}
#liste pour séparer les reads forward et reverse (R1 et R2)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# extraction des noms des fichiers selon le texte avant le premier "_" du nom
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

```{r}
#visualisation des scores de qualité forward
plotQualityProfile(fnFs[1:2])
```

```{r}
#visualisation des scores de qualité reverse
plotQualityProfile(fnRs[1:2])
```

```{r}
#définit le chemin pour les fichiers filtrés dans /filtered
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#associe les fichiers aux échantillons
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

```{r}
#paramètres du filtrage selon les scores de qualité (coupe les reads F à 240 pb et les reads R à 160 pb)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)
head(out)
```


```{r}
#learnErrors fait apprendre les taux d'erreur, pour les visualiser avec plotErrors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #graphique des erreurs observées contre le modèle (ligne noire) et la théorie (ligne rouge)
```




```{r}
#utilisation de dada2 pour séparer les erreurs de séquençage des données
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]] #visualisation du premier échantillon
```


```{r}
#fusion des reads forward et reverse pour en faire une séquence complète
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) 
head(mergers[[1]]) #visualisation des séquences fusionnées
```

```{r}
#création de la table des ASV
seqtab <- makeSequenceTable(mergers) 
dim(seqtab) #nombre d'échantillons et d'ASV
table(nchar(getSequences(seqtab))) #distribution des longueurs des ASV
```

```{r}
#enlève les chimères
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #nombre d'ASV restants
```


```{r}
#proportion des reads conservés
sum(seqtab.nochim)/sum(seqtab)
```

```{r}
#fonction pour compter le nombre total de reads
getN <- function(x) sum(getUniques(x))
#création du tableau pour afficher le nombre de reads à chaque étape du pipeline dada2
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names 
head(track) #visualisation du tableau 
```


```{r}
#assignation taxonomique selon la base de données SILVA 
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE, tryRC=TRUE) #ajout de tryRC=TRUE pour tester le complément inverse
```

```{r}
#affichage des assignations taxonomiques
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```


```{r}
#récupère les ASV de Mock
unqs.mock <- seqtab.nochim["Mock",]
#ne garde que les ASV confirmés dans Mock
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```


```{r}
#compare les vrais ASV de Mock à ceux trouvés par dada2
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref)))) #compare les ASV à la référence
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```



```{r, eval=FALSE}
#installation de phyloseq et igraph
install.packages("igraph", type = "source")
BiocManager::install("phyloseq", update = FALSE)
```

```{r}
#chargement des packages nécessaires
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
```

```{r}
#création du tableau des métadonnées
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

```{r}
#regroupe les ASV avec les métadonnées et la taxonomie
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))
#enlève l’échantillon Mock
ps <- prune_samples(sample_names(ps) != "Mock", ps) 
```

```{r}
#garde les données d'ADN mais raccourcit les noms en ASV_
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

```{r}
# graphique de la diversité alpha
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

```{r}
# Transforme les données pour être appropriées au calcul de distance Bray-Curtis
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

```{r}
# graphique de la diversité beta
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

```{r}
#barplot des ASV les plus communs selon les familles
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```
