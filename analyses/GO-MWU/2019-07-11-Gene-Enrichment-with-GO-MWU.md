Gene Enrichment with GO-MWU
================
Yaamini Venkataraman
7/11/2019

In this R Markdown file, we format files with a gene enrichment with [GO-MWU](https://github.com/z0on/GO_MWU). This is a rank-based gene enrichment method using a Mann-Whitney U test that works well with non-model organisms. There are separate annotations for *Zostera marina* and *Labyrinthula zosterae*, so these files will be formatted separately.

Set up R Markdown document
==========================

``` r
sessionInfo() #Obtain session information.
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.14.5
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.5.0  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
    ##  [5] tools_3.5.0     htmltools_0.3.6 yaml_2.2.0      Rcpp_1.0.0     
    ##  [9] stringi_1.2.4   rmarkdown_1.10  knitr_1.20      stringr_1.3.1  
    ## [13] digest_0.6.15   evaluate_0.11

Format GO-MWU input
===================

GO annotations table
--------------------

The GO annotation tables need to be tab-delimited with 2 columns: gene ID and GOterms separated by semicolons.

### *Z. marina*

``` r
allZosteraAnnotations <- read.delim("../../data/Zostera-blast-annot-withGeneID-noIsoforms-geneCounts.tab", header = FALSE) #Import file with gene counts for annotated genes
colnames(allZosteraAnnotations) <- c("GeneID", "Accession", "Isoform", "E-value", "ProteinN", "GO_BP", "GO_CC", "GO_MF", "GO", "Status", "Organism", "S_10B", "S_9A", "S_13A", "S_42A", "S_46B", "S_47B", "S_48B", "S_2A", "S_2B", "S_7B", "S_8B", "S_33A", "S_36B", "S_38A", "S_40A") #Rename columns
head(allZosteraAnnotations) #Confirm changes
```

    ##                   GeneID Accession                   Isoform  E-value
    ## 1 TRINITY_DN102005_c0_g1    Q1DZQ0 TRINITY_DN102005_c0_g1_i1  1.3e-37
    ## 2 TRINITY_DN102006_c0_g1    Q54NS9 TRINITY_DN102006_c0_g1_i1  3.0e-21
    ## 3 TRINITY_DN102012_c0_g1    Q56872 TRINITY_DN102012_c0_g1_i1 1.9e-125
    ## 4 TRINITY_DN102016_c0_g1    F4I460 TRINITY_DN102016_c0_g1_i1  4.6e-15
    ## 5 TRINITY_DN102024_c0_g1    Q9LK31 TRINITY_DN102024_c0_g1_i1  1.9e-34
    ## 6 TRINITY_DN102050_c0_g1    Q6DIF4 TRINITY_DN102050_c0_g1_i1  3.4e-32
    ##                                                                ProteinN
    ## 1                                       Protein transport protein SEC13
    ## 2                      Apoptosis-inducing factor homolog A (EC 1.-.-.-)
    ## 3 GDP-mannose 4,6-dehydratase (EC 4.2.1.47) (GDP-D-mannose dehydratase)
    ## 4                                        Myosin-8 (Myosin XI B) (AtXIB)
    ## 5                             Kelch repeat-containing protein At3g27220
    ## 6   WD repeat-containing protein 1 (Actin-interacting protein 1) (AIP1)
    ##                                                                                                                                                                                                                        GO_BP
    ## 1                                                                COPII-coated vesicle budding [GO:0090114]; mRNA transport [GO:0051028]; positive regulation of TORC1 signaling [GO:1904263]; protein transport [GO:0015031]
    ## 2                                                                                                                                                                                                                           
    ## 3                                                                          'de novo' GDP-L-fucose biosynthetic process [GO:0042351]; GDP-mannose metabolic process [GO:0019673]; O antigen biosynthetic process [GO:0009243]
    ## 4                                                                                                    actin filament-based movement [GO:0030048]; actin filament organization [GO:0007015]; root hair elongation [GO:0048767]
    ## 5                                                                                                                                                                                         anaerobic respiration [GO:0009061]
    ## 6 actin filament depolymerization [GO:0030042]; actin filament fragmentation [GO:0030043]; locomotion [GO:0040011]; positive regulation of actin filament depolymerization [GO:0030836]; sarcomere organization [GO:0045214]
    ##                                                                                                                                                                                       GO_CC
    ## 1                                                      COPII vesicle coat [GO:0030127]; endoplasmic reticulum membrane [GO:0005789]; Golgi membrane [GO:0000139]; nuclear pore [GO:0005643]
    ## 2                                                                                                                                        cytoplasm [GO:0005737]; lipid droplet [GO:0005811]
    ## 3                                                                                                                                                                                          
    ## 4                                                                                                                                       cytoplasm [GO:0005737]; myosin complex [GO:0016459]
    ## 5                                                        endosome [GO:0005768]; Golgi apparatus [GO:0005794]; integral component of membrane [GO:0016021]; trans-Golgi network [GO:0005802]
    ## 6 actin filament [GO:0005884]; actomyosin, actin portion [GO:0042643]; cortical actin cytoskeleton [GO:0030864]; cytoplasm [GO:0005737]; nucleus [GO:0005634]; plasma membrane [GO:0005886]
    ##                                                                                                                         GO_MF
    ## 1                                                                                   structural molecule activity [GO:0005198]
    ## 2    electron-transferring-flavoprotein dehydrogenase activity [GO:0004174]; flavin adenine dinucleotide binding [GO:0050660]
    ## 3                                               GDP-mannose 4,6-dehydratase activity [GO:0008446]; NADP+ binding [GO:0070401]
    ## 4 actin filament binding [GO:0051015]; ATP binding [GO:0005524]; calmodulin binding [GO:0005516]; motor activity [GO:0003774]
    ## 5                                                                                                                            
    ## 6                                                             actin binding [GO:0003779]; actin filament binding [GO:0051015]
    ##                                                                                                                                                           GO
    ## 1                                                 GO:0000139; GO:0005198; GO:0005643; GO:0005789; GO:0015031; GO:0030127; GO:0051028; GO:0090114; GO:1904263
    ## 2                                                                                                             GO:0004174; GO:0005737; GO:0005811; GO:0050660
    ## 3                                                                                                 GO:0008446; GO:0009243; GO:0019673; GO:0042351; GO:0070401
    ## 4                                                 GO:0003774; GO:0005516; GO:0005524; GO:0005737; GO:0007015; GO:0016459; GO:0030048; GO:0048767; GO:0051015
    ## 5                                                                                                 GO:0005768; GO:0005794; GO:0005802; GO:0009061; GO:0016021
    ## 6 GO:0003779; GO:0005634; GO:0005737; GO:0005884; GO:0005886; GO:0030042; GO:0030043; GO:0030836; GO:0030864; GO:0040011; GO:0042643; GO:0045214; GO:0051015
    ##     Status
    ## 1 reviewed
    ## 2 reviewed
    ## 3 reviewed
    ## 4 reviewed
    ## 5 reviewed
    ## 6 reviewed
    ##                                                                       Organism
    ## 1                       Coccidioides immitis (strain RS) (Valley fever fungus)
    ## 2                                        Dictyostelium discoideum (Slime mold)
    ## 3 Yersinia enterocolitica serotype O:8 / biotype 1B (strain NCTC 13174 / 8081)
    ## 4                                       Arabidopsis thaliana (Mouse-ear cress)
    ## 5                                       Arabidopsis thaliana (Mouse-ear cress)
    ## 6               Xenopus tropicalis (Western clawed frog) (Silurana tropicalis)
    ##   S_10B S_9A S_13A S_42A S_46B S_47B S_48B S_2A S_2B S_7B S_8B S_33A S_36B
    ## 1     0    0     1     1     0     0     0    0    0    4    0     1     3
    ## 2     0    1     0     0     0     1     0    3    1   44    0     0     4
    ## 3     0    1     0     8     2     0     1    1    0    8    0     0     3
    ## 4     0    2     0     0     0     0     0    0    0    0    0     0     0
    ## 5     3    0     0     0     0     0     2    3    0    0    2    18     0
    ## 6     2    4     0     0     0     0     0   11    2    5    0     0    10
    ##   S_38A S_40A
    ## 1     0     0
    ## 2     1     2
    ## 3     0    19
    ## 4     0     0
    ## 5     0     0
    ## 6     0     0

``` r
zosteraGOAnnotationsTable <- data.frame(allZosteraAnnotations$GeneID,
                                        allZosteraAnnotations$GO) #Save gene ID and all gene ontology terms in a new dataframe
head(zosteraGOAnnotationsTable) #Confirm format is good for GO-MWU. It is.
```

    ##   allZosteraAnnotations.GeneID
    ## 1       TRINITY_DN102005_c0_g1
    ## 2       TRINITY_DN102006_c0_g1
    ## 3       TRINITY_DN102012_c0_g1
    ## 4       TRINITY_DN102016_c0_g1
    ## 5       TRINITY_DN102024_c0_g1
    ## 6       TRINITY_DN102050_c0_g1
    ##                                                                                                                                     allZosteraAnnotations.GO
    ## 1                                                 GO:0000139; GO:0005198; GO:0005643; GO:0005789; GO:0015031; GO:0030127; GO:0051028; GO:0090114; GO:1904263
    ## 2                                                                                                             GO:0004174; GO:0005737; GO:0005811; GO:0050660
    ## 3                                                                                                 GO:0008446; GO:0009243; GO:0019673; GO:0042351; GO:0070401
    ## 4                                                 GO:0003774; GO:0005516; GO:0005524; GO:0005737; GO:0007015; GO:0016459; GO:0030048; GO:0048767; GO:0051015
    ## 5                                                                                                 GO:0005768; GO:0005794; GO:0005802; GO:0009061; GO:0016021
    ## 6 GO:0003779; GO:0005634; GO:0005737; GO:0005884; GO:0005886; GO:0030042; GO:0030043; GO:0030836; GO:0030864; GO:0040011; GO:0042643; GO:0045214; GO:0051015

``` r
write.table(zosteraGOAnnotationsTable, "2019-07-11-Zostera-GO-Annotations-Table.tab", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) #Save file
```

``` bash
head 2019-07-11-Zostera-GO-Annotations-Table.tab #Confirm formatting is good.
```

    ## TRINITY_DN102005_c0_g1   GO:0000139; GO:0005198; GO:0005643; GO:0005789; GO:0015031; GO:0030127; GO:0051028; GO:0090114; GO:1904263
    ## TRINITY_DN102006_c0_g1   GO:0004174; GO:0005737; GO:0005811; GO:0050660
    ## TRINITY_DN102012_c0_g1   GO:0008446; GO:0009243; GO:0019673; GO:0042351; GO:0070401
    ## TRINITY_DN102016_c0_g1   GO:0003774; GO:0005516; GO:0005524; GO:0005737; GO:0007015; GO:0016459; GO:0030048; GO:0048767; GO:0051015
    ## TRINITY_DN102024_c0_g1   GO:0005768; GO:0005794; GO:0005802; GO:0009061; GO:0016021
    ## TRINITY_DN102050_c0_g1   GO:0003779; GO:0005634; GO:0005737; GO:0005884; GO:0005886; GO:0030042; GO:0030043; GO:0030836; GO:0030864; GO:0040011; GO:0042643; GO:0045214; GO:0051015
    ## TRINITY_DN102051_c0_g1   GO:0005524; GO:0005886; GO:0016021; GO:0016887; GO:0033344; GO:0042626; GO:0042803; GO:0046982; GO:0055085; GO:1990830
    ## TRINITY_DN102060_c0_g1   GO:0004867; GO:0030212; GO:0062023
    ## TRINITY_DN102062_c0_g1   GO:0005739; GO:0016787; GO:0018773; GO:0046872
    ## TRINITY_DN10208_c0_g1    GO:0003723; GO:0003729; GO:0005737; GO:0005739; GO:0006378; GO:0008033; GO:0009952; GO:0016331; GO:0035770; GO:0042645

### *L. zosterae*

``` r
allnonZosteraAnnotations <- read.delim("../../data/nonZostera-blast-annot-withGeneID-noIsoforms-geneCounts.tab", header = FALSE) #Import file with gene counts for annotated genes
colnames(allnonZosteraAnnotations) <- c("GeneID", "Accession", "Isoform", "E-value", "ProteinN", "GO_BP", "GO_CC", "GO_MF", "GO", "Status", "Organism", "S_10B", "S_9A", "S_13A", "S_42A", "S_46B", "S_47B", "S_48B", "S_2A", "S_2B", "S_7B", "S_8B", "S_33A", "S_36B", "S_38A", "S_40A") #Rename columns
head(allnonZosteraAnnotations) #Confirm changes
```

    ##                   GeneID Accession                   Isoform E-value
    ## 1 TRINITY_DN102004_c0_g1    Q54BM8 TRINITY_DN102004_c0_g1_i1 6.5e-19
    ## 2 TRINITY_DN102011_c0_g1    Q84M24 TRINITY_DN102011_c0_g1_i1 2.5e-14
    ## 3 TRINITY_DN102014_c0_g1    A5DB51 TRINITY_DN102014_c0_g1_i1 2.7e-09
    ## 4 TRINITY_DN102016_c0_g1    F4I460 TRINITY_DN102016_c0_g1_i2 8.7e-06
    ## 5  TRINITY_DN10203_c0_g1    Q19020  TRINITY_DN10203_c0_g1_i1 2.5e-15
    ## 6 TRINITY_DN102051_c0_g1    Q9H172 TRINITY_DN102051_c0_g1_i2 1.5e-06
    ##                                                                                                    ProteinN
    ## 1                                                                                           UPF0652 protein
    ## 2 ABC transporter A family member 1 (ABC transporter ABCA.1) (AtABCA1) (ABC one homolog protein 1) (AtAOH1)
    ## 3                                                                     Riboflavin synthase (RS) (EC 2.5.1.9)
    ## 4                                                                            Myosin-8 (Myosin XI B) (AtXIB)
    ## 5                                    Dynein heavy chain, cytoplasmic (Dynein heavy chain, cytosolic) (DYHC)
    ## 6                                                                ATP-binding cassette sub-family G member 4
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 GO_BP
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        lipid transport [GO:0006869]
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        riboflavin biosynthetic process [GO:0009231]
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             actin filament-based movement [GO:0030048]; actin filament organization [GO:0007015]; root hair elongation [GO:0048767]
    ## 5 anterograde neuronal dense core vesicle transport [GO:1990048]; cytoplasmic microtubule organization [GO:0031122]; establishment of meiotic spindle orientation [GO:0051296]; establishment of spindle localization [GO:0051293]; maintenance of centrosome location [GO:0051661]; microtubule-based movement [GO:0007018]; minus-end-directed vesicle transport along microtubule [GO:0072382]; mitotic cell cycle [GO:0000278]; neuron remodeling [GO:0016322]; nuclear migration [GO:0007097]; positive regulation of dense core granule transport [GO:1904811]; regulation of dendrite morphogenesis [GO:0048814]; regulation of protein localization to synapse [GO:1902473]; retrograde axonal transport [GO:0008090]; synaptic transmission, GABAergic [GO:0051932]; synaptic vesicle transport [GO:0048489]
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 cellular response to leukemia inhibitory factor [GO:1990830]; cholesterol efflux [GO:0033344]; transmembrane transport [GO:0055085]
    ##                                                                                                                                                                                                                                                                                                 GO_CC
    ## 1                                                                                                                                                                                                                                                                                                    
    ## 2                                                                                                                                                                  integral component of membrane [GO:0016021]; intracellular membrane-bounded organelle [GO:0043231]; vacuolar membrane [GO:0005774]
    ## 3                                                                                                                                                                                                                                                                                                    
    ## 4                                                                                                                                                                                                                                                 cytoplasm [GO:0005737]; myosin complex [GO:0016459]
    ## 5 axon cytoplasm [GO:1904115]; cell cortex [GO:0005938]; cytoplasm [GO:0005737]; cytoplasmic dynein complex [GO:0005868]; cytoplasmic microtubule [GO:0005881]; dynein complex [GO:0030286]; kinetochore [GO:0000776]; nuclear envelope [GO:0005635]; spindle [GO:0005819]; spindle pole [GO:0000922]
    ## 6                                                                                                                                                                                                                           integral component of membrane [GO:0016021]; plasma membrane [GO:0005886]
    ##                                                                                                                                                                                                                                                GO_MF
    ## 1                                                                                                                                                                                                                      zinc ion binding [GO:0008270]
    ## 2                                                                                    ATPase activity [GO:0016887]; ATPase-coupled transmembrane transporter activity [GO:0042626]; ATP binding [GO:0005524]; lipid transporter activity [GO:0005319]
    ## 3                                                                                                                                                                                                          riboflavin synthase activity [GO:0004746]
    ## 4                                                                                                                        actin filament binding [GO:0051015]; ATP binding [GO:0005524]; calmodulin binding [GO:0005516]; motor activity [GO:0003774]
    ## 5 ATP binding [GO:0005524]; ATP-dependent microtubule motor activity, minus-end-directed [GO:0008569]; dynein intermediate chain binding [GO:0045505]; dynein light chain binding [GO:0045503]; dynein light intermediate chain binding [GO:0051959]
    ## 6                           ATPase activity [GO:0016887]; ATPase-coupled transmembrane transporter activity [GO:0042626]; ATP binding [GO:0005524]; protein heterodimerization activity [GO:0046982]; protein homodimerization activity [GO:0042803]
    ##                                                                                                                                                                                                                                                                                                                                                                                   GO
    ## 1                                                                                                                                                                                                                                                                                                                                                                         GO:0008270
    ## 2                                                                                                                                                                                                                                                                                     GO:0005319; GO:0005524; GO:0005774; GO:0006869; GO:0016021; GO:0016887; GO:0042626; GO:0043231
    ## 3                                                                                                                                                                                                                                                                                                                                                             GO:0004746; GO:0009231
    ## 4                                                                                                                                                                                                                                                                         GO:0003774; GO:0005516; GO:0005524; GO:0005737; GO:0007015; GO:0016459; GO:0030048; GO:0048767; GO:0051015
    ## 5 GO:0000278; GO:0000776; GO:0000922; GO:0005524; GO:0005635; GO:0005737; GO:0005819; GO:0005868; GO:0005881; GO:0005938; GO:0007018; GO:0007097; GO:0008090; GO:0008569; GO:0016322; GO:0030286; GO:0031122; GO:0045503; GO:0045505; GO:0048489; GO:0048814; GO:0051293; GO:0051296; GO:0051661; GO:0051932; GO:0051959; GO:0072382; GO:1902473; GO:1904115; GO:1904811; GO:1990048
    ## 6                                                                                                                                                                                                                                                             GO:0005524; GO:0005886; GO:0016021; GO:0016887; GO:0033344; GO:0042626; GO:0042803; GO:0046982; GO:0055085; GO:1990830
    ##     Status
    ## 1 reviewed
    ## 2 reviewed
    ## 3 reviewed
    ## 4 reviewed
    ## 5 reviewed
    ## 6 reviewed
    ##                                                                                                                                  Organism
    ## 1                                                                                                   Dictyostelium discoideum (Slime mold)
    ## 2                                                                                                  Arabidopsis thaliana (Mouse-ear cress)
    ## 3 Meyerozyma guilliermondii (strain ATCC 6260 / CBS 566 / DSM 6381 / JCM 1539 / NBRC 10279 / NRRL Y-324) (Yeast) (Candida guilliermondii)
    ## 4                                                                                                  Arabidopsis thaliana (Mouse-ear cress)
    ## 5                                                                                                                  Caenorhabditis elegans
    ## 6                                                                                                                    Homo sapiens (Human)
    ##   S_10B S_9A S_13A S_42A S_46B S_47B S_48B S_2A S_2B S_7B S_8B S_33A S_36B
    ## 1     0    1     0     0     0     0     0    0    0    4    0     0     6
    ## 2     0   11     0     0     0     0     0    0    0    0    0     0     0
    ## 3     0    1     0     0     0     0     0    1    0    3    0     0     0
    ## 4     0    2     0     0     0     0     0    0    0    0    0     0     0
    ## 5     0    0     1     0     0     2     0    0    0    0    0     0     0
    ## 6     0    0     0     1     3     0     0    4    0    0    2     0     1
    ##   S_38A S_40A
    ## 1     0     0
    ## 2     0     0
    ## 3     0     0
    ## 4     0     0
    ## 5     0     0
    ## 6     0     0

``` r
nonZosteraGOAnnotationsTable <- data.frame(allnonZosteraAnnotations$GeneID,
                                           allnonZosteraAnnotations$GO) #Save gene ID and all gene ontology terms in a new dataframe
head(nonZosteraGOAnnotationsTable) #Confirm format is good for GO-MWU. It is.
```

    ##   allnonZosteraAnnotations.GeneID
    ## 1          TRINITY_DN102004_c0_g1
    ## 2          TRINITY_DN102011_c0_g1
    ## 3          TRINITY_DN102014_c0_g1
    ## 4          TRINITY_DN102016_c0_g1
    ## 5           TRINITY_DN10203_c0_g1
    ## 6          TRINITY_DN102051_c0_g1
    ##                                                                                                                                                                                                                                                                                                                                                          allnonZosteraAnnotations.GO
    ## 1                                                                                                                                                                                                                                                                                                                                                                         GO:0008270
    ## 2                                                                                                                                                                                                                                                                                     GO:0005319; GO:0005524; GO:0005774; GO:0006869; GO:0016021; GO:0016887; GO:0042626; GO:0043231
    ## 3                                                                                                                                                                                                                                                                                                                                                             GO:0004746; GO:0009231
    ## 4                                                                                                                                                                                                                                                                         GO:0003774; GO:0005516; GO:0005524; GO:0005737; GO:0007015; GO:0016459; GO:0030048; GO:0048767; GO:0051015
    ## 5 GO:0000278; GO:0000776; GO:0000922; GO:0005524; GO:0005635; GO:0005737; GO:0005819; GO:0005868; GO:0005881; GO:0005938; GO:0007018; GO:0007097; GO:0008090; GO:0008569; GO:0016322; GO:0030286; GO:0031122; GO:0045503; GO:0045505; GO:0048489; GO:0048814; GO:0051293; GO:0051296; GO:0051661; GO:0051932; GO:0051959; GO:0072382; GO:1902473; GO:1904115; GO:1904811; GO:1990048
    ## 6                                                                                                                                                                                                                                                             GO:0005524; GO:0005886; GO:0016021; GO:0016887; GO:0033344; GO:0042626; GO:0042803; GO:0046982; GO:0055085; GO:1990830

``` r
write.table(nonZosteraGOAnnotationsTable, "2019-07-11-nonZostera-GO-Annotations-Table.tab", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) #Save file
```

``` bash
head 2019-07-11-nonZostera-GO-Annotations-Table.tab
```

    ## TRINITY_DN102004_c0_g1   GO:0008270
    ## TRINITY_DN102011_c0_g1   GO:0005319; GO:0005524; GO:0005774; GO:0006869; GO:0016021; GO:0016887; GO:0042626; GO:0043231
    ## TRINITY_DN102014_c0_g1   GO:0004746; GO:0009231
    ## TRINITY_DN102016_c0_g1   GO:0003774; GO:0005516; GO:0005524; GO:0005737; GO:0007015; GO:0016459; GO:0030048; GO:0048767; GO:0051015
    ## TRINITY_DN10203_c0_g1    GO:0000278; GO:0000776; GO:0000922; GO:0005524; GO:0005635; GO:0005737; GO:0005819; GO:0005868; GO:0005881; GO:0005938; GO:0007018; GO:0007097; GO:0008090; GO:0008569; GO:0016322; GO:0030286; GO:0031122; GO:0045503; GO:0045505; GO:0048489; GO:0048814; GO:0051293; GO:0051296; GO:0051661; GO:0051932; GO:0051959; GO:0072382; GO:1902473; GO:1904115; GO:1904811; GO:1990048
    ## TRINITY_DN102051_c0_g1   GO:0005524; GO:0005886; GO:0016021; GO:0016887; GO:0033344; GO:0042626; GO:0042803; GO:0046982; GO:0055085; GO:1990830
    ## TRINITY_DN102053_c0_g1   GO:0004177; GO:0004301; GO:0004463; GO:0005634; GO:0005737; GO:0006691; GO:0008237; GO:0008270; GO:0010043; GO:0019370; GO:0043171; GO:0043434; GO:0044255; GO:0060509
    ## TRINITY_DN10209_c0_g1    GO:0005886; GO:0015297; GO:0016021; GO:0019531
    ## TRINITY_DN102366_c0_g1   GO:0005737; GO:0046872; GO:0070006
    ## TRINITY_DN10241_c0_g1    GO:0008270

Table of significance measures
------------------------------

The table of significance measures needs to be a CSV file with two columns: gene ID and measure of continuous significance. All genes need to be included, not just those deemed significantly different by `edgeR`. Similar to the GO-MWU README, I used signed negative log p-values modified from `edgeR` output.

### *Z. marina*

#### Subset table

``` r
zosteraDEG <- read.delim("../EdgeR/EXP.CON.Zostera.ALL.txt", header = TRUE) #Import all genes used in edgeR output
head(zosteraDEG) #Confirm import
```

    ##                           logFC    logCPM        F       PValue
    ## TRINITY_DN312093_c0_g7 22.11531 10.993989 502.8797 3.799002e-17
    ## TRINITY_DN312726_c1_g3 21.91717 12.272552 406.6525 3.911955e-16
    ## TRINITY_DN311783_c0_g2 21.66304  9.696917 384.3782 7.231585e-16
    ## TRINITY_DN295224_c1_g2 23.52990 10.650603 381.3848 7.874138e-16
    ## TRINITY_DN316951_c1_g2 21.70589  9.761771 368.7151 1.137054e-15
    ## TRINITY_DN276264_c0_g3 21.59372 11.329075 357.0113 1.613840e-15
    ##                                 FDR
    ## TRINITY_DN312093_c0_g7 1.175791e-13
    ## TRINITY_DN312726_c1_g3 6.053751e-13
    ## TRINITY_DN311783_c0_g2 6.092614e-13
    ## TRINITY_DN295224_c1_g2 6.092614e-13
    ## TRINITY_DN316951_c1_g2 7.038364e-13
    ## TRINITY_DN276264_c0_g3 8.151882e-13

``` r
zosteraSigMeasures <- data.frame("gene" = row.names(zosteraDEG),
                                 "logFC" = zosteraDEG$logFC,
                                 "logP" = log(zosteraDEG$PValue)) #Subset row names as gene ID and log fold change. Also include log(p-value)
head(zosteraSigMeasures) #Confirm import
```

    ##                     gene    logFC      logP
    ## 1 TRINITY_DN312093_c0_g7 22.11531 -37.80921
    ## 2 TRINITY_DN312726_c1_g3 21.91717 -35.47732
    ## 3 TRINITY_DN311783_c0_g2 21.66304 -34.86290
    ## 4 TRINITY_DN295224_c1_g2 23.52990 -34.77778
    ## 5 TRINITY_DN316951_c1_g2 21.70589 -34.41034
    ## 6 TRINITY_DN276264_c0_g3 21.59372 -34.06016

#### Correct p-values for upregulated genes

``` r
zosteraSigMeasuresUpregulated <- subset(x = zosteraSigMeasures, subset = zosteraSigMeasures$logFC > 0) #Subset differentially expressed genes that are upregulated by pulling out genes with a positive fold change
tail(zosteraSigMeasuresUpregulated) #Confirm subset
```

    ##                        gene     logFC          logP
    ## 3070 TRINITY_DN274713_c0_g1 4.5350773 -0.0003038720
    ## 3073   TRINITY_DN1723_c0_g1 1.7289533 -0.0002990235
    ## 3075 TRINITY_DN299782_c2_g2 2.5613345 -0.0002967727
    ## 3078 TRINITY_DN294772_c3_g1 1.0233543 -0.0002903804
    ## 3084   TRINITY_DN2516_c0_g1 0.4424923 -0.0002836004
    ## 3095 TRINITY_DN292910_c0_g4 0.0224682 -0.0001045594

``` r
range(zosteraSigMeasuresUpregulated$logP) #Look at the range of p-values for upregulated genes. All values are negative. For signed negative log p-values, the sign of the log p-value should correspond with upregulation for downregulation. All upregulated genes should have positive p-values.
```

    ## [1] -3.780921e+01 -1.045594e-04

``` r
zosteraSigMeasuresUpregulated$signLogP <- -1 * zosteraSigMeasuresUpregulated$logP #Multiply all p-values by -1 so they become positive
head(zosteraSigMeasuresUpregulated) #Confirm changes
```

    ##                     gene    logFC      logP signLogP
    ## 1 TRINITY_DN312093_c0_g7 22.11531 -37.80921 37.80921
    ## 2 TRINITY_DN312726_c1_g3 21.91717 -35.47732 35.47732
    ## 3 TRINITY_DN311783_c0_g2 21.66304 -34.86290 34.86290
    ## 4 TRINITY_DN295224_c1_g2 23.52990 -34.77778 34.77778
    ## 5 TRINITY_DN316951_c1_g2 21.70589 -34.41034 34.41034
    ## 6 TRINITY_DN276264_c0_g3 21.59372 -34.06016 34.06016

#### Correct p-values for downregulated genes

``` r
zosteraSigMeasuresDownregulated <- subset(x = zosteraSigMeasures, subset = zosteraSigMeasures$logFC < 0) #Subset differentially expressed genes that are downregulated by pulling out genes with a negative fold change
head(zosteraSigMeasuresDownregulated) #Confirm subset
```

    ##                       gene      logFC      logP
    ## 108 TRINITY_DN292378_c0_g2 -12.699776 -20.99580
    ## 112 TRINITY_DN293499_c0_g1 -11.744522 -20.38603
    ## 113 TRINITY_DN295071_c1_g1  -7.504449 -20.21674
    ## 121 TRINITY_DN311070_c0_g3  -9.955325 -17.56470
    ## 122 TRINITY_DN312001_c0_g1  -7.636863 -17.44754
    ## 123 TRINITY_DN298476_c0_g1  -6.179492 -17.36296

``` r
range(zosteraSigMeasuresDownregulated$logP) #All downregulated genes should have negative p-values. These p-values are already negative, so no more modification needs to happen.
```

    ## [1] -2.099580e+01 -2.111773e-04

``` r
zosteraSigMeasuresDownregulated$signLogP <- zosteraSigMeasuresDownregulated$logP #Copy entries to a new column to be consistent with formatting for table with upregulated genes
```

#### Combine tables with corrected p-values

``` r
zosteraSigMeasuresCorrected <- rbind(zosteraSigMeasuresUpregulated, zosteraSigMeasuresDownregulated) #Combine tables by columns
zosteraSigMeasuresCorrected <- zosteraSigMeasuresCorrected[,-c(2:3)] #Remove unnecessary columns
head(zosteraSigMeasuresCorrected) #Confirm changes
```

    ##                     gene signLogP
    ## 1 TRINITY_DN312093_c0_g7 37.80921
    ## 2 TRINITY_DN312726_c1_g3 35.47732
    ## 3 TRINITY_DN311783_c0_g2 34.86290
    ## 4 TRINITY_DN295224_c1_g2 34.77778
    ## 5 TRINITY_DN316951_c1_g2 34.41034
    ## 6 TRINITY_DN276264_c0_g3 34.06016

``` r
write.csv(zosteraSigMeasuresCorrected, "2019-07-11-Zostera-Table-of-Significance-Measures.csv", quote = FALSE, row.names = FALSE) #Save file
```

``` bash
head 2019-07-11-Zostera-Table-of-Significance-Measures.csv #Check that row name column was not included in output
```

    ## gene,signLogP
    ## TRINITY_DN312093_c0_g7,37.8092082631101
    ## TRINITY_DN312726_c1_g3,35.477324117117
    ## TRINITY_DN311783_c0_g2,34.8629032073826
    ## TRINITY_DN295224_c1_g2,34.7777778071615
    ## TRINITY_DN316951_c1_g2,34.4103356722594
    ## TRINITY_DN276264_c0_g3,34.0601601687074
    ## TRINITY_DN255876_c0_g1,33.8660382291265
    ## TRINITY_DN277915_c0_g1,33.793459082109
    ## TRINITY_DN312752_c3_g1,33.4013612718826

### *L. zosterae*

#### Subset table

``` r
nonzosteraDEG <- read.delim("../EdgeR/EXP.CON.NonZostera.all.txt", header = TRUE) #Import all genes used in edgeR output
head(nonzosteraDEG) #Confirm import
```

    ##                           logFC   logCPM        F       PValue
    ## TRINITY_DN296708_c0_g1 21.68291 14.66492 140.7744 1.089995e-13
    ## TRINITY_DN312737_c1_g1 22.74768 13.85166 127.8846 1.391260e-12
    ## TRINITY_DN294773_c0_g2 19.48342 14.07991 112.5191 2.278414e-12
    ## TRINITY_DN293894_c0_g2 18.76497 13.89717 111.8221 2.473772e-12
    ## TRINITY_DN317297_c2_g1 21.33990 12.54172 122.0424 2.508026e-12
    ## TRINITY_DN233011_c0_g1 19.06413 14.00250 108.5450 3.661602e-12
    ##                                 FDR
    ## TRINITY_DN296708_c0_g1 1.396284e-10
    ## TRINITY_DN312737_c1_g1 6.017860e-10
    ## TRINITY_DN294773_c0_g2 6.017860e-10
    ## TRINITY_DN293894_c0_g2 6.017860e-10
    ## TRINITY_DN317297_c2_g1 6.017860e-10
    ## TRINITY_DN233011_c0_g1 6.017860e-10

``` r
nonzosteraSigMeasures <- data.frame("gene" = row.names(nonzosteraDEG),
                                    "logFC" = nonzosteraDEG$logFC,
                                    "logP" = log(nonzosteraDEG$PValue)) #Subset row names as gene ID and log fold change. Also include log(p-value)
head(nonzosteraSigMeasures) #Confirm changes
```

    ##                     gene    logFC      logP
    ## 1 TRINITY_DN296708_c0_g1 21.68291 -29.84743
    ## 2 TRINITY_DN312737_c1_g1 22.74768 -27.30081
    ## 3 TRINITY_DN294773_c0_g2 19.48342 -26.80754
    ## 4 TRINITY_DN293894_c0_g2 18.76497 -26.72528
    ## 5 TRINITY_DN317297_c2_g1 21.33990 -26.71153
    ## 6 TRINITY_DN233011_c0_g1 19.06413 -26.33312

#### Correct p-values for upregulated genes

``` r
nonzosteraSigMeasuresUpregulated <- subset(x = nonzosteraSigMeasures, subset = nonzosteraSigMeasures$logFC > 0) #Subset differentially expressed genes that are upregulated by pulling out genes with a positive fold change
tail(nonzosteraSigMeasuresUpregulated) #Confirm subset
```

    ##                        gene     logFC          logP
    ## 1252 TRINITY_DN293683_c0_g2 2.1907131 -0.0003349213
    ## 1259 TRINITY_DN227741_c0_g1 4.0761051 -0.0003225405
    ## 1266 TRINITY_DN231560_c0_g1 2.0066553 -0.0003208705
    ## 1270    TRINITY_DN341_c0_g1 1.4630195 -0.0003128309
    ## 1277 TRINITY_DN222992_c0_g1 0.5253773 -0.0002930425
    ## 1280 TRINITY_DN344645_c0_g1 1.2791514 -0.0002879332

``` r
range(nonzosteraSigMeasuresUpregulated$logP) #Look at the range of p-values for upregulated genes. All values are negative. For signed negative log p-values, the sign of the log p-value should correspond with upregulation for downregulation. All upregulated genes should have positive p-values.
```

    ## [1] -2.984743e+01 -2.879332e-04

``` r
nonzosteraSigMeasuresUpregulated$signLogP <- -1 * nonzosteraSigMeasuresUpregulated$logP #Multiply all p-values by -1 so they become positive
head(nonzosteraSigMeasuresUpregulated) #Confirm changes
```

    ##                     gene    logFC      logP signLogP
    ## 1 TRINITY_DN296708_c0_g1 21.68291 -29.84743 29.84743
    ## 2 TRINITY_DN312737_c1_g1 22.74768 -27.30081 27.30081
    ## 3 TRINITY_DN294773_c0_g2 19.48342 -26.80754 26.80754
    ## 4 TRINITY_DN293894_c0_g2 18.76497 -26.72528 26.72528
    ## 5 TRINITY_DN317297_c2_g1 21.33990 -26.71153 26.71153
    ## 6 TRINITY_DN233011_c0_g1 19.06413 -26.33312 26.33312

#### Correct p-values for downregulated genes

``` r
nonzosteraSigMeasuresDownregulated <- subset(x = nonzosteraSigMeasures, subset = nonzosteraSigMeasures$logFC < 0) #Subset differentially expressed genes that are downregulated by pulling out genes with a negative fold change
head(nonzosteraSigMeasuresDownregulated) #Confirm subset
```

    ##                      gene      logFC       logP
    ## 52 TRINITY_DN310761_c3_g3  -7.856980 -12.227691
    ## 53 TRINITY_DN310761_c3_g1  -8.549530 -12.051182
    ## 56 TRINITY_DN294779_c0_g2  -7.295120 -11.616013
    ## 58 TRINITY_DN314888_c2_g4 -10.446246 -10.506591
    ## 59 TRINITY_DN311486_c4_g1  -7.339634 -10.033101
    ## 60 TRINITY_DN298167_c1_g1  -6.416088  -9.615547

``` r
range(nonzosteraSigMeasuresDownregulated$logP) #All downregulated genes should have negative p-values. These p-values are already negative, so no more modification needs to happen.
```

    ## [1] -1.222769e+01 -2.820227e-04

``` r
nonzosteraSigMeasuresDownregulated$signLogP <- nonzosteraSigMeasuresDownregulated$logP #Copy entries to a new column to be consistent with formatting for table with upregulated genes
```

#### Combine tables with corrected p-values

``` r
nonzosteraSigMeasuresCorrected <- rbind(nonzosteraSigMeasuresUpregulated, nonzosteraSigMeasuresDownregulated) #Combine tables by columns
nonzosteraSigMeasuresCorrected <- nonzosteraSigMeasuresCorrected[,-c(2:3)] #Remove unnecessary columns
head(nonzosteraSigMeasuresCorrected) #Confirm changes
```

    ##                     gene signLogP
    ## 1 TRINITY_DN296708_c0_g1 29.84743
    ## 2 TRINITY_DN312737_c1_g1 27.30081
    ## 3 TRINITY_DN294773_c0_g2 26.80754
    ## 4 TRINITY_DN293894_c0_g2 26.72528
    ## 5 TRINITY_DN317297_c2_g1 26.71153
    ## 6 TRINITY_DN233011_c0_g1 26.33312

``` r
write.csv(nonzosteraSigMeasuresCorrected, "2019-07-11-nonZostera-Table-of-Significance-Measures.csv", quote = FALSE, row.names = FALSE) #Save file
```

``` bash
head 2019-07-11-nonZostera-Table-of-Significance-Measures.csv #Confirm row names are not included in the file
```

    ## gene,signLogP
    ## TRINITY_DN296708_c0_g1,29.8474329572232
    ## TRINITY_DN312737_c1_g1,27.3008115384771
    ## TRINITY_DN294773_c0_g2,26.8075413257451
    ## TRINITY_DN293894_c0_g2,26.7252768716133
    ## TRINITY_DN317297_c2_g1,26.7115250404414
    ## TRINITY_DN233011_c0_g1,26.3331202997808
    ## TRINITY_DN295219_c1_g1,26.277331944757
    ## TRINITY_DN255273_c0_g1,26.2610866873347
    ## TRINITY_DN271642_c0_g1,26.1892909277872
