Data and supplementary tables for Audrey Ward's 2025 dissertation "Developing the wine yeast _Lachancea thermotolerans_ as a model for evolutionary genomics" at the University of Georgia under the direction of Douda Bensasson.

**Citation:** Ward, A.K. (2025). Developing the wine yeast _Lachancea thermotolerans_ as a model for evolutionary genomics (Publication No. 31846584)[Doctoral dissertation, University of Georgia]. ProQuest Dissertations and Theses Global.

## Supplementary data tables (within WardDissertation25_supplementaryTables.xlsx):

### Chapter 2 - Geography, domestication, and population genomics of a wine yeast, _Lachancea thermotolerans_

Table 2.S1. Overview of _Lachancea thermotolerans_ and _Lachancea quebecensis_ strains used in the study, including compiled metadata from publicly available and new whole-genome data. Below is a description of each column:
- Strain Name: Strain name, as defined in the original supplementary material for each cited study.
- Secondary Strain Name: Strain name with alternate capitalization, if any, that appeared during analyses.
- Abbreviated Strain Name: Strain name without dot, dash, underscore, or space characters.
- Species: Identified species (_L. thermotolerans_ or _L. quebecensis_)
- Genome Data: NCBI Bioproject identifier.
- Geographic Origin: Location where each strain was isolated, as defined in original supplementary material for each cited study.
- Region: Continent of origin, based on geographic location information.
- Subregion: Assigned area of origin based on geographic location, region, and proximity to other isolate sites represented.
- Ecological Origin: General ecology of isolate origin, as defined in original supplementary material.
- Mean Genome Coverage: Calculated mean genome coverage across the whole genome after read mapping using average read depth from consensus call files.
- Published Clade Assignment: Clade and citation for previously defined lineages.
- This Study Clade Assignment: Lineage as defined within this study.
- Study: Isolate and genome source citations for each strain, if different. 

Table 2.S2. Library preparation methods and sequencing platforms used during this study.
Below is a description of each column:
- Strain Name: Strain name, as defined in the original supplementary material for each cited study.
- Secondary Strain Name: Strain name with alternate capitalization, if any, that appeared during analyses.
- Abbreviated Strain Name: Strain name without dot, dash, underscore, or space characters.
- Species: Identified species (_L. thermotolerans_ or _L. quebecensis_)
- Sequencing Order: Sequencing center and associated project number.
- Library Prep: Preparation kit used for sequencing libraries.
- Sequencing Platform: Brand and sequencing platform used for each sequence.
- Runs by Read Length: Number of sequencing runs and size (base pairs) of each read.
- Paired End: Whether sequences were paired or not (yes/no).

Table 2.S3. Strains used in population structure analyses. Below is a description of each column:
- Strain Name: Strain name, as defined in the original supplementary material for each cited study.
- Secondary Strain Name: Strain name with alternate capitalization, if any, that appeared during analyses.
- Abbreviated Strain Name: Strain name without dot, dash, underscore, or space characters.
- Geographic Origin: Location where each strain was isolated, as defined in original supplementary material for each cited study.
- Region: Continent of origin, based on geographic location information.
- Subregion: Assigned area of origin based on geographic location, region, and proximity to other isolate sites represented.
- Ecological Origin: General ecology of isolate origin, as defined in original supplementary material. 
- Study: Isolate and genome source citations for each strain, if different. 
- Published Clade Assignment: Clade and citation for previously defined lineages.
- Admixed: A strain is defined as being admixed if percent ancestry from a single population is < 90% from ADMIXTURE analysis when K = 15.

Table 2.S4. f4 test results for significant migration events. Table shows f4 test ((A,B);(C,D)) results for significant (P < 0.01) edges across TreeMix graphs for all values of M (M0-M10). Below is a description of each column:
- Migration Event: Direction of gene flow (SourcePopulation_TargetPopulation).
- Pop A: Population in the A position in the f4 test described above.
- Pop B: Population in the B position in the f4 test described above.
- Pop C: Population in the C position in the f4 test described above.
- Pop D: Population in the D position in the f4 test described above.
- f4 Stat: TreeMix calculated f4 value.
- f4 SE: Standard error for calculated f4 value.
- Z Score: Calculated Z Score for each f4 value.
- Significant Z Score: Z Score was defined as being significant if its value was ≤ -3.00 or ≥ 3.00.

Table 2.S5. Genes used in copy-number variation analyses. List of gene names and positions within the _L. thermotolerans_ reference genome with names of corresponding _Saccharomyces cerevisiae_ homologs. Below is a description of each column:
- Gene Name: _Lachancea thermotolerans_ systematic name for each gene of interest.
- Sc Standard Name: _Saccharomyces cerevisiae_ standard name for each gene of interest.
- Sc Systematic Name: _Saccharomyces cerevisiae_ systematic name for each gene of interest.
- Reference: Citation originally mentioning gene as showing copy number variation within _L. thermotolerans_.
- Chr Letter: Identifying letter (A-H) for each _L. thermotolerans_ chromosome.
- Chr Number: Corresponding number for each _L. thermotolerans_ chromosome; needed as input for analysis.
- Position Start (bp): Coordinates for gene starting position in chromosome.
- Position End (bp): Coordinates for gene ending position in chromosome.

Table 2.S6. Copy number of genes of interest per strain. Copy number of 16 _Saccharomyces cerevisiae_ homologs in _Lachancea thermotolerans_ strains and description of strain clade and ecological origin. Below is a description of each column:
- Strain: Strain name, as defined in the original supplementary material for each cited study and used here.
- Clade: Lineage as defined within this study.
- Ecological Category: From original ecological information (Table 2.S1), categorizing origin as being from animal, cactus, crop, fermentation, flower, food, fruit, grape or wine associated, or trees.
- Domestic: Defined based on proximity of ecological category to human activity.
- Admixture: A strain is defined as being admixed if percent ancestry from a single population is < 90% from ADMIXTURE analysis when K = 15.
- ATO3 - YKL107W: Number of copy numbers identified for each gene using our methods.

 Table 2.S7. Marginal likelihood summary of ancestral character estimation model selection. Summary of marginal likelihood (ml) results from the equal rates (ERM) and free rates (freeK) models of evolution for ancestral character estimation in RevBayes calculated via stepping-stone (ss) and path sampling (ps) methods. Includes calculated Bayes Factors (BF). Below is a description of each column:
- Gene: Homologous gene of interest.
- ERM Convergence: Equal Rates Model calculated in RevBayes converged.
- Free K Convergence: Free Rates model calculate in RevBayes converged.
- ERM (ml, ss): Marginal likelihood results for stepping-stone analysis using the ERM.
- ERM (ml, ps): Marginal likelihood results for path sampling analysis using the ERM.
- Free K (ml, ss): Marginal likelihood results for stepping-stone analysis using freeK.
- Free K (ml, ps): Marginal likelihood results for path sampling analysis using freeK.
- BF (ml, ss): Calculated Bayes Factor between the stepping-stone analyses for marginal likelihoods of both models. 
- BF (ml, ps): Calculated Bayes Factor between the path sampling analyses for marginal likelihoods of both models.
- Model Pref BF: Preferred model based on calculated Bayes Factors.
- Model Pref Rates Graph: Preferred model based on comparisons of modeled rates of evolution (Figure 2.S13).

Table 2.S8. Strains used in lineage divergence or nucleotide diversity estimates. Summary of strains used to calculated nucleotide diversity estimates in MEGA11 for _L. thermotolerans_, _S. cerevisiae_, and _S. paradoxus_. Below is a description of each column:
- Strain: Strain name, as defined in the original supplementary material for each cited study.
- Species: Identified species.
- Clade: Defined lineage for each strain.
- Clade Source: Citation for strain clade assignments.
- Genome Source: Citation for strain sequence origin.

### Chapter 3 - Growth rate differences at high temperatures are driven by genetic background in wild _Lachancea thermotolerans_

Table 3.S1. _Lachancea thermotolerans_ strains used in this study. Metadata was compiled from publicly available data. Includes geographic and ecological origin, as well as genetic background (clade, admixture status, divergence event) and the phenotyping plate for each strain. Below is a description of each column:
- Strain Name: Strain name, as defined in the original supplementary material for each cited study.
- Geographic Origin: Location where each strain was isolated, see ‘Geographic Origin’ in Table 2.S1.
- Continent: Continent of isolate origin, see ‘Region’ in Table 2.S1.
- Subregion: Assigned area of origin; see Table 2.S1.
- Ecological Origin: General ecology of isolate origin, as defined in original supplementary material.
- Clade: Genetic lineage as defined in Chapter 2.
- Admixed: A strain is defined as being admixed if percent ancestry from a single population is < 90% from ADMIXTURE analysis when K = 15 (Chapter 2).
- Divergence Event: Categorization of tree branching; see Chapter 2.
- Phenotype Plate: 96-well plate strain in which strain was included for high-throughput thermal growth assay.
- Study Source: Isolate source citations for each strain.

Table 3.S2. Calculated growth rate data from all replicates. Calculated growth rate data from the linear regression analysis of growth curves across all replicates and temperatures. Below is a description of each column:
- Strain: Name of strain being experimented upon.
- Clade: Genetic lineage as defined in Chapter 2.
- Clade Simplified: Name of genetic lineage after statistical simplification.
- Continent: Continent of isolate origin.
- Subregion: Assigned area of origin.
- Tmax: Maximum temperature (degrees Celsius) during the hottest month for each isolate location according to WorldClim2 dataset.
- Precip Warmest Quarter: Precipitation (mm) during the warmest quarter of the year for each isolate location according to WorldClim2 dataset.
- Plate: Phenotyping plate strain included within for specific replicate.
- Temp: Experimental temperature condition.
- Rep: Replicate number.
- mslope: Slope at exponential growth.
- lag: Lag to growth.
- max: Maximum growth value (optical density or OD).
- rsq: R-squared as a measure of model fit.
- Any Growth: Binary character assigned according to whether there was growth based on maximum growth value (‘max’) of OD > 0.36.

Table 3.S3. Mean growth information. Mean value of maximum growth (‘meanMax’) and growth rate (‘meanMSlope’) across all replicates of temperatures 30C-39C for each strain, excluding negatives. Below is a description of each column:
- Strain: Name of strain being experimented upon.
- Clade: Genetic lineage as defined in Chapter 2.
- Clade Simplified: Name of genetic lineage after statistical simplification.
- Temp: Experimental temperature condition.
- Mean Max: Mean value of maximum growth across all replicates for temperatures 30-39C.
- Mean M Slope: Mean value of growth rate across all replicates for temperatures 30-39C.
- Any Growth: Binary character assigned according to whether mean maximum growth value OD > 0.36.

Table 3.S4:  Climatic variables do not significantly impact growth rate in _Lachancea thermotolerans_. P-values for Spearman correlation tests for relationships between growth rate and two climatic variables: max summer temperature (Tmax) and precipitation during the warmest quarter (precipitation) when grouped by genetic lineage.

### Chapter 4 - Low level contamination confounds population genomic analysis

Table 4.S1. Differences among studies in rates of intra-species _Saccharomyces cerevisiae_ contamination. Data from Peña et al (2025) after excluding genomes with low read depth and studies with fewer than 10 high depth genomes. Rates appear different (Fisher's exact test, P = 2 × 10-6) even after excluding PRJEB11698 (Fisher's exact test, P = 0.002).

Table 4.S2. Intra-species contamination does not lower the quality of base calls.

Table 4.S3. Summary of the strains used for _in silico_ contamination and population genomic analyses.
