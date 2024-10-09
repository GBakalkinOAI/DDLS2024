# DDLS2024
Repository for the course https://ddls.aicell.io/course/ddls-2024

For assignment 1 I tried to find DDLS paper about scRNA-seq using chatGPT (chat link)[https://whatever] and Perplexity.

Hackaton: all files starting with CELLxGENE... are for the hackaton. 
I am not good at navigating "commits", so it is easier for me to write _v01, _v02, etc. in case I will need to borrow code from the earlier version to the latest one.

I am statistician / bioinformatician at Dept IGP of Uppsala University ( see https://www.uu.se/en/contact-and-organisation/staff?query=N12-1805 ).
My background is Appled Mathematics. 
Our group studies chromosomal alterations, like Loss of Y chromosome, which is detectable in every old male and affects his immune response.

My Final Project code will grow further after the deadline to reach these milestones:

- Learn to use harmonization and embeddings achieved by CELLxGENE Census scRNA-seq database. They have integrated much of published datasets anyway, including, I suspect, human scRNAseq datasets published by our group.
- Start by processing Monocyte-like cell types only to avoid crashing Colab out of RAM.
- Detect LOY cells by absense of Y chromosome genes expression. One has to threshold, because usually some chrX reads are mismapped to chrY.
- Plot scVI (Geneformer, etc.) embedding, color it with LOY detection, select LOY-enriched cell types.
- If time permits, resclaster LOY-enriched cell types to get finer subtypes.
- Detect LOY-related differentially expressed genes _withing each donor_ and theck how concordant they are across donors from different datasets. In https://pubmed.ncbi.nlm.nih.gov/33837451/ we called LOY DEGs as LOY-associated transcriptional effect (LATE) genes.
- Correct for LOY confounders, since LOY is strongly associated with age, and many diseases (cancers, AD, etc.), and environmental exposures (smoking, etc.)
- If time permits, check how these associations / confounding effects manifest themselves in LATE genes.
- Quick&dirty start in Google Colab.
- When the code stabilizes, create singularity container to share reproducible results with our collaborators and readers of the future papers.
