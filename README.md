# DDLS2024
Final Project (by Daniil Sarkisyan) for the course https://ddls.aicell.io/course/ddls-2024

__Working project's part:__ Run these self-contained .ipynb in Google Colab in order listed

- `CellxGene_Census_scVI_Monocytes_read_08.ipynb`  - read metadata, filter cells, read scRNAseq and pre-trained embeddings

- `CellxGene_Census_scVI_Monocytes_process_02.ipynb` - process ingested metadata

__Failed to do:__ Teach chatGPT about breaking API changes introduced since GPT-4o (o1-PREVIEW, o1-mini) knowledge cut-off

- Latest release (v.1.16.2) was downloaded from https://github.com/chanzuckerberg/cellxgene-census

- `get_ipynb_py_md.sh` - converted .ipynb into .py, copied files into flat folder (see `cellxgene-census-1.16.2-inject-tutorials-as-py`), dropping identical files caused by links

- `combine_py_md.sh` - combine .md and .py tutorials to obey 10 files limit of chatGPT (see `cellxgene-census-1.16.2-attach-tutorials-to-chat`)

- Attaching `CELLxGENE_Census_combined.md` and `CELLxGENE_Census_combined.py` to chat did not stop chatGPT from using obsolete API syntax and functions

- `00_prompt_v01.txt` I asked o1-PREVIEW to syntesize plain-text API tutorial, to be injected into future prompts, by ingesting separate CELLxGENE documentation.
  Will probably work in the future, with bigger context windows and higher usage limits available.
  For now I would either run out of usage allowance, or GPT-o1-PREVIEW will take longer and longer to answer until it will fail

- `00_prompt_v02.txt` I tried both o1-PREVIEW and o1-mini (supposedly better at coding) to synthesize tutorial "inside", responding with just a summary of changes they had to make.
  They both did not complain until I will ask to print the full tutorial - then they will both refuse the request, citing the response size limit set by OpenAI.
  Asking them to break response into parts and output each part failed (only the first part was OK, but they did not continue).

- `00_prompt_v03.txt` Investigation of OpenAI developer forum found a lot of developer struggling (and failing) with adding "custom knowledge" to chatGPT...

- https://github.com/openai/openai-cookbook/blob/main/examples/Question_answering_using_embeddings.ipynb from OpenAI wise advice: _Although fine-tuning can feel like the more natural option — training on data is how GPT learned all of its other knowledge, after all — we generally do not recommend it as a way to teach the model knowledge. Fine-tuning is better suited to teaching specialized tasks or styles, and is less reliable for factual recall._

__Summarizes several years of my own ML failures all too well__

- `bitter_lesson.pdf` - I wish I would learn it sooner


I am not good at navigating commits / switching branches / temporarily returning to old stages, so it is easier for me to write _v01, _v02, etc. in case I will need to borrow code from the earlier version to the latest one.

I am statistician / bioinformatician at Dept IGP of Uppsala University ( see https://www.uu.se/en/contact-and-organisation/staff?query=N12-1805 ).
My background is Appled Mathematics. 
Our group studies chromosomal alterations, like Loss of Y chromosome, which is detectable in every old male and affects his immune response.

Final Project code will grow, extending beyond Project Proposal after the 2024-10-27 deadline to reach these milestones:

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
