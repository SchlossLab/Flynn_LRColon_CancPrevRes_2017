# Spatial variation of the native colon microbiota in healthy adults

The microbiome has been implicated in the development of colorectal cancer (CRC) and inflammatory bowel diseases (IBD). The specific traits of these diseases vary along the axis of the digestive tract. Further, variation in the structure of the gut microbiota has been associated with both diseases. Here we profiled the microbiota of the healthy proximal and distal mucosa and lumen to better understand how bacterial populations vary along the colon. We used a two-colonoscope approach to sample proximal and distal mucosal and luminal contents from the colons of 20 healthy subjects that had not undergone any bowel preparation procedure. The biopsies and home-collected stool were subjected to 16S rRNA gene sequencing and Random Forest classification models were built using taxa abundance and location to identify microbiota specific to each site. The right mucosa and lumen had the most similar community structures of the five sites we considered from each subject. The distal mucosa had higher relative abundance of _Finegoldia, Murdochiella, Peptoniphilus, Porphyromonas_ and _Anaerococcus_. The proximal mucosa had more of the genera _Enterobacteriaceae, Bacteroides_ and _Pseudomonas_. The classification model performed well when classifying mucosal samples into proximal or distal sides (AUC = 0.808). Separating proximal and distal luminal samples proved more challenging (AUC = 0.599) and specific microbiota that differentiated the two were hard to identify. By sampling the unprepped colon, we identified distinct bacterial populations native to the proximal and distal sides. Further investigation of these bacteria may elucidate if and how these groups contribute to different disease processes on their respective sides of the colon.

Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            # raw and primary data, are not changed once created
    |  |- references/  # reference files to be used in analysis
    |  |- raw/         # raw data, will not be altered
    |  |- mothur/      # mothur processed data
    |  +- process/     # cleaned data, will not be altered once created;
    |                  # will be committed to repo
    |
    |- code/           # any programmatic code
    |  |- old_code     # old or prelimary branches of analysis
    |  |- PBS_scripts  # scripts for running on the HPCC
    |
    |- doc/            # protocol, IRB and other documents
    |
    |- submission/     # submission documents, including figures and responses to reviewers
    |  |- manuscript.* # executable Rmarkdown and pdf for this study
    |
