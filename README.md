# sMRD Code

This code is used for analysis and visualization of data from the sMRD manuscript by Ransohoff et al from the laboratory of Ash Alizadeh at Stanofrd University. 

# Installation

Install [anaconda or miniconda.](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)

Then, create and activate a conda environment:
```
conda env create -f environment.yml
conda activate smrd
```

Installation should complete within a couple of minutes on a standard desktop computer.

# Usage

Synthetic example data is provided under:
```
data/sample_info/ex_sample_info.csv -- Sample metadata
data/features/ex_cnv.csv -- CNV calls
data/features/ex_snv.csv -- SNV calls
data/features/ex_enrichment.csv -- Raw Variant Calls
```

To generate intermediate data and figures, run the following script:
```
./run_analysis.sh
```

To perform variant enrichment analysis, run the R Markdown notebook under `analysis/enrichment_analysis.Rmd`.

To reproduce the manuscript figures or to run a custom analysis, replace the provided synthetic data with a new dataset.

# Output

Tabular output, phylogenetic trees, and statistics are written to the `data/output` directory. Plots are written to the `data/plots` directory. Runtime for example data should not exceed a couple of minutes on a standard desktop computer.

# License Terms of Use
The Board of Trustees of the Leland Stanford Junior University (“Stanford”) provides sMRD features and services (“Service”) free of charge for non-commercial use only. By using the Service, you agree to be bound by the terms of this Agreement. Please read it carefully.
1. You agree to use the Service solely for research purposes. You agree not to use the Service for diagnosis, treatment, cure, prevention, or mitigation of disease or other conditions in man or other animals.
2. You acknowledge that the Service and any information obtained therefrom is not intended to substitute for care by a licensed healthcare professional.
3. You represent and warrant that you are not an employee of a for-profit entity, and that your use of the Service is not directed towards commercial advantage or for-profit activities. Commercial entities wishing to use this Service should contact Stanford University’s Office of Technology Licensing.
4. THE SERVICE IS OFFERED “AS IS”, AND, TO THE EXTENT PERMITTED BY LAW, STANFORD MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED. STANFORD SHALL NOT BE LIABLE FOR ANY CLAIMS OR DAMAGES WITH RESPECT TO ANY LOSS OR OTHER CLAIM BY YOU OR ANY THIRD PARTY ON ACCOUNT OF, OR ARISING FROM THE USE OF THE SERVICE. YOU HEREBY AGREE TO DEFEND AND INDEMNIFY STANFORD, ITS TRUSTEES, EMPLOYEES, OFFICERS, STUDENTS, AGENTS, FACULTY, REPRESENTATIVES, AND VOLUNTEERS (“STANFORD INDEMNITEES”) FROM ANY LOSS OR CLAIM ASSERTED AGAINST STANFORD INDEMNITEES ARISING FROM YOUR USE OF THE SERVICE.
5. You agree not to use (i) Stanford’s name or other trademarks, or (ii) the name of any Stanford faculty member, employee, student or volunteer without the prior written consent of Stanford. Permission may be withheld at Stanford’s sole discretion. This prohibition includes, but is not limited to, use in press releases, advertising, marketing materials, other promotional materials, presentations, case studies, reports, websites, application or software interfaces, and other electronic media.
6. All rights not expressly granted to you in this Agreement are reserved and retained by Stanford or its licensors or content providers.
7. You agree that this Agreement and any dispute arising under it is governed by the laws of the State of California, United States of America, applicable to agreements negotiated, executed, and performed within California.
8. Subject to your compliance with the terms and conditions set forth in this Agreement, Stanford grants you a royalty-free, non-exclusive, non-sublicensable, and non-transferable license to access and make use of the Service.
