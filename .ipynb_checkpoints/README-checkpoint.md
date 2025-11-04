# scRNA-seq and Xenium spatial transcriptomic analysis of human carotid artery plaque samples

- List of python packages and Python version can be found in ```environment.yaml``` file
- At the beginning of each .Rmd file, the folder where R packages were installed is set. Change this path accordingly to your R-packge folder.
- Folder structure:
    - scRNA-seq alignment should be done on an HPC cluster. Change the the paths of ```raw_data_path``` and  ```output_path```  in ```run_cellranger_count.sh``` to your path if wish to run the alignment. 
    - Set up a ```base_dir```, which contains:
        - The raw Xenium data in a folder called ```xenium_data```
        - The root folder of the code (preferentially named ```src```)
        - The folder ```aligned_data```, which contains the output of the scRNA-seq alignment script, i.e. the aligned data

1. scRNA-seq
   
    1. Alignment
       - Code for the alignment is in the ```alignment_code``` folder.
       - First create a referemce genome with ```run_cellranger_mkref.sh```. If you wish, change the parameters of reference genome creation within the code.
       - Then run ```run_cellranger_count.sh```. If you wish, change the parameters of alignment within the code.
         
    2. Analysis of aligned data
        - Code of the scRNA-seq analysis is contained in ```scrna_analysis``` folder
        - Analysis of cell counts: notebooks starting with "ath" should be executed in the order of number indicated in the prefix of their names
        - The two Notebooks analysing the dataset of GSE159677 dataset have to be executed before starting  ath_8_1... .ipnyb Notebook
        - Information about the analysis steps can be found in the individual notebook

2. Xenium spatial transcriptomics
   
    1. Custom cell segmentation
        - Code is in the ```cell_segmentation``` folder with own Readme.me file
        - Detailed Readme file contains information about cell segmentation steps
         
    2. Processed spatial transcriptomic analysis
        - Code of spatiel transciptomic analysis can be found in  ```xenium_analysis``` folder
        - Notebooks starting with "xen" should be executed in the order of number indicated in the prefix of their names
        - Information about the analysis steps can be found in the individual notebook

