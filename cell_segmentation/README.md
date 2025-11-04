# SCRIPT FOR RUNNING BAYSOR 0.6.2 ON XENIUM ONBOARD ANALYZER (XOA) OUTPUT

## General overview:
1. Run cell-presegmentation on stained Xenium-slides. (optional)
   ==> creates binary masks for cell, that can be used by Baysor as input in addition to transcript coordinates

2. Run Baysor cell segmentation (with or without cell presegmentation masks) using XOA output (transcript coordinates)


## Steps:
1. Set up base directory, which should contain following folders:
   - xenium_data: directory containing RAW Xenium data and high-quality staining of Xenium slides in .TIFF format, unzipped
   - code_folder: root folder for code (name can be changed)
   - Set up the correct base_dir variable based on this in the script ```cs_run_cell_segm.sh``` and 

2. Preprocessing:
 - Run ```pp_rename_raw_folder_names.ipynb``` to rename raw folder names to reflect patient number and condition

3. Pull docker image containing Baysor 0.6.2 + pyometiff for reading ome.tiff files
   - Run ```docker pull gargerd/baysor:v0.6.2```
   - Extract image ID, and set it as the "dock_img_name" variable in start_baysor_docker_container() function, in the bash script ```bays_run_baysor_within_docker.sh```
   - NOTE: Baysor released a newer version (0.7.0) in September 2024, which is much faster and has  sligthly different hyperparameters. if you wish to use that version with this code, you have to build your own docker image. You can do it by changing the Baysor version installed in the Docker image. Just change the Baysor version in the ```dockerfile.txt```, and build your own Docker image.

4. Run presegmentation and create binary masks of cells in Xenium slide (optional), 3 methods:
   - Presegmentation with Cellpose (cyto or nucleus mode) ==> large GPU RAM requirements!
   - Presegmentation with default 10x XOA cell segmentation (no heavy computation involved)
   - No presegmentation 

   - Set up base_dir correctly + set up which methods should be used for presegmentation within ```cs_run_cell_segm.sh```
   - And then ```cs_run_cell_segm.sh```

5. Run Baysor using XOA input (+ presegmentation masks)
    - Within ```bays_run_baysor_within_docker.sh```:
       - Set up ```base_dir``` and ```dock_img_name``` variables 
       - Set up which presegmentation methods should be used as input
       - Set up hyperparameters of Baysor (util_run_baysor.py contains help with hyperparameters)
       - Currently loop is set up for multiple scale hyperparameter values ==> this can be also estimated by Baysor
    - Run ```bays_run_baysor_within_docker.sh```

   


