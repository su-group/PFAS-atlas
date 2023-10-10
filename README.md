
 # An Artificial Intelligence Platform for Automated PFAS Subgroup Classification: a Discovery Tool for PFAS Screening.
 * This repository has the data and code for our manuscript "An Artificial Intelligence Platform for Automated PFAS Subgroup Classification: a Discovery Tool for PFAS Screening". The manuscript is currently under peer review. PLEASE DO NOT SHARE THE CODE BEFORE OUR MANUSCRIPT IS PUBLISHED. Thank you!


 
# PFAS_Atlas

   - **Website link**
   - **How to cite**
   - **Installation**
   - **How to Use PFAS_Atlas**

**Website link**
     
    http://116.62.166.157/   or  www.pfas-atlas.net(which may not be accessible during the review process)

**How to cite**

    Please cite this paper as:Su, An & Cheng, Yingying & Zhang, Chengwei & Yang, Yun-Fang & She, Yuan-Bin & Rajan, Krishna. (2023). An Artificial Intelligence Platform for Automated PFAS Subgroup Classification: A Discovery Tool for PFAS Screening. 10.26434/chemrxiv-2023-4m96k. 



# If you want to run it locally, perform the following steps to deploy the environment.

**Installation**
    
    1.Do not worry about operating systems as PFAS_Atlas has been tested on Linux(Ubuntu LTS 20.04).
    2.It is essential to install Pycharm first if Pycharm is not installed on your computer.
      https://www.jetbrains.com/pycharm/download/#section=linux
      Click on the download page of the community edition, it is free and for pure Python development.
    3.Install Anaconda for Python 3.6 if your computer do not have Anaconda for Python 3.6.
      https://www.anaconda.com/products/individual
    4.Opening the linux Terminal,input:
      "conda create -n pfas python=3.6 -y"
      This will setup a conda virtual environment "pfas" that you will use to run PFAS_Atlas.
    5.When the pfas environment is properly setup, input:
      "conda activate pfas"
      Now you are in the pfas virtual environment.Let's install the basic machine learning and visualization libraries.
    6.To install rdkit, input:
      "conda install -c rdkit rdkit=2020.03.3 -y"
    7.To install tmap, input:
      “conda install -c tmap tmap -y"
    8.To install rxnfp, input:
      "pip install rxnfp"
    9.To install pytorch, input:
      "conda install pytorch==1.10.1 torchvision==0.11.2 torchaudio==0.10.1 cudatoolkit=11.3 -c pytorch -c conda-forge"
    10.To install mhfp, input:
      "pip install mhfp"
    11.Manually move "best_model" from the "outputs" directory to python's environment directory, such as ../anaconda3/envs/pfas/lib/python3.6/site-packages/rxnfp/models/transformers
    12.Start exploring PFAS Atlas! After entering PyCharm, click File→Open and select the folder to import the project 
      in the popup window.After the Python project is started, you need to configure the Python corresponding to the project
      to run properly.Follow these steps:
      "File → settings→ → Project → Python Interpreter → Add → conda environments → Existing environment → ok
    13.Start exploring PFAS_Atlas!If there is any problem in the installation, please contact ansu@zjut.edu.cn. 
     
    Note: This project can be used not only with Pycharm but also with another IDE or Command Line, but it is important to 
          switch to the Conda Environment first!

  
**How to Use PFAS_Atlas**
    
    1."step1_process.py": In order to pre-data processing.
    2."step2_get_mhfp.py": In order to compute molecular fingerprints.
    3."step3_classify.py": Classify OECD PFAS.
    4."step4_unsupervised_training.py": Unsupervised pre-training phase of bert model.
    5."step5_visualization_OECD_PFAS.py": Use tmap to reduce dimension and visualize OECD PFAS classification results.
    6. "step6_visualization_toxicity.py": Use tmap to explore hit ratio of PFAS toxicity.
       "step6_visualization_VDR.py": Use tmap to explore the binding affinity of PFAS to vitamin D receptors.
 
