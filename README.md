# Fine-mapping pipeline

This repo can be used to perform fine-mapping with SuSiE and CARMA on a large, shared computing cluster.

### Important Folders:
	•	src: Contains scripts for running the pipeline.
	•	output: Contains the output of the finemapping process.
	•	data: Will contain summary statistics and lead SNPs.

### Installation
The pipeline requires the use of hail and Google Cloud SDK. We can install these in a conda environment by installing from the file hail_env.yml located in the base folder:

mamba env create -f hail_env.yml

This requires that you have conda and mamba installed (mamba will drastically improve installation times). 

Once the hail environment is installed you can activate it using 

conda activate hail_env

Once activated you will have to install CARMA directly in an R session:

R -e "remotes::install_github('ZikunY/CARMA')"

You will have to define project and region, and install gcloud beta in order to create a virtual machine on Google Cloud:

gcloud config set project $PROJECTID
gcloud config set compute/region us-central1
gcloud components install beta


### Running the pipeline

#### Activate conda environment

If you haven’t already, make sure to activate the conda environment for running the pipeline:

conda activate hail_env

Google Cloud Authentication

Before starting the pipeline you will have to login to gcloud:

gcloud auth login --no-launch-browser

You will receive a link to login to Google Cloud. Copy it into your browser and follow the instructions. Return the code after you have authenticated your login.

#### Pipeline

We should now be ready to run the actual pipeline. The pipeline consists of three scripts: 

1.0.0_prep_locus_ss.sh - Prepares loci and uploads them to cloud \
2.0.0_get_LDmat.py - Runs on the cloud to generate LD matrices \
3.1.0_run_FM_per_locus.sh - Performs finemapping with Susie and CARMA

##### 1.0.0_prep_locus_ss.sh

We can run the first script by navigating to the src folder and running:

sbatch src/1.0.0_prep_locus_ss.sh [/path/to/lead_SNP_file] [gcloud_sumstats_name] [/path/to/sumstats] [Ancestry] [Region size]

Arguments:
Path to lead snp file
Name of project to be created on google cloud
Path to summary statistics
Ancestry (we tend to use three-letter abbreviation like EUR for european)
Region size in Mbp (for 1.5Mbp region around each lead SNP you should put in 1.5, leading to a total window size of 3Mb)

Example:

sbatch src/1.0.0_prep_locus_ss.sh data/2016_27723757_VIT_EUR_leadSNPs.tsv 2016_27723757_VIT_EUR data/2016_27723757_VIT_EUR.assoc EUR 1.5

##### 2.0.0_get_LDmat.py

The second script needs to be run  on a virtual machine (VM) on Google Cloud. To initiate this we run:

CLOUD_PROJECT=nygc-comp-d-95c4
REGION=us-central1
NAME="$(whoami)-cluster" 
SERVICE_ACCOUNT=project-service-account@nygc-comp-d-95c4.iam.gserviceaccount.com

hailctl dataproc start $NAME \
  --project $CLOUD_PROJECT \
  --region $REGION \
  --no-address \
  --max-age 8h \
  --num-preemptible-workers 2 \
  --num-workers 2 \
  --max-idle 30m \
  --service-account $SERVICE_ACCOUNT \
  --subnet subnet-10-128

Once you have run this it will need a few minutes to initiate the VM. When this is done you can run the script to obtain the LD Matrix by running:

hailctl dataproc submit $NAME src/2.0.0_get_LDmat.py --ss_name [gcloud_sumstats_name] --window_mb [Region size]

The gcloud_sumstats_name must match the name you chose when running script 1.0.0_prep_locus_ss.sh


Example:
hailctl dataproc submit $NAME src/2.0.0_get_LDmat.py --ss_name 2016_27723757_VIT_EUR --window_mb 1.5


##### 3.0.0_run_FM.sh

Final step is to run the finemapping pipeline using both CARMA and SuSiE. The script can be initiated from the src directory using:

sbatch src/3.0.0_run_FM.sh [gcloud_sumstats_name] [Ancestry] [n samples] [n_cases] [Region size]

Example:
sbatch src/3.0.0_run_FM.sh 2016_27723757_VIT_EUR EUR 40258 2853 1.5

Here n samples refers to the total number of samples used in the GWAS study while n cases refers to only the number of cases. Both these numbers can be found in the GWAS catalog entry for the study of use.



## Cite

## Maintainer

Sophia Gunn @ sgunn@nygenome.org

## Acknowledgements

Thank you to Nathanael Andrews for assembling README!

## Release Notes
