#ensure you are in your home directory - this script assumes starting in /home/user
# Submit with: sbatch charr_cutadapt.sh

## Load python and launch virtual environment
module load python
source ~/python_env/bin/activate
cd /home/tkess/scratch/CharrGDL_WG

### 1. TRIM EVERYTHING
ls HI.*R1.fastq.gz | \
  sed 's/R1.fastq.gz//'| \
  parallel -j 20 'cutadapt \
    -u 15 \
    --minimum-length 40 \
    -q 10 \
    -A CTGTCTCTTATACACA \
    -a CTGTCTCTTATACACA \
    -G GATGTGTATAAGAGACAG  \
    -g GATGTGTATAAGAGACAG \
    -o {}R1_trim.fastq.gz \
    -p {}R2_trim.fastq.gz \
    {}R1.fastq.gz \
    {}R2.fastq.gz '


###2.  THEN STOP
deactivate



