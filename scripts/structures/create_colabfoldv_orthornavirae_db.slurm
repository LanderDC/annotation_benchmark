#!/bin/bash -l
#SBATCH --nodes="1"
#SBATCH -c "96"
#SBATCH --cluster=wice
#SBATCH --time="2-00:00:00"
#SBATCH -A lp_jm_virome_group
#SBATCH -p batch_sapphirerapids

cd $VSC_SCRATCH/colabfoldDBs/viral_db

conda activate colabfoldv_env

THREADS=96
mmseqs createdb colabfoldv_orthornavirae.faa  colabfoldv_orthornavirae_DB
mmseqs linclust colabfoldv_orthornavirae_DB  colabfoldv_orthornavirae_DB100 tmp --min-seq-id 1.0  -c 1.0 --threads $THREADS
mmseqs createclusearchdb  colabfoldv_orthornavirae_DB  colabfoldv_orthornavirae_DB100  colabfoldv_orthornavirae_DB100_search  --threads $THREADS
mmseqs linclust colabfoldv_orthornavirae_DB100_search  colabfoldv_orthornavirae_DB100_clu tmp --min-seq-id 0.3  -c 0.9 --threads $THREADS
mmseqs createclusearchdb  colabfoldv_orthornavirae_DB100_search  colabfoldv_orthornavirae_DB100_clu  colabfoldv_orthornavirae_DB100_clu_searchDB  --threads $THREADS
mmseqs align colabfoldv_orthornavirae_DB100_search  colabfoldv_orthornavirae_DB100_search colabfoldv_orthornavirae_DB100_clu  colabfoldv_orthornavirae_DB100_aln -e inf -a --threads $THREADS
