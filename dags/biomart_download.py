from datetime import timedelta
from pathlib import Path

import subprocess
from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.dummy_operator import DummyOperator
from airflow.operators.python_operator import PythonOperator

from consts.biomart import BIOMART_BLAST_PATH, BIOMART_DATA_PATH, DATASET, FORCE_BIOMART_DOWNLOAD, MINIMAL_FILE_SIZE, \
    REGION_LIST, SERVER, \
    lncRNA_FILTER
from consts.dag import  default_args
from utils.logger import logger
from pandas import DataFrame

from utils.utils import DirectorySpecificBashOperator, fasta_to_dataframe, filter_Sequenceunavailable_from_fasta


def biomart_query(server: str, dataset: str, region: str, fout: str) -> str:
    def get_filter(attribute: str) -> str:
        if attribute == "lncRNA":
            return lncRNA_FILTER
        return ""
    def get_region(attribute: str) -> str:
        if attribute == "lncRNA":
            return "transcript_exon_intron"
        return attribute

    return f"wget -O {fout} \'{server}/biomart/martservice?query=" \
           f"<?xml version=\"1.0\" " \
           f"encoding=\"UTF-8\"?> " \
           f"<!DOCTYPE Query> <Query virtualSchemaName = \"default\" " \
           f"formatter = \"FASTA\" header = \"0\" uniqueRows = \"0\" " \
           f"count = \"\" datasetConfigVersion = \"0.6\" >  " \
           f"<Dataset name = \"{dataset}\" " \
           f"interface = \"default\" > " \
           f"{get_filter(region)} " \
           f"<Attribute name = \"ensembl_gene_id\" /> " \
           f"<Attribute name = \"ensembl_transcript_id\" />  " \
           f"<Attribute name = \"{get_region(region)}\" /> " \
           f"</Dataset> </Query>\'"


biomart_download = DAG(
    "biomart_download",
    default_args=default_args,
    description='download sequences from biomart',
    # schedule_interval=timedelta(days=30),
)

def run_biomart_query(organism: str, region: str, force: bool=True):
    output_file: Path = BIOMART_DATA_PATH / f"{organism}_{region}.fa"
    if not force:
        if output_file.exists():
            if output_file.stat().st_size > MINIMAL_FILE_SIZE:
                return
    cmd = biomart_query(SERVER, DATASET[organism], region, str(output_file))
    print(biomart_query(SERVER, DATASET[organism], region, BIOMART_DATA_PATH / f"{organism}_{region}.fa"))
    return subprocess.call(cmd, shell=True)


def create_biomart_python_operator(organism: str, region: str) -> PythonOperator:
    return PythonOperator(
        task_id=f"{organism}_{region}",
        python_callable=run_biomart_query,
        op_kwargs={'organism' : organism,
                   'region': region,
                   'force': FORCE_BIOMART_DOWNLOAD},
        dag=biomart_download)


def create_filter_operator(organism: str, region: str) -> PythonOperator:
    return PythonOperator(
        task_id=f"Fasta_Filter_{organism}_{region}",
        python_callable=filter_Sequenceunavailable_from_fasta,
        op_kwargs={'fasta_file': BIOMART_DATA_PATH / f"{organism}_{region}.fa"},
        dag=biomart_download
    )


def create_blast_db(organism: str, region: str) -> PythonOperator:
    fasta_in = BIOMART_DATA_PATH / f"{organism}_{region}.fa"
    db_title = f"{organism}_{region}"
    return DirectorySpecificBashOperator(
        task_id=f"blast_db_{organism}_{region}",
        cmd=f"makeblastdb -in {fasta_in} -parse_seqids -dbtype nucl -out {db_title}",
        dag=biomart_download,
        cwd=BIOMART_BLAST_PATH)
   

def create_organism_dag(organism: str):
    tasks = [create_biomart_python_operator(organism, region) >>
             create_filter_operator(organism, region) >>
             create_blast_db(organism, region)
             for region in REGION_LIST]

    done =  DummyOperator(
        task_id=f"{organism}_done",
        dag=biomart_download,
        )
    return tasks, done


all_done = DummyOperator(
        task_id="all_done",
        dag=biomart_download,
        )


human_tasks, human_done = create_organism_dag("human")
human_tasks >> human_done

mouse_tasks, mouse_done = create_organism_dag("mouse")
mouse_tasks >> mouse_done


celegans_tasks, celegans_done = create_organism_dag("celegans")
celegans_tasks >> celegans_done


cattle_tasks, cattle_done = create_organism_dag("cattle")
cattle_tasks >> cattle_done


def csvConverter(dir: Path):
    for f in dir.glob("*.fa"):
        logger.info(f"convert fasta to csv: {f}")
        df: DataFrame = fasta_to_dataframe(f)
        try:
            df["sequence length"] = df["sequence"].apply(lambda x: len(x))
        except KeyError:
            #empty file
            pass
        outfile = f.parent / (f.stem + ".csv")
        df.to_csv(outfile)

#
# def fastaCleaner(dir: Path):
#     for f in dir.glob("*.fa"):
#         logger.info(f"filter fasta file: {f}")
#         filter_Sequenceunavailable_from_fasta(f)
#
#
# fastaFilter = PythonOperator(
#     task_id='Fasta_Filter',
#     python_callable=fastaCleaner,
#     op_kwargs={'dir': BIOMART_DATA_PATH},
#     dag=biomart_download,
# )


convert_to_csv = PythonOperator(
    task_id='convert_to_csv',
    python_callable=csvConverter,
    op_kwargs={'dir': BIOMART_DATA_PATH},
    dag=biomart_download,
)


#
# [human_done, mouse_done, celegans_done, cattle_done] >> fastaFilter >> convert_to_csv >> all_done


[human_done, mouse_done, celegans_done, cattle_done] >> convert_to_csv >> all_done

