from datetime import timedelta
# The DAG object; we'll need this to instantiate a DAG
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.python_operator import PythonOperator
from pandas import DataFrame

from consts.dag import DYNAMIC_FLY_GAMBIAE, default_args
import numpy as np

#instantiates a directed acyclic graph
from utils.utils import read_csv, to_csv

dynamic_fly_dag = DAG(
    DYNAMIC_FLY_GAMBIAE,
    default_args=default_args,
    description='Parse the paper: Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae',
    # schedule_interval=timedelta(days=30),
)


file_name = "dynamic_fly.csv"
path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"

paper_read = BashOperator(
    task_id='paper_read',
    bash_command="python " + path + "papers/Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae/reader.py " + file_name,
    dag=dynamic_fly_dag,
)


mirna_seq_insertion = BashOperator(
    task_id='mirna_seq_insertion',
    bash_command="python " + step_path + "mirna_seq_insertion.py mirna-seq-insertion " + file_name,
    dag=dynamic_fly_dag,
)



mRNA_insertion = BashOperator(
    task_id='mRNA_insertion',
    bash_command=f"python {step_path}mRNA_seq_insertion.py insert-from-fasta "
                 f"{path}data/gambiae/gambiae_utr.fa "
                 f"{data_step_path}mirna_sequence/{file_name} "
                 f"{data_step_path}mrna_sequences/{file_name} ",
    dag=dynamic_fly_dag,
)

#
# rna_site_insertion = BashOperator(
#     task_id='rna_site_insertion',
#     bash_command=f"python {step_path}rna_site_insertion.py insert-site-by-coordinates "
#                  f"{data_step_path}mrna_sequences/{file_name} "
#                  f"{data_step_path}site/{file_name} ",
#     dag=dynamic_fly_dag,
# )

rna_region_insertion = BashOperator(
    task_id='rna_region_insertion',
    bash_command=f"python {step_path}rna_region_insertion.py gambiae-run "
                 f"{data_step_path}mrna_sequences/{file_name} "
                 f"{data_step_path}region/{file_name} ",
    dag=dynamic_fly_dag,
)

def df_col_rename(fname: Path):
    df: DataFrame = read_csv(fname)
    df.rename(columns={'region_sequence' : 'sequence',
               # 'chimera_start' : 'start',
               # 'chimera_end' : 'end'
                       },
              inplace=True)
    print(df.columns)
    df.insert(0, 'region count', np.nan)
    df.insert(0, 'identity', np.nan)
    df.insert(0, 'coverage', np.nan)
    to_csv(df, fname)


df_rename = PythonOperator(
    task_id=f"df_col_rename",
    python_callable=df_col_rename,
    op_kwargs={'fname': Path(f"{data_step_path}region/{file_name}")},
    dag=dynamic_fly_dag)

normalization = BashOperator(
    task_id='normalization',
    bash_command=f"python {step_path}normalization_final_step.py finalize "
                 f"{data_step_path}region/{file_name} "
                 f"{data_step_path}normalization_final/{file_name} ",
    dag=dynamic_fly_dag,
)



#sets the ordering of the DAG. The >> directs the 2nd task to run after the 1st task. This means that
#download images runs first, then train, then serve.
paper_read >> mirna_seq_insertion >> mRNA_insertion >> rna_region_insertion >> df_rename >> normalization



# paper_read >> mirna_seq_insertion >> mRNA_insertion >> rna_site_insertion >> rna_region_insertion >> normalization

