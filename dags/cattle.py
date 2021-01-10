from datetime import timedelta
# The DAG object; we'll need this to instantiate a DAG
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from consts.dag import CATTLE, default_args


#instantiates a directed acyclic graph
from consts.pipeline_steps import CONCAT_BLAST, MIRNA_SEQ_PATH, NORMALIZATION_PATH, REGION_PATH
from utils.blast import get_blast_graph

cattle_dag = DAG(
    CATTLE,
    default_args=default_args,
    description=f"Parse the paper: {CATTLE}",
    # schedule_interval=timedelta(days=30),
)


file_name = "cattle.csv"
path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"

paper_read = BashOperator(
    task_id='paper_read',
    bash_command=f"python {path}papers/{CATTLE}/reader.py {file_name}",
    dag=cattle_dag,
)


mirna_seq_insertion = BashOperator(
    task_id='mirna_seq_insertion',
    bash_command="python " + step_path + "mirna_seq_insertion.py " + file_name,
    dag=cattle_dag,
)

# rna_site_insertion = BashOperator(
#     task_id='rna_site_insertion',
#     bash_command=f"python {step_path}rna_site_insertion.py insert-site-from-chromosome "
#                  f"{MIRNA_SEQ_PATH}/{file_name} "
#                  f"{data_step_path}site/{file_name} "
#                  f"{path}data/genome/BosTau7 ",
#     dag=cattle_dag,
# )

blast_start, blast_tasks, blast_end = get_blast_graph(cattle_dag,
                                                      (MIRNA_SEQ_PATH / file_name),
                                                      "cattle")

concat_blast = BashOperator(
    task_id='concat_blast',
    bash_command=f"python {step_path}concat_blast_result.py concat-blast-result "
                 f"{REGION_PATH}/ "
                 f"cattle "
                 f"{MIRNA_SEQ_PATH / file_name} "
                 f"{CONCAT_BLAST / file_name} ",
    dag=cattle_dag,
)

normalization = BashOperator(
    task_id='normalization',
    bash_command=f"python {step_path}normalization_final_step.py finalize "
                 f"{CONCAT_BLAST / file_name} "
                 f"{NORMALIZATION_PATH / file_name} ",
    dag=cattle_dag,
)




#
# rna_region_insertion_using_join = BashOperator(
#     task_id='rna_region_insertion_using_join',
#     bash_command=f"python {step_path}rna_region_insertion.py human-mapping-run "
#                  f"{data_step_path}read/{file_name} "
#                  f"{data_step_path}region/{file_name} ",
#     dag=cattle_dag,
# )
#
# rna_site_insertion = BashOperator(
#     task_id='rna_site_insertion',
#     bash_command=f"python {step_path}rna_site_insertion.py get-site-from-extended-site "
#                  f"{data_step_path}region/{file_name} "
#                  f"{data_step_path}site/{file_name} ",
#     dag=cattle_dag,
# )
#
# normalization = BashOperator(
#     task_id='normalization',
#     bash_command=f"python {step_path}normalization_final_step.py finalize "
#                  f"{data_step_path}site/{file_name} "
#                  f"{data_step_path}normalization_final/{file_name} ",
#     dag=cattle_dag,
# )

paper_read >> mirna_seq_insertion >> blast_start
blast_end >> concat_blast >> normalization
# paper_read >> mirna_seq_insertion >> rna_site_insertion >> blast_start

# rna_region_insertion_using_join >> rna_site_insertion >> normalization

