from datetime import timedelta
# The DAG object; we'll need this to instantiate a DAG
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from consts.dag import QCLASH_MELANOMA, default_args


#instantiates a directed acyclic graph
from consts.pipeline_steps import CONCAT_BLAST, MIRNA_SEQ_PATH, NORMALIZATION_PATH, REGION_PATH, SITE_PATH
from utils.blast import get_blast_graph

qclash_melanoma_dag = DAG(
    QCLASH_MELANOMA,
    default_args=default_args,
    description=f"Parse the paper: {QCLASH_MELANOMA}",
    # schedule_interval=timedelta(days=30),
)



path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"

file_name = f"qclash_melanoma_human.csv"


paper_read = BashOperator(
    task_id=f"paper_read",
    bash_command=f" python /home/local/BGU-USERS/benorgi/pipeline/papers/Cross-Linking_Ligation_and_Sequencing_of_Hybrids_\(qCLASH\)_Reveals_an_Unpredicted_miRNA_Targetome_in_Melanoma_Cells/reader.py "
                 f"{file_name}",
    dag=qclash_melanoma_dag,
)

mirna_seq_insertion = BashOperator(
    task_id=f"mirna_seq_insertion",
    bash_command="python " + step_path + "mirna_seq_insertion.py qclash-melanoma-mirna-seq-insertion " + file_name,
    dag=qclash_melanoma_dag,
)

blast_start, blast_tasks, blast_end = get_blast_graph(qclash_melanoma_dag,
                                                      (MIRNA_SEQ_PATH / file_name),
                                                      organism="human")

concat_blast = BashOperator(
    task_id=f"concat_blast",
    bash_command=f"python {step_path}concat_blast_result.py concat-blast-result "
                 f"{REGION_PATH}/ "
                 f"{Path(file_name).stem} "
                 f"{MIRNA_SEQ_PATH / file_name} "
                 f"{CONCAT_BLAST / file_name} ",
    dag=qclash_melanoma_dag,
)

normalization = BashOperator(
    task_id=f"normalization",
    bash_command=f"python {step_path}normalization_final_step.py finalize "
                 f"{CONCAT_BLAST / file_name} "
                 f"{NORMALIZATION_PATH / file_name} ",
    dag=qclash_melanoma_dag,
)

paper_read >> mirna_seq_insertion >> blast_start
blast_end >> concat_blast >> normalization

