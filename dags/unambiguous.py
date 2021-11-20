from datetime import timedelta
# The DAG object; we'll need this to instantiate a DAG
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from consts.dag import UNAMBIGUOUS, default_args


#instantiates a directed acyclic graph
from consts.pipeline_steps import CONCAT_BLAST, MIRNA_SEQ_PATH, NORMALIZATION_PATH, READ_PATH, REGION_PATH, SITE_PATH
from utils.blast import get_blast_graph

unambiguous_dag = DAG(
    UNAMBIGUOUS,
    default_args=default_args,
    description=f"Parse the paper: {UNAMBIGUOUS}",
    # schedule_interval=timedelta(days=30),
)



path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"

genome_path = {"mouse": f"{path}data/genome/mm9",
               "human": f"{path}data/genome/hg18",
               }

def get_unambiguous_graph(organism: str):
    file_name = f"unambiguous_{organism}.csv"


    paper_read = BashOperator(
        task_id=f"paper_read_{organism}",
        bash_command=f"python {path}papers/{UNAMBIGUOUS}/reader.py {organism} {file_name}",
        dag=unambiguous_dag,
    )

    mirna_seq_insertion = BashOperator(
        task_id=f'miRNA_id_fix_{organism}',
        bash_command="python " + step_path + "mirna_seq_insertion.py mirnaid-fix " + file_name,
        dag=unambiguous_dag,
    )


    # rna_site_insertion = BashOperator(
    #     task_id=f"rna_site_insertion_{organism}",
    #     bash_command=f"python {step_path}rna_site_insertion.py insert-site-from-chromosome "
    #                  f"{MIRNA_SEQ_PATH}/{file_name} "
    #                  f"{SITE_PATH}/{file_name} "
    #                  f"{genome_path[organism]} ",
    #     dag=unambiguous_dag,
    # )
    #
    blast_start, blast_tasks, blast_end = get_blast_graph(unambiguous_dag,
                                                          (MIRNA_SEQ_PATH / file_name),
                                                          organism)


    concat_blast = BashOperator(
        task_id=f"concat_blast_{organism}",
        bash_command=f"python {step_path}concat_blast_result.py concat-blast-result "
                     f"{REGION_PATH}/ "
                     f"{Path(file_name).stem} "
                     f"{MIRNA_SEQ_PATH / file_name} "
                     f"{CONCAT_BLAST / file_name} ",
        dag=unambiguous_dag,
    )

    normalization = BashOperator(
        task_id=f"normalization_{organism}",
        bash_command=f"python {step_path}normalization_final_step.py finalize "
                     f"{CONCAT_BLAST / file_name} "
                     f"{NORMALIZATION_PATH / file_name} ",
        dag=unambiguous_dag,
    )

    paper_read >> mirna_seq_insertion >> blast_start
    blast_end >> concat_blast >> normalization


for organism in ["celegans", "human", "mouse"]:
    get_unambiguous_graph(organism)

