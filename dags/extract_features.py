from datetime import timedelta
# The DAG object; we'll need this to instantiate a DAG
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.dummy_operator import DummyOperator

from consts.dag import FILE_NAMES, default_args


#instantiates a directed acyclic graph
from consts.global_consts import DUPLEX_DICT
from consts.pipeline_steps import CONCAT_BLAST, MIRNA_SEQ_PATH, NORMALIZATION_PATH, READ_PATH, REGION_PATH, SITE_PATH
from utils.blast import get_blast_graph

extract_features_dag = DAG(
    "extract_features",
    default_args=default_args,
    description=f"extract_features",
    # schedule_interval=timedelta(days=30),
)



path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"


def duplex_graph(method: str, dag: DAG):
    duplex_start = DummyOperator(
        task_id=f"{method}_duplex_start",
        dag=dag)

    duplex_done = DummyOperator(
        task_id=f"{method}_duplex_done",
        dag=dag)

    for f in FILE_NAMES:
        t = BashOperator(
            task_id=f"{method}_{f}",
            bash_command=f"python {step_path}duplex_step.py duplex {method} "
                         f"{data_step_path}normalization_final/{f}.csv {data_step_path}duplex/{f}_{method}.csv",
            dag=dag,
        )
        duplex_start >> t >> duplex_done

    return duplex_start, duplex_done

duplex_start = DummyOperator(
    task_id=f"duplex_start",
    dag=extract_features_dag)

duplex_done = DummyOperator(
    task_id=f"duplex_done",
        dag=extract_features_dag)

duplex_dags = []
for duplex_method in DUPLEX_DICT.keys():
    start, done = duplex_graph(duplex_method, extract_features_dag)

    duplex_start >> start
    done >> duplex_done



