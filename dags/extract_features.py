# The DAG object; we'll need this to instantiate a DAG
from itertools import product
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.dummy_operator import DummyOperator
from airflow.operators.python_operator import PythonOperator

from consts.dag import FILE_NAMES, default_args
from consts.global_consts import DUPLEX_DICT
from utils.utils import split_file

extract_features_dag = DAG(
    "extract_features",
    default_args=default_args,
    description=f"extract_features",
    # schedule_interval=timedelta(days=30),
)



path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"
#
#
# def duplex_graph(method: str, dag: DAG):
#     duplex_start = DummyOperator(
#         task_id=f"{method}_duplex_start",
#         dag=dag)
#
#     duplex_done = DummyOperator(
#         task_id=f"{method}_duplex_done",
#         dag=dag)
#
#     for f in FILE_NAMES:
#         t = BashOperator(
#             task_id=f"{method}_{f}",
#             bash_command=f"python {step_path}duplex_step.py duplex {method} "
#                          f"{data_step_path}normalization_final/{f}.csv {data_step_path}duplex/{f}_{method}.csv",
#             dag=dag,
#         )
#         duplex_start >> t >> duplex_done
#
#     return duplex_start, duplex_done
#
# duplex_start = DummyOperator(
#     task_id=f"duplex_start",
#     dag=extract_features_dag)
#
# duplex_done = DummyOperator(
#     task_id=f"duplex_done",
#         dag=extract_features_dag)
#
# duplex_dags = []
# for duplex_method in DUPLEX_DICT.keys():
#     start, done = duplex_graph(duplex_method, extract_features_dag)
#
#     duplex_start >> start
#     done >> duplex_done
#

#################################################
# Feature extraction
#################################################
files_for_features = ["_".join(t) for t in product(FILE_NAMES, DUPLEX_DICT.keys())]
files_for_features.sort()

NUMBER_OF_CHUNKS = 5

split_start = DummyOperator(
    task_id=f"split_start",
        dag=extract_features_dag)


split_done = DummyOperator(
    task_id=f"split_done",
        dag=extract_features_dag)



for f in files_for_features:
    t = PythonOperator(
    task_id=f"split_{f}",
    python_callable=split_file,
    op_kwargs={'infile': f"{data_step_path}duplex/{f}.csv",
               'dir': Path(data_step_path) / "duplex"/ "split_files",
               'number_of_chunks': NUMBER_OF_CHUNKS},
    dag=extract_features_dag)

    split_start >> t >> split_done


features_start = DummyOperator(
    task_id=f"features_start",
        dag=extract_features_dag)

split_done >> features_start

features_done = DummyOperator(
    task_id=f"features_done",
        dag=extract_features_dag)

for f in files_for_features:
    for i in range(NUMBER_OF_CHUNKS):
        t = BashOperator(
            task_id=f"features_{f}_{i}",
            bash_command=f"python {step_path}feature_extraction.py feature-extraction "
                         f"{data_step_path}duplex/split_files/{f}{i}.csv {data_step_path}features/{f}_features{i}.csv",
            dag=extract_features_dag,
        )
        features_start >> t >> features_done


