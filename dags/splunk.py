# The DAG object; we'll need this to instantiate a DAG
from itertools import product
from pathlib import Path

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.dummy_operator import DummyOperator
from airflow.operators.python_operator import PythonOperator

from consts.dag import FILE_NAMES, SPLUNK_DICT, default_args
from consts.global_consts import DUPLEX_DICT
from utils.utils import split_file

splunk_dag = DAG(
    "splunk",
    default_args=default_args,
    description=f"prepare for splunk",
    # schedule_interval=timedelta(days=30),
)



path = "/home/local/BGU-USERS/benorgi/pipeline/"
step_path = path + "pipeline_steps/"
data_step_path = path + "data/pipeline_steps/"
feature_dir = Path(data_step_path)/ "features"
splunk_dir_full = Path(data_step_path)/ "splunk_full"
splunk_dir_lite = Path(data_step_path)/ "splunk_lite"

#################################################
# Feature extraction
#################################################
NUMBER_OF_CHUNKS = 5


start = DummyOperator(
    task_id=f"start",
        dag=splunk_dag)


done = DummyOperator(
    task_id=f"done",
        dag=splunk_dag)

for fname, (host, source) in SPLUNK_DICT.items():
    for duplex in DUPLEX_DICT.keys():
        for i in range(NUMBER_OF_CHUNKS):
            fin = feature_dir / f"{fname}_{duplex}_features{i}.csv"
            fout = splunk_dir_full / duplex / host / source / f"{fname}_{duplex}_features{i}.csv"
            fout_lite = splunk_dir_lite / duplex / host / source / f"{fname}_{duplex}_features{i}.csv"


            t = BashOperator(
                task_id=f"splunk{fin.stem}",
                bash_command=f"python {step_path}splunk.py splunk "
                             f"{fin} {fout}",
                dag=splunk_dag,
            )
            y = BashOperator(
                task_id=f"splunk-mirnaid-fix_{fin.stem}",
                bash_command=f"python {step_path}splunk.py mirnaid-fix "
                             f"{fout}",
                dag=splunk_dag,
            )
            z = BashOperator(
                task_id=f"splunk-lite_{fin.stem}",
                bash_command=f"python {step_path}splunk.py splunk-lite "
                             f"{fout} {fout_lite}",
                dag=splunk_dag,
            )


            start >> t >> y >> z >> done
