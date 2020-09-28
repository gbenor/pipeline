from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from datetime import timedelta

from consts.dag import MIRBASE_DOWNALOD, default_args


#instantiates a directed acyclic graph
mirbase_download_dag = DAG(
    MIRBASE_DOWNALOD,
    default_args=default_args,
    description='download mirbase mature.fa of several versions and create mirna sequences database',
    # schedule_interval=timedelta(days=30),
)


script = "python /home/local/BGU-USERS/benorgi/pipeline/mirna_utils/mirbase.py "


download_mirbase = BashOperator(
    task_id='download_from_mirbase',
    bash_command=script + 'mirbase-download 17',
    dag=mirbase_download_dag,
)


read_fasta_files = BashOperator(
    task_id='read_and_parse_fasta_files',
    bash_command=script + 'read-fasta-files 17',
    dag=mirbase_download_dag,
)


read_aga_fasta = BashOperator(
    task_id='read_aga_fasta',
    bash_command="python /home/local/BGU-USERS/benorgi/pipeline/mirna_utils/gambiae.py",
    dag=mirbase_download_dag,
)


postprocessing = BashOperator(
    task_id='post-processing',
    bash_command=script + 'postprocess',
    dag=mirbase_download_dag,
)


#sets the ordering of the DAG. The >> directs the 2nd task to run after the 1st task. This means that
#download images runs first, then train, then serve.

download_mirbase >> read_fasta_files >> read_aga_fasta >> postprocessing


