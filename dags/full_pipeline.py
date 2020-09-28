from datetime import timedelta

from airflow import DAG
from airflow.operators.dummy_operator import DummyOperator
from airflow.operators.subdag_operator import SubDagOperator

from consts.dag import DYNAMIC_FLY_GAMBIAE, FULL_PIPELINE_NAME, MIRBASE_DOWNALOD, default_args
from dags.dynamic_fly_gambiae import dynamic_fly_dag
from dags.mirbase_download import mirbase_download_dag

def subdag(parent_dag_name, child_dag_name):
    dag_subdag = DAG(
            dag_id='%s.%s' % (parent_dag_name, child_dag_name),
            default_args=default_args,
            # schedule_interval=timedelta(days=30),
    )
    return dag_subdag




dag = DAG(
    FULL_PIPELINE_NAME,
    default_args=default_args,
    description='full pipeline',
    # schedule_interval=timedelta(days=30),
)


pipeline_start = DummyOperator(
        task_id='pipeline_start',
        dag=dag,
        )

mirbase_download = SubDagOperator(
        task_id=MIRBASE_DOWNALOD,
        subdag=subdag(FULL_PIPELINE_NAME, MIRBASE_DOWNALOD),
        dag=dag,
        )

dynamic_fly = SubDagOperator(
        task_id=DYNAMIC_FLY_GAMBIAE,
        subdag=subdag(FULL_PIPELINE_NAME, DYNAMIC_FLY_GAMBIAE),
        dag=dag,
        )


pipeline_end = DummyOperator(
        task_id='pipeline_end',
        dag=dag,
        )



pipeline_start >> mirbase_download >> dynamic_fly >> pipeline_end
