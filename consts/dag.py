from airflow.utils.dates import days_ago
from datetime import timedelta

# Dag Names
MIRBASE_DOWNALOD = "mirbase_downalod"
FULL_PIPELINE_NAME = "full_pipeline"
DYNAMIC_FLY_GAMBIAE = "Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae"





default_args = {
    'owner': 'Gilad Ben Or',
    'depends_on_past': False,
    'start_date': days_ago(31),
    'email': ['benorgi@post.bgu.ac.il'],
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=1),
    "schedule_interval" : 'None'
    # 'queue': 'bash_queue',
    # 'pool': 'backfill',
    # 'priority_weight': 10,
    # 'end_date': datetime(2016, 1, 1),
    # 'wait_for_downstream': False,
    # 'dag': dag,
    # 'sla': timedelta(hours=2),
    # 'execution_timeout': timedelta(seconds=300),
    # 'on_failure_callback': some_function,
    # 'on_success_callback': some_other_function,
    # 'on_retry_callback': another_function,
    # 'sla_miss_callback': yet_another_function,
    # 'trigger_rule': 'all_success'
}
