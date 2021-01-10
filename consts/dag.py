from airflow.utils.dates import days_ago
from datetime import timedelta

FILE_NAMES = ['dynamic_fly',
              'cattle',
              'human_mapping',
              'darnell_mouse',
              'darnell_human',
              'unambiguous_human',
              'unambiguous_celegans',
              'unambiguous_mouse',
              'pairing_beyond']


# Dag Names
MIRBASE_DOWNALOD = "mirbase_downalod"
FULL_PIPELINE_NAME = "full_pipeline"
DYNAMIC_FLY_GAMBIAE = "Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae"
HUMAN_MAPPING = "Mapping_the_Human_miRNA_Interactome_by_CLASH_Reveals_Frequent_Noncanonical_Binding"
CATTLE = "Global_mapping_of_miRNA-target_interactions_in_cattle_Bos_taurus"
DARNELL = "Darnell_miRNA_targe_chimeras_reveal"
UNAMBIGUOUS = "Unambiguous_Identification_of_miRNA_Target_Site_Interactions"
PAIRING_BEYOND = "Pairing_beyond_the_Seed_Supports_MicroRNA_Targeting_Specificity"



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
