
# Pipeline: A Computational Framework for miRNA–Target Interaction Analysis

## Overview

Pipeline is a GitHub project leveraging Apache Airflow for managing compute-intensive tasks in parsing miRNA interaction datasets. Its primary goal is to create duplexes and generate both positive and negative samples for comprehensive analysis. This project requires an SQL server for running Airflow.

The foundation of Pipeline's methodology is inspired by the research presented in the paper "Comprehensive machine-learning-based analysis of microRNA–target interactions reveals variable transferability of interaction rules across species" published in BMC Bioinformatics. [Link to the paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04164-x#Sec23).

## Methodology

The project employs a sophisticated approach for data processing and interaction analysis, as outlined in the referenced paper:

1. **Data Acquisition**: Collection of high-throughput chimeric miRNA–target datasets from various species.
2. **Standardization of Data**: Conversion of datasets into a uniform format, incorporating metadata, miRNA names and sequences, target site sequences, and more.
3. **miRNA Sequence Retrieval**: Extraction of missing miRNA sequences based on their names from miRBase.
4. **Target Sequence Extraction**: Identification of target sequences based on genomic coordinates, focusing on the most functional sites within the 3’UTRs.
5. **Interaction Analysis**: Utilization of the ViennaRNA suite for calculating interaction duplexes, classifying them based on seed type, and distinguishing between positive and negative interactions.
6. **Negative Interaction Generation**: Creation of negative interactions through a synthetic method to balance the dataset.
7. **miRNA Distribution Calculation**: Analysis of miRNA sequence occurrence and cumulative distribution function generation.
8. **Feature Representation**: Utilization of 490 expert-designed features classified into categories and subcategories for representing miRNA–target interactions.

## Setup and Installation

1. **Requirements**:
    - Apache Airflow
    - SQL Server
    - Additional Python dependencies as listed in `requirements.txt`.

2. **Installation**:
    - Clone the repository: `git clone [repository URL]`.
    - Set up an SQL server for Airflow metadata.
    - Configure Apache Airflow to interface with the SQL server.
    - Install required Python packages: `pip install -r requirements.txt`.

3. **Configuration**:
    - Follow the configuration instructions in `config.json` to set up data paths and processing parameters.

4. **Usage**:
    - Use the Airflow web server to manage and monitor your workflows.
    - Detailed instructions on running specific pipelines are provided in the `docs` folder.

## Contribution Guidelines

Contributors are welcome to propose enhancements, report bugs, and suggest improvements via pull requests and issues. Please refer to the `CONTRIBUTING.md` file for detailed guidelines.

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## Acknowledgments

This work is a computational interpretation and implementation of methods described in the referenced paper. We acknowledge the authors and contributors of the original research for their groundbreaking work in the field of miRNA–target interaction analysis.
