# Define variables
datasets = ["hub_01", "hub_02"]
output_files = [f"QC_dashboard_{dataset}.html" for dataset in datasets]

# Rule to render QC_dashboard.qmd for each dataset
rule all:
    input:
        output_files

rule render_qmd:
    input:
        "QC_dashboard.qmd"
    output:
        "QC_dashboard_{dataset}.html"
    conda: 
        "envs/QC_dashboard_environment.yaml"
    shell:
        """
        Rscript -e "library(quarto);
        quarto_render(
          input = '{input}',
          output_file = '{output}',
          execute_params = list(dataset_name = '{wildcards.dataset}')
        )"
        """
