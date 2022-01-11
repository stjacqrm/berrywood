#!/usr/bin/env nextflow

//Description: Analyzing SC2 genomes
//Author: Rachael St. Jacques
//email: rachael.stjacques@dgs.virginia.gov

//starting parameters
//params.report = workflow.launchDir + '/report/report_template.Rmd'
//params.bash = workflow.launchDir + '/bash/create_berrywood_report.sh'

//setup channel to read in the fasta files
Channel
    .fromPath("${params.assemblies}/*.fasta")
    .ifEmpty { exit 1, "Cannot find any fasta assemblies in ${params.assemblies}" }
    .set { assemblies}

Channel
    .fromPath("$baseDir/report/report_template.Rmd")
    .set { report }

Channel
    .fromPath("$baseDir/bash/create_berrywood_report.sh")
    .set { bash }


//Step 1: collect the SC2 assembled_genomes
process collect_fasta {
  publishDir "${params.outdir}/logs/multi_fasta/", mode: 'copy'


  input:
  file(assembly) from assemblies.collect()

  output:
  file("multi.fasta") into annotate_multi_fasta,split_multi,nextclade_genomes,pangolin_genomes

  script:
  """
  cat *.fasta > multi.fasta;
  """

}

//Step 2: annotate SC2 genomes
process annotate {
  publishDir "${params.outdir}/results/vadr_results/annotations", mode: 'copy', pattern:'*'

  input:
  file(assembly) from annotate_multi_fasta

  output:
  file("*vadr_out.vadr.sqa") into convert_tsv

  script:
  """
  v-annotate.pl --split --cpu ${task.cpus} --glsearch -s -r -f --nomisc --mkey sarscov2 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models ${assembly} vadr_out
  cd vadr_out && cp * ..
  """

}

//Step 3: convert sqa to a usable csv file
process convert_sqa {
  publishDir "${params.outdir}/results/vadr_results/csvs", mode: 'copy' , pattern:"*"

  input:
  file(tsv) from convert_tsv

  output:
  file("*vadr_results.csv") into vadr_results, compile_vadr,render_vadr

  script:
  """
  awk -v OFS='\t' '{\$1=\$1; print}' vadr_out.vadr.sqa > pass_fail_1.tsv
  sed -i 1d pass_fail_1.tsv
  sed -i 2d pass_fail_1.tsv
  awk -F'\t' '{sub(","," ",\$14)}1' OFS='\t' pass_fail_1.tsv > pass_fail_2.tsv
  awk -F'\t' '{sub(","," ",\$14)}1' OFS='\t' pass_fail_2.tsv > pass_fail_3.tsv
  sed 's/\t/,/g' pass_fail_3.tsv > pass_fail_3.csv
  cut -d, -f1 --complement pass_fail_3.csv > vadr_results.csv
  sed -i '1s/.*/isolate,length,status,annotated,model,subgenera,species,num_features_annotated,features_not_ann,num_feat_5,num_feat_3,feat_alerts,seq_alerts/' vadr_results.csv
  """
}

//Step 4: nextclade for mutations
process nextlclade_mutations {
  publishDir "${params.outdir}/results/nextclade_results", mode: 'copy' , pattern:"*_results.csv"
  publishDir "${params.outdir}/results/nextclade_results", mode: 'copy', pattern:"*_complete_version.csv"

  input:
  file(assembly) from nextclade_genomes

  output:
  file("*_results.csv") into nextclade_results
  file("*complete_version.csv") into nextclade_version

  script:
  """
  nextclade --input-fasta "${assembly}" --output-csv "nextclade_results.csv"
  nextclade --version > nextclade_version.csv
  echo "version" >> nextclade_complete_version.csv
  cat nextclade_version.csv >> nextclade_complete_version.csv
  """
}

//Step 5: pangolin for lineages
process lineage {
  publishDir "${params.outdir}/results/pangolin_results", mode: 'copy' , pattern:"*_report.csv"

  input:
  file(assembly) from pangolin_genomes

  output:
  file("*report.csv") into lineage_results
  file("*complete_lineage_version.csv") into pangolin_version

  script:
  """
  pangolin ${assembly} --outfile lineage_report.csv
  pangolin --version > pango_version.csv
  echo "version" >> complete_lineage_version.csv
  cat pango_version.csv >> complete_lineage_version.csv
  """
}

//Step 6: build the report
process build_report {
  beforeScript 'ulimit -s unlimited'
  tag "$name"
  publishDir "${params.outdir}", mode: 'copy', pattern:"*.csv"

  input:
  file(nextclade_clades) from nextclade_results
  file(pangolin_lineage) from lineage_results
  file(nextclade_version) from nextclade_version
  file(pangolin_version) from pangolin_version
  file(vadr_results) from vadr_results

  output:
  file("*.csv") into render_berrywood

  script:
"""
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
from datetime import datetime
today = datetime.today()
today = today.strftime("%m%d%y")
nextclade_df = pd.read_csv('nextclade_results.csv',sep=';')
nextclade_df = nextclade_df[['seqName','clade','insertions','deletions','substitutions','aaSubstitutions']]
nextclade_df.columns = ['sample','clade','insertions','deletions','nucleotide_subs','amino_acid_subs']
nextclade_df['sample'] = pd.Series(nextclade_df['sample'], dtype="string")
pangolin_df = pd.read_csv('lineage_report.csv',sep=',')
pangolin_df = pangolin_df[['taxon','lineage']]
pangolin_df.columns = ['sample','pangolin_lineage']
pangolin_df['sample'] = pd.Series(pangolin_df['sample'], dtype="string")
vadr_df = pd.read_csv('vadr_results.csv',sep=',')
vadr_df = vadr_df[['isolate','status']]
vadr_df.columns = ['sample','vadr_status']
vadr_df['sample'] = pd.Series(vadr_df['sample'], dtype="string")
nextclade_version_df = pd.read_csv('nextclade_complete_version.csv',sep=',')
nextclade_version_df = nextclade_version_df[['version']]
nextclade_version_df.columns = ['nextclade_version']
pangolin_version_df = pd.read_csv('complete_lineage_version.csv',sep=',')
pangolin_version_df = pangolin_version_df[['version']]
pangolin_version_df.columns = ['pangolin_version']
dataframes = [nextclade_df,pangolin_df,vadr_df]
total = pd.merge(nextclade_df,pangolin_df, on='sample')
total = pd.merge(total,vadr_df,on='sample')
total['pangolin_version'] = pangolin_version_df._get_value(0,'pangolin_version')
total['nextclade_version'] = nextclade_version_df._get_value(0,'nextclade_version')
total = total[['sample', 'vadr_status', 'pangolin_lineage','pangolin_version','clade','nextclade_version','insertions','deletions','nucleotide_subs','amino_acid_subs']]
title = ["berrywood_", today,".csv"]
pd.DataFrame.to_csv(total, "".join(title), sep=';',index=False)
"""
}

//Step 7: generate the pdf report
process render{
  publishDir "${params.outdir}/", mode: 'copy', pattern:'Berrywood-report.pdf'
  publishDir "${params.outdir}/", mode: 'copy', pattern:"berrywood*"
  publishDir "${params.outdir}/results/berrywood_results/csvs", mode: 'copy' , pattern:"*csv"
  beforeScript 'ulimit -s unlimited'

  input:
  file(berrywood) from render_berrywood
  file(vadr) from render_vadr
  file(rmd) from report
  file(bash) from bash

  output:
  file('*.csv')
  file "Berrywood-report.pdf"

  shell:
  """
  cp ${rmd} ./report_template.Rnw
  cp ${bash} ./create_report.sh
  chmod +x create_report.sh
  bash create_report.sh -p "Berrywood-report" -t "${params.title}" -T report_template.Rnw -o . -a ${berrywood} -v ${vadr}
  """

}
