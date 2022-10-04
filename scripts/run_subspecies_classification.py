#!/usr/bin/env python

import argparse
from datetime import datetime
import json
import os
import shutil
import subprocess
import sys

#
# Determine path to our alignments.
#
top = os.getenv("KB_TOP")
tree_deployed = os.path.join(top, "lib", "ref-tree-alignment");
tree_dev = os.path.join(top, "modules", "bvbrc_subspecies_classification", "lib", "ref-tree-alignment");
if os.path.exists(tree_deployed):
  ALIGNMENT_PATH = tree_deployed
else:
  ALIGNMENT_PATH = tree_dev

report_deployed = os.path.join(top, "lib", "classification_report.html");
report_dev = os.path.join(top, "modules", "bvbrc_subspecies_classification", "lib", "classification_report.html");
if os.path.exists(tree_deployed):
  REPORT_TEMPLATE_PATH = report_deployed
else:
  REPORT_TEMPLATE_PATH = report_dev 

rota_genotyper_deployed = os.path.join(top, "lib", "rota-a-genotyper");
rota_genotyper_dev = os.path.join(top, "modules", "bvbrc_subspecies_classification", "lib", "rota-a-genotyper");
if os.path.exists(rota_genotyper_deployed):
  ROTA_GENOTYPER_PATH = rota_genotyper_deployed
else:
  ROTA_GENOTYPER_PATH = rota_genotyper_dev

GENOTYPER_JAR_NAME = "StandAloneRtvAGenotyper.jar"
GENOTYPER_CONFIG_NAME = "rotaAGenotyper.config"

if "P3_BASE_URL" in os.environ:
  BASE_URL = os.environ["P3_BASE_URL"]
else:
  BASE_URL = "https://www.bv-brc.org"

MAFFT_OUTPUT_F_NAME = "ref_MSA_with_query_seqs.fasta"
PPLACER_OUTPUT_F_NAME = "out.json"
GUPPY_OUTPUT_F_NAME = "out.sing.tre"
CLADINATOR_OUTPUT_F_NAME = "cladinator_results.tsv"

CLADE_DELIMITER = "'.+\{(.+)\}'"
REPORT_DATE = "<span>Report Date:</span> %s"
TABLE_ROW = "<td class=\"dgrid-cell dgrid-cell-padding\">%{data}</td>"
TREE_LINK = "<a href=\"%s/view/PhylogeneticTree2/?wsTreeFile=%s/.%s/details/%s.tre&fileType=nwk&isClassification=1&initialValue=%s\" target=\"_blank\">VIEW TREE</a>"

if __name__ == "__main__" :
  parser = argparse.ArgumentParser(description="SubSpecies Classification Script")
  parser.add_argument("-j", "--jfile", help="json file for job", required=True)
  parser.add_argument("-o", "--output", help="Output directory. defaults to current directory", required=False, default=".")

  args = parser.parse_args()

  #Load job data
  job_data = None
  try:
    with open(args.jfile, "r") as j:
      job_data = json.load(j)
  except Exception as e:
    print("Error in opening job file:\n %s" %(e))
    sys.exit(-1)

  if not job_data:
    print("job_data is null")
    sys.exit(-1) 
  print(job_data)

  #Setup output directory
  output_dir = args.output
  output_dir = os.path.abspath(output_dir)
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  os.chdir(output_dir)

  output_file = os.path.join(output_dir, job_data["output_file"] + ".txt")

  #Define input file
  input_file = os.path.join(output_dir, "input.fasta")
  if job_data["input_source"] == "fasta_file":
    #Fetch input file from workspace
    try:
      fetch_fasta_cmd = ["p3-cp", "ws:%s" %(job_data["input_fasta_file"]), input_file]
      subprocess.check_call(fetch_fasta_cmd, shell=False)
    except Exception as e:
      print("Error copying fasta file from workspace:\n %s" %(e))
      sys.exit(-1)
  elif job_data["input_source"] == "fasta_data":
    #Copy user data to input file
    try:
      with open(input_file, "w+") as input:
        input.write(job_data["input_fasta_data"])
    except Exception as e:
      print("Error copying fasta data to input file:\n %s" %(e))
      sys.exit(-1)

  if os.path.getsize(input_file) == 0:
    print("Input fasta file is empty")
    sys.exit(-1)

  virus_type = job_data["virus_type"]
  if virus_type == 'ROTAA':
    rota_genotyper_cmd = ["java", "-jar", os.path.join(ROTA_GENOTYPER_PATH, GENOTYPER_JAR_NAME), os.path.join(ROTA_GENOTYPER_PATH, GENOTYPER_CONFIG_NAME), input_file]
    try:
      subprocess.check_call(rota_genotyper_cmd, shell=False)
    except Exception as e:
      print("Error running rotavirus a genotyper:\n %s" %(e))
      sys.exit(-1)
  else:
    #Find reference paths for the virus type
    ref_folder_path = os.path.join(ALIGNMENT_PATH, virus_type)
    if not os.path.exists(ref_folder_path):
      print("Cannot find reference folder path for virus type %s" %(virus_type))
      sys.exit(-1)
  
    reference_mfa_file = None
    reference_tree_file = None
    stats_file = None
    for file in os.listdir(ref_folder_path):
      if file.endswith(".aln"):
        reference_mfa_file = os.path.join(ref_folder_path, file)
      elif file.endswith(".nh"):
        reference_tree_file = os.path.join(ref_folder_path, file)
      elif file.find("info.") != -1:
        stats_file = os.path.join(ref_folder_path, file)
  
    #Get reference MSA file from workspace if selected
    if "ref_msa_fasta" in job_data:
      reference_mfa_file = os.path.join(output_dir, "ref_mfa.fasta")
      try:
        fetch_msa_cmd = ["p3-cp", "ws:%s" %(job_data["ref_msa_fasta"]), reference_mfa_file]
        subprocess.check_call(fetch_msa_cmd, shell=False)
      except Exception as e:
        print("Error copying mfa fasta file from workspace:\n %s" %(e))
        sys.exit(-1)
  
    #Define input file
    input_file = os.path.join(output_dir, "input.fasta")
    if job_data["input_source"] == "fasta_file":
      #Fetch input file from workspace
      try:
        fetch_fasta_cmd = ["p3-cp", "ws:%s" %(job_data["input_fasta_file"]), input_file] 
        subprocess.check_call(fetch_fasta_cmd, shell=False)
      except Exception as e:
        print("Error copying fasta file from workspace:\n %s" %(e))
        sys.exit(-1)
    elif job_data["input_source"] == "fasta_data":
      #Copy user data to input file
      try:
        with open(input_file, "w+") as input:
          input.write(job_data["input_fasta_data"])
      except Exception as e:
        print("Error copying fasta data to input file:\n %s" %(e))
        sys.exit(-1)  
  
    if os.path.getsize(input_file) == 0:
      print("Input fasta file is empty")
      sys.exit(-1)
        
    #Aligning of the query sequence(s) to the reference multiple sequence alignment
    mafft_cmd = ["mafft", "--add", input_file, reference_mfa_file]
    mafft_output = os.path.join(output_dir, MAFFT_OUTPUT_F_NAME)
    try:
      with open(mafft_output, "w+") as o:
        subprocess.check_call(mafft_cmd, shell=False, stdout=o)
    except Exception as e:
      print("Error running mafft:\n %s" %(e))
      sys.exit(-1)
  
    #Pplacer
    pplacer_output = os.path.join(output_dir, PPLACER_OUTPUT_F_NAME)
    pplacer_cmd = ["pplacer", "-m", "GTR", "-t", reference_tree_file, "-s", stats_file, mafft_output, "-o", pplacer_output] 
    try:
      subprocess.check_call(pplacer_cmd, shell=False)
    except Exception as e:
      print("Error running pplacer:\n %s" %(e))
      sys.exit(-1)
  
    #Guppy to generate a tree for every query sequence
    guppy_cmd = ["guppy", "sing", pplacer_output]
    try:
      subprocess.check_call(guppy_cmd, shell=False)
    except Exception as e:
      print("Error running guppy:\n %s" %(e))
      sys.exit(-1)
  
    #Cladinator
    guppy_output = os.path.join(output_dir, GUPPY_OUTPUT_F_NAME)
    cladinator_output = os.path.join(output_dir, CLADINATOR_OUTPUT_F_NAME)
    cladinator_cmd = ["cladinator", "-S=(.+?)\|.+", guppy_output, output_file]
    #if "clade_mapping" in job_data:
    #  mapping_cmd = ["-m=%s" %(job_data["clade_mapping"]), "-S=%s" %(CLADE_DELIMITER)]
    #  cladinator_cmd.extend(mapping_cmd)
    try:
      subprocess.check_call(cladinator_cmd, shell=False)
  
      #Generate query dictionary to match with tre files in next step
      query_dict = {}
      with open(output_file, "r") as f:
        next(f) 
        for line in f:
          split = line.split('\t')
          query = split[0]
          if split[1] == "Matching Clades":
            query_dict[query] = split[2]
          elif split[1] == "Matching Down-tree Bracketing Clades" and query_dict[query] == "?":
            #Use down-tree classification value if matching clade is ?
            if split[2] == "?":
              query_dict[query] = "Sequence cannot be classified based on the reference tree"
            else:
              query_dict[query] = split[2]
      
      #Generate tre files for each query
      file_name_syntax = "%s.tre"
      with open(guppy_output, "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
          for key in query_dict:
            if key in lines[i]:
              file_name = file_name_syntax %(key) 
              break
          if "file_name" in dir():
            with open(os.path.join(output_dir, file_name), "w") as q:
              q.write(lines[i])
  
      #Create summary file
      result_file = os.path.join(output_dir, "result.tsv")
      with open(result_file, "w") as s:
        s.write("Query Identifier\tClade Classification\n")
        for key, value in query_dict.iteritems():
          s.write(key + "\t" + value + "\n")
  
      #Create html summary file
      report_file = os.path.join(output_dir, job_data["output_file"] + "_classification_report.html")
      with open(REPORT_TEMPLATE_PATH, 'r') as f:
        html_data = f.readlines()
  
      rows = ""
      for key, value in query_dict.iteritems():
        rows += "<tr>"
        rows += TABLE_ROW.replace("%{data}", key)
        rows += TABLE_ROW.replace("%{data}", value + "-like")
        initial_value = value.startswith("Sequence") and "" or value
        rows += TABLE_ROW.replace("%{data}", TREE_LINK %(BASE_URL, job_data["output_path"], job_data["output_file"], key, initial_value))
        rows += "</tr>"
  
      html_data[60] = REPORT_DATE %(datetime.now().strftime("%B %d, %Y %H:%M:%S"))
      html_data[72] = rows
      with open(report_file, 'w') as f:
        f.writelines(html_data)
  
      #Remove initial result 
      #os.remove(output_file)
  
      #Move data files to details folde
      details_folder = os.path.join(output_dir, "details")
      if not os.path.exists(details_folder):
        os.mkdir(details_folder)
  
      files = os.listdir(output_dir)
      for f in files:
        if not f.endswith("html"):
          shutil.move(os.path.join(output_dir, f), details_folder)
    except Exception as e:
      print("Error running cladinator:\n %s" %(e))
      sys.exit(-1)
