#!/usr/bin/env python

import argparse
from datetime import datetime
import json
import re
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
tree_local = os.path.join("/home", "ac.mkuscuog", "git", "bvbrc_subspecies_classification", "lib", "ref-tree-alignment");
if os.path.exists(tree_deployed):
  ALIGNMENT_PATH = tree_deployed
elif os.path.exists(tree_dev):
  ALIGNMENT_PATH = tree_dev
else:
  ALIGNMENT_PATH = tree_local

report_deployed = os.path.join(top, "lib", "classification_report.html");
report_dev = os.path.join(top, "modules", "bvbrc_subspecies_classification", "lib", "classification_report.html");
report_local = os.path.join("/home", "ac.mkuscuog", "git", "bvbrc_subspecies_classification", "lib", "classification_report.html");
if os.path.exists(report_deployed):
  REPORT_TEMPLATE_PATH = report_deployed
elif os.path.exists(report_dev):
  REPORT_TEMPLATE_PATH = report_dev 
else:
  REPORT_TEMPLATE_PATH = report_local

rota_genotyper_deployed = os.path.join(top, "lib", "rota-a-genotyper");
rota_genotyper_dev = os.path.join(top, "modules", "bvbrc_subspecies_classification", "lib", "rota-a-genotyper");
rota_genotyper_local = os.path.join("/home", "ac.mkuscuog", "git", "bvbrc_subspecies_classification", "lib", "rota-a-genotyper");
if os.path.exists(rota_genotyper_deployed):
  ROTA_GENOTYPER_PATH = rota_genotyper_deployed
elif os.path.exists(rota_genotyper_dev):
  ROTA_GENOTYPER_PATH = rota_genotyper_dev
else:
  ROTA_GENOTYPER_PATH = rota_genotyper_local

if "P3_BASE_URL" in os.environ:
  BASE_URL = os.environ["P3_BASE_URL"]
else:
  BASE_URL = "https://www.bv-brc.org"

MAFFT_OUTPUT_F_NAME = "ref_MSA_with_query_seqs.fasta"
PPLACER_OUTPUT_F_NAME = "out.json"
GUPPY_OUTPUT_F_NAME = "out.sing.tre"
CLADINATOR_OUTPUT_F_NAME = "cladinator_results.tsv"

GENOTYPER_RESULT_F_NAME = "input.fasta.result"
GENOTYPER_ERROR_F_NAME = "input.fasta.err"

CLADE_DELIMITER = "(.+?)\|.+"
CLADE_DELIMITER_INFLUENZAH5 = ".+_\{(.+)\}"
REPORT_DATE = "<span>Report Date:</span> %s"
TABLE_HEADER_C = "<th class=\"dgrid-cell dgrid-cell-padding\">Query Identifier</th><th class=\"dgrid-cell dgrid-cell-padding\">Clade Classification</th><th class=\"dgrid-cell dgrid-cell-padding\">Tree Link</th>"
TABLE_HEADER_R_RESULT = "<th class=\"dgrid-cell dgrid-cell-padding\">Input FASTA unique ID</th><th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Segment number</th><th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Genotype</th><th class=\"dgrid-cell dgrid-cell-padding\">Best hit accession</th><th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:13%\">Query coverage %</th><th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Ident %</th><th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">E Value</th>"
TABLE_HEADER_R_ERR = "<th class=\"dgrid-cell dgrid-cell-padding\">Input FASTA unique ID</th><th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Segment number</th><th class=\"dgrid-cell dgrid-cell-padding\">Error Description</th>"
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
    rota_genotyper_cmd = ["ss-rotaA-genotyper", input_file]
    try:
      subprocess.check_call(rota_genotyper_cmd, shell=False)

      #Create html summary file from result and error files
      report_file = os.path.join(output_dir, job_data["output_file"] + "_classification_report.html")
      with open(REPORT_TEMPLATE_PATH, 'r') as f:
        html_data = f.readlines()

      html_data[60] = REPORT_DATE %(datetime.now().strftime("%B %d, %Y %H:%M:%S"))

      result_rows = ""
      result_file = os.path.join(output_dir, GENOTYPER_RESULT_F_NAME)
      with open(result_file, "r") as f:
        lines = f.readlines()[1:]
        if len(lines) > 0:
          html_data[68] = TABLE_HEADER_R_RESULT
          for line in lines:
            split = line.split('\t')
            result_rows += "<tr>"
            for data in split:
              result_rows += TABLE_ROW.replace("%{data}", data)
            result_rows += "</tr>"
      html_data[70] = result_rows

      error_rows = ""
      error_file = os.path.join(output_dir, GENOTYPER_ERROR_F_NAME)
      with open(error_file, "r") as f:
        lines = f.readlines()[1:]
        if len(lines) > 0:
          html_data[76] = TABLE_HEADER_R_ERR
          for line in lines:
            split = line.split('\t')
            error_rows += "<tr>"
            for data in split:
              error_rows += TABLE_ROW.replace("%{data}", data)
            error_rows += "</tr>"
      html_data[78] = error_rows

      with open(report_file, 'w') as f:
        f.writelines(html_data)
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
      if file.endswith(".aln") or file.endswith(".fasta"):
        reference_mfa_file = os.path.join(ref_folder_path, file)
      elif file.endswith(".nh"):
        reference_tree_file = os.path.join(ref_folder_path, file)
      elif file.find("info.") != -1:
        stats_file = os.path.join(ref_folder_path, file)
      elif file.endswith(".tsv"):
        mapping_file = os.path.join(ref_folder_path, file)
  
    #Get reference MSA file from workspace if selected
    if "ref_msa_fasta" in job_data:
      reference_mfa_file = os.path.join(output_dir, "ref_mfa.fasta")
      try:
        fetch_msa_cmd = ["p3-cp", "ws:%s" %(job_data["ref_msa_fasta"]), reference_mfa_file]
        subprocess.check_call(fetch_msa_cmd, shell=False)
      except Exception as e:
        print("Error copying mfa fasta file from workspace:\n %s" %(e))
        sys.exit(-1)
        
    #Aligning of the query sequence(s) to the reference multiple sequence alignment
    mafft_cmd = ["mafft", "--keeplength", "--add", input_file, reference_mfa_file]
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

    cladinator_cmd = ["cladinator", guppy_output, output_file]
    if virus_type == "INFLUENZAH5":
      cladinator_cmd.insert(1, "-S=%s" %(CLADE_DELIMITER_INFLUENZAH5))
    if virus_type == "INFLUENZAH5" or virus_type == "SWINEH1" or virus_type == "SWINEH3" or virus_type == "SWINEH1US":
      cladinator_cmd.insert(1, "-m=%s" %(mapping_file))
    try:
      subprocess.check_call(cladinator_cmd, shell=False)
    except Exception as e:
      print("Error running cladinator:\n %s" %(e))
      sys.exit(-1)
 
    try:
      #Generate query dictionary to match with tre files in next step
      query_dict = {}
      with open(output_file, "r") as f:
        next(f) 
        for line in f:
          split = line.split('\t')
          query = split[0]
          if split[1] == "Matching Clades" and (not query_dict.has_key(query) or query_dict[query] == "?"):
            query_dict[query] = split[2]
          elif split[1] == "Matching Down-tree Bracketing Clades" and query_dict[query] == "?":
            #Use down-tree classification value if matching clade is ?
            if split[2] == "?":
              query_dict[query] = "Sequence cannot be classified based on the reference tree"
            else:
              query_dict[query] = split[2] + "-like"
      
      if virus_type == "INFLUENZAH5" or virus_type == "SWINEH1" or virus_type == "SWINEH3" or virus_type == "SWINEH1US":
        decorator_output = os.path.join(output_dir, "outtree.tre")
        decorator_cmd = ["decorator", "-f=n", "-nh", "INPUT_TRE_FILE", mapping_file, decorator_output]
      #Generate tre files for each query
      file_name_syntax = "%s.tre"
      with open(guppy_output, "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
          for key in query_dict:
            if key in lines[i]:
              key = re.sub("[^a-zA-Z0-9 \n\.]", "_", key) 
              if len(key) > 250:
                print('File name is too long. Truncating it')
                key = key[0:200] + '_'
              file_name = file_name_syntax %(key)
              break
          if "file_name" in dir():
            file_path = os.path.join(output_dir, file_name)
            with open(file_path, "w") as q:
              q.write(lines[i])

            #Update tre file for influenza to display labels in phylogenetic tree
            if virus_type == "INFLUENZAH5" or virus_type == "SWINEH1" or virus_type == "SWINEH3" or virus_type == "SWINEH1US":
              try:
                decorator_cmd[3] = file_path
                subprocess.check_call(decorator_cmd, shell=False)
                shutil.move(decorator_output, file_path)
              except Exception as e:
                print("Error running decorator for %s:\n %s" %(file_name, e))
                sys.exit(-1)
  
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
        rows += TABLE_ROW.replace("%{data}", value)
        initial_value = value.startswith("Sequence") and "" or re.sub("\-like$", "", value) 
        #TODO: handle file name truncating in a function during dict creation
        if len(key) > 250:
          key = key[0:200] + '_'
        rows += TABLE_ROW.replace("%{data}", TREE_LINK %(BASE_URL, job_data["output_path"], job_data["output_file"], re.sub("[^a-zA-Z0-9 \n\.]", "_", key), initial_value))
        rows += "</tr>"
  
      html_data[60] = REPORT_DATE %(datetime.now().strftime("%B %d, %Y %H:%M:%S"))
      html_data[68] = TABLE_HEADER_C
      html_data[70] = rows
      with open(report_file, 'w') as f:
        f.writelines(html_data)
  
      #Remove initial result 
      #os.remove(output_file)
    except Exception as e:
      print("Error generating result files:\n %s" %(e))
      sys.exit(-1)

  try:
    #Move data files to details folder
    details_folder = os.path.join(output_dir, "details")
    if not os.path.exists(details_folder):
      os.mkdir(details_folder)

    files = os.listdir(output_dir)
    for f in files:
      if not f.endswith("html"):
        shutil.move(os.path.join(output_dir, f), details_folder)
  except Exception as e:
    print("Error moving data files to details folder:\n %s" %(e))
    sys.exit(-1)
