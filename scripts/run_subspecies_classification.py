#!/usr/bin/env python

import argparse
import json
import os
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

MAFFT_OUTPUT_F_NAME = "ref_MSA_with_query_seqs.fasta"
PPLACER_OUTPUT_F_NAME = "out.json"
GUPPY_OUTPUT_F_NAME = "out.sing.tre"
CLADINATOR_OUTPUT_F_NAME = "cladinator_results.tsv"

CLADE_DELIMITER = "'.+\{(.+)\}'"

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

  output_file = os.path.join(output_dir, job_data["output_file"] + ".tsv")

  #Find reference paths for the virus type
  virus_type = job_data["virus_type"]
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
  cladinator_cmd = ["cladinator", guppy_output, output_file]
  if "clade_mapping" in job_data:
    mapping_cmd = ["-m=%s" %(job_data["clade_mapping"]), "-S=%s" %(CLADE_DELIMITER)]
    cladinator_cmd.extend(mapping_cmd)
  try:
    subprocess.check_call(cladinator_cmd, shell=False)

    #Generate query dictionary to match with tre files in next step
    query_dict = {}
    with open(output_file, "r") as f:
      next(f) 
      previous_query = ""
      for line in f:
        split = line.split('\t')
        if previous_query != split[0]:
          previous_query = split[0]
          query_dict[split[0]] = split[0] + "\t" + split[1] + "\t" + split[2] + "\n"
    
    #Generate tre files for each query
    file_name_syntax = "%s.tre"
    with open(guppy_output, "r") as f:
      lines = f.readlines()
      for i in range(len(lines)):
        for key in query_dict:
          if key in lines[i]:
            file_name = file_name_syntax %(key) 
            break
        with open(os.path.join(output_dir, file_name), "w") as q:
          q.write(lines[i])

    #Create summary file
    result_file = os.path.join(output_dir, "result.tsv")
    with open(result_file, "w") as s:
      s.write("Query Identifier\tClade Classification\tLink\n")
      for key, value in query_dict.iteritems():
        s.write(value)

    #Remove initial output file
    os.remove(output_file)
  except Exception as e:
    print("Error running cladinator:\n %s" %(e))
    sys.exit(-1)
