#!/usr/bin/env python

import argparse
from datetime import datetime
import json
import re
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

# ===========================
# Constants
# ===========================

MAFFT_OUTPUT_F_NAME = "ref_MSA_with_query_seqs.fasta"
PPLACER_OUTPUT_F_NAME = "out.json"
GUPPY_OUTPUT_F_NAME = "out.sing.tre"
CLADINATOR_OUTPUT_F_NAME = "cladinator_results.tsv"
GENOTYPER_RESULT_F_NAME = "input.fasta.result"
GENOTYPER_ERROR_F_NAME = "input.fasta.err"

CLADE_DELIMITER = "([A-Za-z0-9._]+)\|.+"
CLADE_DELIMITER_INFLUENZAH5 = ".+_\{(.+)\}"
CLADE_DELIMITER_MPOX = "([A-Za-z0-9.]+)\|.+"
CLADE_DELIMITER_DENGUE = "(\d+).+"
# Add default CLADE_DELIMITER for these virus types only
CLADE_DELIMITER_VIRUS_TYPES = [
    'BOVDIARRHEA1', 'HCV', 'JAPENCEPH', 'MURRAY',
    'TKBENCEPH', 'YELLOWFEVER', 'ZIKA', 'STLOUIS', 'WESTNILE'
]

BAD_CHARS_REGEX = r"[\[\'\"(),;_|:\]]"

REPORT_DATE = "<span>Report Date:</span> %s"
TABLE_HEADER_C = (
    "<th class=\"dgrid-cell dgrid-cell-padding\">Query Identifier</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\">Clade Classification</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\">Tree Link</th>"
)
TABLE_HEADER_R_RESULT = (
    "<th class=\"dgrid-cell dgrid-cell-padding\">Input FASTA unique ID</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Segment number</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Genotype</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\">Best hit accession</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:13%\">Query coverage %</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Ident %</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">E Value</th>"
)
TABLE_HEADER_R_ERR = (
    "<th class=\"dgrid-cell dgrid-cell-padding\">Input FASTA unique ID</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\" style=\"width:10%\">Segment number</th>"
    "<th class=\"dgrid-cell dgrid-cell-padding\">Error Description</th>"
)
TABLE_ROW = "<td class=\"dgrid-cell dgrid-cell-padding\">%{data}</td>"
TREE_LINK = (
    "<a href=\"%s/view/PhylogeneticTree2/?wsTreeFile=%s/.%s/details/%s&fileType=nwk&isClassification=1&initialValue=%s\" target=\"_blank\">VIEW TREE</a>"
)
TREE_LINK_ALL = (
    "<a href=\"%s/view/PhylogeneticTree2/?wsTreeFile=%s/.%s/details/%s.tre&fileType=nwk&isClassification=1\" target=\"_blank\">VIEW TREE FOR ALL</a>"
)


class SubspeciesClassification:
    def __init__(self, job_file: Path, output_dir: Path):
        """Initialize classifier with job config and prepare paths."""
        self.job_file = job_file
        self.output_dir = output_dir.resolve()
        self.job_data = self.load_job_file()
        self.virus_type = self.job_data["virus_type"]
        self.input_file = self.output_dir / "input.fasta"
        self.output_file_base = self.output_dir / f"{self.job_data['output_file']}"
        self.report_template_path = self.find_lib_path("classification_report.html")
        self.alignment_path = self.find_lib_path("ref-tree-alignment")
        self.rota_genotyper_path = self.find_lib_path("rota-a-genotyper")
        self.id_map = {}  # clean_id -> original_id
        self.env_vars()

    def env_vars(self):
        """Set up environment-dependent variables."""
        self.BASE_URL = os.environ.get("P3_BASE_URL", "https://www.bv-brc.org")

    def find_lib_path(self, name: str) -> Path:
        """Determine library path based on environment."""
        top = os.getenv("KB_TOP")
        paths = [
            Path(top) / "lib" / name,
            Path(top) / "modules" / "bvbrc_subspecies_classification" / "lib" / name,
            Path.home() / "git" / "bvbrc_subspecies_classification" / "lib" / name,
        ]
        for path in paths:
            if path.exists():
                return path
        raise FileNotFoundError(f"Template or directory not found: {name}")

    def load_job_file(self):
        """Read the job configuration JSON file."""
        try:
            with self.job_file.open() as f:
                return json.load(f)
        except Exception as e:
            print(f"Error opening job file: {e}")
            traceback.print_exc()
            sys.exit(-1)

    def setup_output_directory(self):
        """Ensure the output directory exists and set it as working directory."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(self.output_dir)

    def fetch_or_write_input_fasta(self):
        """Retrieve input FASTA from workspace or write it from provided string."""
        if self.job_data["input_source"] == "fasta_file":
            try:
                subprocess.check_call([
                    "p3-cp",
                    f"ws:{self.job_data['input_fasta_file']}",
                    str(self.input_file)
                ])
            except Exception as e:
                print(f"Error copying fasta file from workspace: {e}")
                traceback.print_exc()
                sys.exit(-1)
        elif self.job_data["input_source"] == "fasta_data":
            try:
                self.input_file.write_text(self.job_data["input_fasta_data"])
            except Exception as e:
                print(f"Error writing fasta data into input file: {e}")
                traceback.print_exc()
                sys.exit(-1)
        self.clean_fasta_headers()

    def clean_fasta_headers(self):
        """Sanitize FASTA headers by removing bad characters and replacing spaces."""
        if self.input_file.stat().st_size == 0:
            print("Input fasta file is empty")
            sys.exit(-1)

        cleaned = []
        with self.input_file.open() as f:
            for line in f:
                if line.startswith(">"):
                    original_id = line[1:].strip()
                    clean_id = re.sub(BAD_CHARS_REGEX, "", original_id)
                    clean_id = re.sub(" ", "", clean_id)

                    self.id_map[clean_id] = original_id

                    line = f">{clean_id}\n"
                cleaned.append(line)

        self.input_file.write_text("".join(cleaned))

    def run(self):
        self.setup_output_directory()
        self.fetch_or_write_input_fasta()

        if self.virus_type == "ROTAA":
            self.run_rota_genotyper()
        else:
            self.run_tree_based_classification()

        self.organize_outputs()

    def run_rota_genotyper(self):
        """Run the rotavirus A genotyper and generate HTML report."""

        try:
            subprocess.check_call(["ss-rotaA-genotyper", str(self.input_file)], cwd=self.output_dir)
            report_path = self.output_file_base.with_name(self.output_file_base.name + "_classification_report.html")
            with self.report_template_path.open() as f:
                html_data = f.readlines()

            html_data[60] = REPORT_DATE % datetime.now().strftime("%B %d, %Y %H:%M:%S")

            result_rows = ""
            with (self.output_dir / GENOTYPER_RESULT_F_NAME).open() as f:
                lines = f.readlines()[1:]  # Skip header
                if lines:
                    html_data[68] = TABLE_HEADER_R_RESULT
                    for line in lines:
                        result_rows += "<tr>"
                        for data in line.strip().split("\t"):
                            result_rows += TABLE_ROW.replace("%{data}", data)
                        result_rows += "</tr>"
            html_data[70] = result_rows

            error_rows = ""
            with (self.output_dir / GENOTYPER_ERROR_F_NAME).open() as f:
                lines = f.readlines()[1:]  # Skip header
                if lines:
                    html_data[76] = TABLE_HEADER_R_ERR
                    for line in lines:
                        result_rows += "<tr>"
                        for data in line.strip().split("\t"):
                            error_rows += TABLE_ROW.replace("%{data}", data)
                        error_rows += "</tr>"
            html_data[78] = error_rows

            with report_path.open("w") as f:
                f.writelines(html_data)

        except Exception as e:
            print(f"Error running rotaA genotyper: {e}")
            traceback.print_exc()
            sys.exit(-1)

    def run_tree_based_classification(self):
        # Step 1: Locate reference files in alignment_path/virus_type
        ref_folder = self.alignment_path / self.virus_type
        if not ref_folder.exists():
            print(f"Cannot find reference folder path for virus type {self.virus_type}")
            sys.exit(-1)

        reference_mfa_file = None
        reference_tree_file = None
        stats_file = None
        mapping_file = None

        for file in ref_folder.iterdir():
            if file.suffix in [".aln", ".fasta"]:
                reference_mfa_file = file
            elif file.suffix == ".nh":
                reference_tree_file = file
            elif "info." in file.name:
                stats_file = file
            elif file.suffix == ".tsv":
                mapping_file = file

        # Step 2: Override reference MSA fasta if job_data has 'ref_msa_fasta'
        if "ref_msa_fasta" in self.job_data:
            reference_mfa_file = self.output_dir / "ref_mfa.fasta"
            try:
                subprocess.check_call([
                    "p3-cp",
                    f"ws:{self.job_data['ref_msa_fasta']}",
                    str(reference_mfa_file)
                ])
            except Exception as e:
                print(f"Error copying mfa fasta file from workspace:\n{e}")
                traceback.print_exc()
                sys.exit(-1)

        # Step 3: Run MAFFT to add query sequences to reference MSA
        mafft_output = self.output_dir / MAFFT_OUTPUT_F_NAME
        mafft_cmd = [
            "mafft", "--keeplength", "--add",
            str(self.input_file), str(reference_mfa_file)
        ]
        try:
            with mafft_output.open("w") as o:
                subprocess.check_call(mafft_cmd, stdout=o)
        except Exception as e:
            print(f"Error running mafft:\n{e}")
            traceback.print_exc()
            sys.exit(-1)

        # Step 4: Run pplacer for phylogenetic placement
        pplacer_output = self.output_dir / PPLACER_OUTPUT_F_NAME
        pplacer_cmd = [
            "pplacer", "-m", "GTR",
            "-t", str(reference_tree_file),
            "-s", str(stats_file),
            str(mafft_output),
            "-o", str(pplacer_output)
        ]
        try:
            subprocess.check_call(pplacer_cmd)
        except Exception as e:
            print(f"Error running pplacer:\n{e}")
            traceback.print_exc()
            sys.exit(-1)

        # Step 5: Run guppy to generate trees for query sequences
        guppy_output = self.output_dir / GUPPY_OUTPUT_F_NAME
        guppy_cmd = ["guppy", "sing", str(pplacer_output)]
        try:
            subprocess.check_call(guppy_cmd)
        except Exception as e:
            print(f"Error running guppy:\n{e}")
            traceback.print_exc()
            sys.exit(-1)

        # Step 6: Run cladinator with correct options
        cladinator_output = self.output_dir / CLADINATOR_OUTPUT_F_NAME
        cladinator_cmd = ["cladinator3", str(guppy_output), str(cladinator_output)]

        is_ortho = self.virus_type in ["INFLUENZAH5", "SWINEH1", "SWINEH3", "SWINEH1US"]
        is_adeno = self.virus_type in ["MASTADENOA", "MASTADENOB", "MASTADENOC", "MASTADENOE", "MASTADENOF"]
        is_paramyxo = self.virus_type in ["MEASLES", "MUMPS"]
        is_pox = (self.virus_type == "MPOX")

        if self.virus_type == "INFLUENZAH5":
            cladinator_cmd.insert(1, f"-S={CLADE_DELIMITER_INFLUENZAH5}")
        elif self.virus_type == "MPOX":
            cladinator_cmd.insert(1, f"-S={CLADE_DELIMITER_MPOX}")
        elif self.virus_type == "DENGUE":
            cladinator_cmd.insert(1, f"-S={CLADE_DELIMITER_DENGUE}")
        elif self.virus_type in CLADE_DELIMITER_VIRUS_TYPES:
            cladinator_cmd.insert(1, f"-S={CLADE_DELIMITER}")

        if is_ortho or is_adeno or is_paramyxo or is_pox:
            cladinator_cmd.insert(1, f"-m={str(mapping_file)}")
        if is_adeno or is_paramyxo:
            cladinator_cmd.insert(1, "-x")

        try:
            subprocess.check_call(cladinator_cmd)
        except Exception as e:
            print(f"Error running cladinator:\n{e}")
            traceback.print_exc()
            sys.exit(-1)

        # Step 7: Parse cladinator output to locate and name individual .tre files later
        # "#Tree # \t Query \t Assignment \t Confidence \t Brackets \t Conclusion \t Placement count"
        # We extract "Query" and "Assignment" columns by name
        query_dict = {}
        with cladinator_output.open() as f:
            header = f.readline().strip().split("\t")
            if "Query" not in header or "Assignment" not in header:
                raise ValueError("Unexpected cladinator output format: missing Query/Assignment columns")

            q_idx = header.index("Query")
            a_idx = header.index("Assignment")

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) <= max(q_idx, a_idx):
                    continue  # Skip malformed lines

                query = parts[q_idx].strip()
                assignment = parts[a_idx].strip()
                query_dict[query] = assignment

        # Step 8: Generate per-query .tre files
        id_file_map = {}
        with guppy_output.open() as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                for key in query_dict:
                    if f"{key}_" in line:
                        safe_key = re.sub(r"[^a-zA-Z0-9_.-]", "_", key)
                        if len(safe_key) > 250:
                            print(f"Filename too long for {key}, truncating")
                            safe_key = safe_key[:200] + "_"

                        tre_file = self.output_dir / f"{safe_key}.tre"
                        id_file_map[key] = tre_file.name
                        with tre_file.open("w") as q:
                            q.write(line)
                        break

        # Step 9: Decorate ,tre files to display labels in phylogenetic tree
        if is_ortho or is_adeno or is_paramyxo or is_pox:
            decorator_output = self.output_dir / "outtree.tre"
            decorator_cmd = [
                "decorator", "-f=n", "-nh",
                "INPUT_TRE_FILE", str(mapping_file), str(decorator_output)
            ]
            global_tree = self.output_dir / "out.tree.tre"

            for key, tre_filename in id_file_map.items():
                tre_path = self.output_dir / tre_filename
                try:
                    decorator_cmd[3] = str(tre_path)
                    subprocess.check_call(decorator_cmd)
                    shutil.move(decorator_output, tre_path)

                    # Append to global tree file
                    with tre_path.open() as infile, global_tree.open("a") as outfile:
                        shutil.copyfileobj(infile, outfile)
                        outfile.write("\n")
                except Exception as e:
                    print(f"Error decorating {tre_filename}: {e}")
                    traceback.print_exc()
                    sys.exit(-1)

        # Step 10: Create summary file
        result_path = self.output_dir / "result.tsv"
        with result_path.open("w") as f:
            f.write("Query Identifier\tClade Classification\n")
            for k, v in query_dict.items():
                display_id = self.id_map.get(k, k)
                f.write(f"{display_id}\t{v}\n")

        # Step 11: Create final report file
        report_path = self.output_file_base.with_name(self.output_file_base.name + "_classification_report.html")
        with self.report_template_path.open() as f:
            html_data = f.readlines()

        html_data[60] = REPORT_DATE % datetime.now().strftime("%B %d, %Y %H:%M:%S")
        html_data[68] = TABLE_HEADER_C
        html_data[70] = ""
        for key, val in query_dict.items():
            display_id = self.id_map.get(key, key) # original ID for display
            html_data[70] += "<tr>"
            html_data[70] += TABLE_ROW.replace("%{data}", display_id)
            html_data[70] += TABLE_ROW.replace("%{data}", val)
            initial_val = "" if val.startswith("Sequence") else re.sub("-like$", "", val)
            html_data[70] += TABLE_ROW.replace("%{data}", TREE_LINK % (
                self.BASE_URL, self.job_data["output_path"], self.output_file_base.name, id_file_map[key], initial_val
            ))
            html_data[70] += "</tr>"

        with report_path.open("w") as f:
            f.writelines(html_data)

    def organize_outputs(self):
        """Move data files to 'details' folder."""
        try:
            details_folder = self.output_dir / "details"
            details_folder.mkdir(exist_ok=True)
            for file in self.output_dir.iterdir():
                if file.is_file() and not file.name.endswith("html"):
                    shutil.move(str(file), details_folder / file.name)
        except Exception as e:
            print(f"Error moving files: {e}")
            traceback.print_exc()
            sys.exit(-1)


def parse_args():
    parser = argparse.ArgumentParser(description="Subspecies Classification Script")
    parser.add_argument("-j", "--jfile", required=True, type=Path, help="JSON job file")
    parser.add_argument("-o", "--output", required=False, type=Path, default=Path("."), help="Output directory")
    return parser.parse_args()


def main():
    args = parse_args()
    classifier = SubspeciesClassification(args.jfile, args.output)
    classifier.run()


if __name__ == "__main__":
    main()
