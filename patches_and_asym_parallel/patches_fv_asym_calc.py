import os
import json
import ray
import tempfile
import argparse
import logging
from pathlib import Path
from propermab import defaults
from propermab.features import feature_utils
from propermab.features import StructFeaturizer
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import subprocess
import tempfile
import time

# Setup logging
logging.basicConfig(filename="pdb_processing.log", level=logging.ERROR, 
                    format="%(asctime)s - %(levelname)s - %(message)s")

# Initialize Ray with all available CPUs
ray.init(num_cpus=os.cpu_count())

# # Load system config once (reduces overhead)
# defaults.system_config.update_from_json('./default_config.json')

def get_pdb_files(input_path):
    """Retrieve PDB files from a directory or return the file if it's a single PDB."""
    input_path = Path(input_path)
    if input_path.is_dir():
        return list(input_path.glob("*.pdb"))
    elif input_path.is_file() and input_path.suffix == ".pdb":
        return [input_path]
    else:
        raise ValueError("Invalid input path. Must be a PDB file or a directory containing PDB files.")
    

def add_hydrogens(pdb_file_path, pH=7.0):
    """Fixes a PDB file by adding missing hydrogens at a specified pH and returns the path to the fixed file."""
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_pdb:
        fix_name = temp_pdb.name  # Temporary file path

    fixer = PDBFixer(filename=pdb_file_path)
    fixer.addMissingHydrogens(pH=pH)
    
    with open(fix_name, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    return fix_name

def calculate_features(fixed_pdb_path):
    """Computes molecular features from a fixed PDB file."""
    mol_feature = feature_utils.calculate_patch_features(fixed_pdb_path)
    featurizer = StructFeaturizer(fixed_pdb_path)
    mol_feature['charge_asym'] = featurizer.fv_chml()

    return mol_feature

def run_immunopdb(fixed_pdb_path, use_temp=False, log_output=False):
    """Runs ImmunoPDB.py on a fixed PDB file and captures output/errors.

    Args:
        fixed_pdb_path (str): Path to the fixed PDB file.
        use_temp (bool): Whether to use a temporary file for output. Default is False.

    Returns:
        str or None: Path to the ImmunoPDB output file if successful, otherwise None.
    """
    if use_temp:
        temp_file = tempfile.NamedTemporaryFile(suffix="_immuno_output.pdb", delete=False)
        immuno_output_file = temp_file.name
        temp_file.close()  # Close immediately, as subprocess will write to it
    else:
        immuno_output_file = fixed_pdb_path.replace(".pdb", "_immuno_output.pdb")

    try:
        # Run ImmunoPDB
        result = subprocess.run(
            ["python", "ImmunoPDB.py", "-i", fixed_pdb_path, "-o", immuno_output_file],
            capture_output=True,
            text=True
        )

        # Log output and errors
        with open("immunopdb.log", "a") as log_file:
            log_file.write(f"Processing {fixed_pdb_path}:\n")
            if log_output:
                log_file.write(result.stdout)
            if result.stderr:
                log_file.write("\nErrors:\n" + result.stderr)
                logging.error(f"ImmunoPDB error for {fixed_pdb_path}: {result.stderr}")

        return immuno_output_file if result.returncode == 0 else None

    except Exception as e:
        logging.error(f"Failed to run ImmunoPDB on {fixed_pdb_path}: {e}")
        return None


@ray.remote
def process_pdb_file(pdb_file_path, repeats=1, use_temp=False, pH=7.0):
    """Processes a PDB file multiple times and computes molecular features."""
    pdb_file_path = str(pdb_file_path)
    all_results = []
    defaults.system_config.update_from_json('../default_config.json')

    # Initialize variables to None
    fixed_pdb_path = None
    immuno_output = None

    for i in range(repeats):
        try:
            fixed_pdb_path = add_hydrogens(pdb_file_path, pH=pH)  # Now always assigned before use
            immuno_output = run_immunopdb(fixed_pdb_path, use_temp=use_temp)

            if immuno_output:
                mol_feature = calculate_features(immuno_output)
                all_results.append(mol_feature)
            else:
                logging.error(f"Skipping feature calculation for {pdb_file_path} due to ImmunoPDB failure.")

        except Exception as e:
            logging.error(f"Error processing {pdb_file_path} (Iteration {i+1}): {e}")

        finally:
            # Ensure only defined variables are used
            for file in [fixed_pdb_path, immuno_output]:
                if file and os.path.exists(file):
                    os.remove(file)

    return {pdb_file_path: all_results}


def run_processing(input_path, repeats=1, pH=7.0, use_ray_wait=False, output_file="molecular_features.json"):
    """Runs the PDB processing pipeline with optional `ray.wait()`."""
    start_time = time.time()  # Capture the start time

    pdb_files = get_pdb_files(input_path)
    tasks = [process_pdb_file.remote(str(pdb), repeats, pH) for pdb in pdb_files]

    try:
        if use_ray_wait:
            results = []
            while tasks:
                done, tasks = ray.wait(tasks, num_returns=1)
                results.extend(ray.get(done))
        else:
            results = ray.get(tasks)

        # Save results to JSON
        with open(output_file, "w") as f:
            json.dump(results, f, indent=4)

        elapsed_time = time.time() - start_time
        print(f"Processing complete in {elapsed_time} seconds. Results saved to {output_file}")

    except Exception as fatal_error:
        print(f"ERROR: {fatal_error}")
        logging.critical(f"ERROR: {fatal_error}")
        exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files using Ray for parallel execution.")
    
    parser.add_argument("input", type=str, help="Path to a PDB file or a directory containing PDB files.")
    parser.add_argument("-r", "--repeats", type=int, default=1, help="Number of times to repeat feature calculations for each structure.")
    parser.add_argument("-w", "--wait", action="store_true", help="Use ray.wait() to process tasks incrementally.")
    parser.add_argument("-o", "--output", type=str, default="molecular_features.json", help="Output file for the processed results.")
    parser.add_argument("--pH", type=float, default=7.0, help="pH value for adding missing hydrogens.")

    args = parser.parse_args()

    run_processing(args.input, repeats=args.repeats, use_ray_wait=args.wait, output_file=args.output, pH=args.pH)
