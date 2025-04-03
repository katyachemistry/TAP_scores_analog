# Parallelized TAP-scores Analog Based on Propermab

`propermab` is a Python package for calculating and predicting molecular features and properties of monoclonal antibodies (mAbs).

This fork introduces a **parallelized** version for computing hydrophobic, negative, and positive patches, as well as Vh-Vl charge asymmetry from antibody PDB structures—analogous to **Therapeutic Antibody Profiler (TAP) scores**.

This implementation has been tested on **ABodyBuilder3**-generated structures and requires PDB files to have chain names **H** (heavy) and **L** (light). **IMGT format numbering is not required**, as it is handled in-place using the **ImmunoPDB** module from `anarci`.

---
## Installation (Linux)

Follow the installation instructions from the original `propermab` repository:

### 1. Set up a Conda environment
```bash
conda env create -f propermab/conda_env.yml
conda activate propermab
```

### 2. Install `propermab`
```bash
pip install -e propermab/
```

---
## APBS Installation
The **Adaptive Poisson-Boltzmann Solver (APBS) v3.0.0** is used by `propermab` to calculate electrostatic potentials. Download and extract it:
```bash
wget https://github.com/Electrostatics/apbs/releases/download/v3.0.0/APBS-3.0.0_Linux.zip -O apbs.zip
unzip apbs.zip
```
Record the installation path for use in the configuration step.

---
## Configuration
**Rename** the default_config_copy.json to **default_config.json** and edit it to specify the correct paths:
```json
{
    "hmmer_binary_path": "",
    "nanoshaper_binary_path": "/APBS_PATH/APBS-3.0.0.Linux/bin/NanoShaper",
    "apbs_binary_path": "/APBS_PATH/APBS-3.0.0.Linux/bin/apbs",
    "pdb2pqr_path": "pdb2pqr",
    "multivalue_binary_path": "/APBS_PATH/APBS-3.0.0.Linux/share/apbs/tools/bin/multivalue",
    "atom_radii_file": "",
    "apbs_ld_library_paths": ["LIB_PATH", "/APBS_PATH/APBS-3.0.0.Linux/lib/"]
}
```

### Finding Required Paths
- **`hmmer_binary_path`**: Run the following command to find the correct path:
  ```bash
  dirname $(which hmmscan)
  ```
- **`atom_radii_file`**: This should point to `amber.siz`, found in the `pdb2xyzr.zip` archive, available at: [Electrostatics Zone](https://electrostaticszone.eu/downloads/scripts-and-utilities.html). It is also included in this repo. 

### Setting Up `LIB_PATH`
To determine the correct `LIB_PATH`, install `readline 7.0` in a separate Conda environment:
```bash
conda deactivate
conda env create --name readline python=3.8
conda activate readline
conda install readline=7.0
```
Find the path by running:
```bash
echo ${CONDA_PREFIX}/lib/
```
Replace `APBS_PATH` in `default_config.json` with your actual APBS installation directory.

### Final Steps
Deactivate the `readline` environment, reactivate the `propermab` environment, and install `ray` for parallel execution:
```bash
conda deactivate
conda activate propermab
pip install ray==2.10.0
```

---
## Example Usage
Navigate to the `patches_and_asym_parallel` directory and run an example script:
```bash
python patches_fv_asym_calc.py example/ --repeats 3 --wait --output molecular_features_example.json --pH 7.5
```

### Arguments:
- **`--repeats`**: Number of repeats for each structure. Outputs may vary due to randomness in hydrogen addition and minimization.
- **`--wait`**: Uses `ray.wait()` to process tasks incrementally.
- **`--pH`**: Specifies the pH for hydrogen addition.

The PDB file should have:
- **Heavy chain** named as **H**
- **Light chain** named as **L**

---
## Third-party Software Dependencies
`propermab` relies on third-party software that must be installed separately. These dependencies may have their own license requirements, which users should review before installation and usage.

---
### Notes
- Ensure all paths are correctly configured in `default_config.json` before running calculations.
- `ray` enables efficient parallel execution for large-scale computations.

---