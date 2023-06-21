# DFT Training Database Generation for RANN

## Steps in the `vasp_data_generate` Folder:

1. Replicate the unit cell (POSCAR) using the `replicate.py` code. This will create a file called `replicated_poscar`.

2. Create your INCAR and KPOINTS files.

3. Run the following command: `python3 thermal_data.py`

4. The `thermal_data.py` script will read the `replicated_poscar` file and generate `n` files based on user inputs. The generated files will have varying temperature, strain from thermal effects, length and angle increments.

5. Each generated folder will contain a copy of the INCAR and KPOINTS files. You can also add a POTCAR file to the desired directory.

