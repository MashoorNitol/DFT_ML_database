import os
import random
import numpy as np
import shutil



def create_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        
def perturb_coordinates(coordinates,coordchange):
    perturbed_coords = []
    for coord in coordinates:
        perturbed_coord = [num + random.uniform(-coordchange, coordchange) for num in coord]
        perturbed_coords.append(perturbed_coord)
    return perturbed_coords

def perturb_lattice_vectors(vector,lmax,amax):

    v1 = vector[0]
    v2 = vector[1]
    v3 = vector[2]
    p_a = np.linalg.norm(v1)
    p_b = np.linalg.norm(v2)
    p_c = np.linalg.norm(v3)
    alpha = np.arccos(np.dot(v2, v3) / (p_b * p_c))
    beta = np.arccos(np.dot(v1, v3) / (p_a * p_c))
    gamma = np.arccos(np.dot(v1, v2) / (p_a * p_b))

    lx = p_a
    p_xy = p_b * np.cos(gamma)
    p_xz = p_c * np.cos(beta)

    ly = np.sqrt(p_b ** 2 - p_xy ** 2)
    p_yz = (p_b * p_c * (np.cos(alpha)) - p_xy * p_xz) / (ly)

    lz = np.sqrt(p_c ** 2 - p_xz ** 2 - p_yz ** 2)

    cellMat_randomstrain = np.zeros((3, 3), dtype=np.double)
    cellMat_randomstrain[0, 0] = lx + lmax * (2 * np.random.random(1)[0] - 1.0)
    cellMat_randomstrain[0, 1] = 0.0
    cellMat_randomstrain[0, 2] = 0.0

    cellMat_randomstrain[1, 0] = p_xy + amax * (2 * np.random.random(1)[0] - 1.0)
    cellMat_randomstrain[1, 1] = ly + lmax * (2 * np.random.random(1)[0] - 1.0)
    cellMat_randomstrain[1, 2] = 0.0

    cellMat_randomstrain[2, 0] = p_xz + amax * (2 * np.random.random(1)[0] - 1.0)
    cellMat_randomstrain[2, 1] = p_yz + amax * (2 * np.random.random(1)[0] - 1.0)
    cellMat_randomstrain[2, 2] = lz + lmax * (2 * np.random.random(1)[0] - 1.0)

    newvol = np.linalg.det(cellMat_randomstrain)
    return cellMat_randomstrain, newvol

def read_simulation_box(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    lattice_vectors = [list(map(float, lines[i].split()[:3])) for i in range(2, 5)]
    return lattice_vectors

def read_poscar(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    index_direct = lines.index("Direct\n")
    num_atoms = int(lines[index_direct - 1].strip())
    element_type = lines[index_direct - 2].strip()
    coordinates = [list(map(float, line.split()[:3])) for line in lines[index_direct + 1: index_direct + 1 + num_atoms]]

    return element_type, num_atoms, coordinates

def write_poscar(file_path, element_type, num_atoms, coordinates, lattice_vectors):
    with open(file_path, 'w') as f:
        f.write("Perturbed " + "\n" + "1.0" + "\n")
        for vector in lattice_vectors:
            f.write(f"{vector[0]:.9f} {vector[1]:.9f} {vector[2]:.9f}\n")
        f.write(element_type + "\n")
        f.write(str(num_atoms) + "\n")
        f.write("Direct\n")
        for coord in coordinates:
            f.write(f"{coord[0]:.9f} {coord[1]:.9f} {coord[2]:.9f}\n")

if __name__ == "__main__":
    input_file_path = "initial_poscar"  # Replace with the path to your input POSCAR file
    coordchange_values = np.linspace(0.001, 0.05, 10)
    lmax = 0.3
    amax = 0.1

    element_type, num_atoms, coordinates = read_poscar(input_file_path)
    lattice_vectors = read_simulation_box(input_file_path)

    all_input_folder = "all_input"
    create_folder(all_input_folder)

    oldvol = np.linalg.det(lattice_vectors)

    with open("volumechanges.txt", "w") as volume_file:
        volume_file.write("CoordChange VolumeChange(%)\n")

        for i, coordchange in enumerate(coordchange_values):
            input_folder = f"input_{i:04d}"
            create_folder(os.path.join(all_input_folder, input_folder))

            perturbed_coordinates = perturb_coordinates(coordinates, coordchange)
            perturbed_lattice_vectors, newvol = perturb_lattice_vectors(lattice_vectors, lmax, amax)

            output_file_path = os.path.join(all_input_folder, input_folder, "Perturbed_POSCAR")
            write_poscar(output_file_path, element_type, num_atoms, perturbed_coordinates, perturbed_lattice_vectors)

            changevol = 100 * ((abs(newvol - oldvol)) / newvol)
            volume_file.write(f"{coordchange:.6f} {changevol:.6f}\n")


if __name__ == "__main__":
    input_file_path = "initial_poscar"  # Replace with the path to your input POSCAR file
    coordchange_values = np.linspace(0.0001, 0.05, 1000)
    lmax = 0.3
    amax = 0.1

    element_type, num_atoms, coordinates = read_poscar(input_file_path)
    lattice_vectors = read_simulation_box(input_file_path)

    all_input_folder = "all_input"
    create_folder(all_input_folder)

    # Create a folder for storing the copied POSCAR files
    poscars_folder = os.path.join(all_input_folder, "all_POSCARS")
    create_folder(poscars_folder)

    oldvol = np.linalg.det(lattice_vectors)

    with open("volumechanges.txt", "w") as volume_file:
        volume_file.write("CoordChange VolumeChange(%)\n")

        for i, coordchange in enumerate(coordchange_values):
            input_folder = f"input_{i:04d}"
            create_folder(os.path.join(all_input_folder, input_folder))

            perturbed_coordinates = perturb_coordinates(coordinates, coordchange)
            perturbed_lattice_vectors, newvol = perturb_lattice_vectors(lattice_vectors, lmax, amax)

            output_file_path = os.path.join(all_input_folder, input_folder, "POSCAR")
            write_poscar(output_file_path, element_type, num_atoms, perturbed_coordinates, perturbed_lattice_vectors)

            changevol = 100 * ((abs(newvol - oldvol)) / newvol)
            volume_file.write(f"{coordchange:.6f} {changevol:.6f}\n")

            # Copy the Perturbed_POSCAR file to the "all_POSCARS" folder with a unique identifier
            unique_file_name = f"POSCAR_{i:04d}"
            shutil.copy(output_file_path, os.path.join(poscars_folder, unique_file_name))
            
            # Copy INCAR and KPOINTS files to the "input_{i:04d}" folder
            shutil.copy("INCAR", os.path.join(all_input_folder, input_folder, "INCAR"))
            shutil.copy("KPOINTS", os.path.join(all_input_folder, input_folder, "KPOINTS"))
