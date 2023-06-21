import os
import shutil
import numpy as np

def vec2angle(vec1, vec2):
    angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle = np.arccos(angle)
    return angle

def generate_POSCARs(filename, lmax, amax, tempstart, tempincr, filestart, fileend, incr, num_iterations):
    rad2deg = 180 / np.pi

    natoms = np.genfromtxt(filename, dtype='i8', skip_header=6, max_rows=1)
    scaling = np.genfromtxt(filename, dtype='f8', skip_header=1, max_rows=1)
    basevect = np.genfromtxt(filename, dtype='f8', skip_header=2, max_rows=3)
    ntype = np.genfromtxt(filename, dtype='<U11', skip_header=5, max_rows=1)
    x, y, z = np.genfromtxt(filename, skip_header=8, unpack=True, max_rows=natoms, usecols=[0, 1, 2])
    axis_matrix = np.array(basevect)
    oldvol = np.linalg.det(axis_matrix)
    v1 = basevect[:, 0]
    v2 = basevect[:, 1]
    v3 = basevect[:, 2]
    p_a = np.linalg.norm(v1)
    p_b = np.linalg.norm(v2)
    p_c = np.linalg.norm(v3)
    alpha = vec2angle(v2, v3)
    beta = vec2angle(v1, v3)
    gamma = vec2angle(v1, v2)

    lx = p_a
    p_xy = p_b * np.cos(gamma)
    p_xz = p_c * np.cos(beta)

    ly = np.sqrt(p_b ** 2 - p_xy ** 2)
    p_yz = (p_b * p_c * (np.cos(alpha)) - p_xy * p_xz) / (ly)

    lz = np.sqrt(p_c ** 2 - p_xz ** 2 - p_yz ** 2)
    cellMat_fixed = np.zeros((3, 3), dtype=np.double)
    cellMat_fixed[0, 0] = lx
    cellMat_fixed[0, 1] = 0.0
    cellMat_fixed[0, 2] = 0.0

    cellMat_fixed[1, 0] = p_xy
    cellMat_fixed[1, 1] = ly
    cellMat_fixed[1, 2] = 0.0

    cellMat_fixed[2, 0] = p_xz
    cellMat_fixed[2, 1] = p_yz
    cellMat_fixed[2, 2] = lz

    # Create directory for all-poscars
    if os.path.exists("all-poscars"):
        shutil.rmtree("all-poscars")
    os.mkdir("all-poscars")
    for k in range(num_iterations):
        j = tempstart + tempincr * k

        for i in range(filestart, fileend + 1):
            dirname = f"input_{i}"
            if os.path.exists(dirname):
                shutil.rmtree(dirname)
            os.mkdir(dirname)
            shutil.copy("inputs/INCAR", dirname)
            shutil.copy("inputs/KPOINTS", dirname)
            # Generate new lattice parameters
            delta_a = np.random.uniform(-lmax, lmax)
            delta_b = np.random.uniform(-lmax, lmax)
            delta_c = np.random.uniform(-lmax, lmax)

            delta_alpha = np.random.uniform(-amax, amax)
            delta_beta = np.random.uniform(-amax, amax)
            delta_gamma = np.random.uniform(-amax, amax)

            lx_new = lx + delta_a
            ly_new = ly + delta_b
            lz_new = lz + delta_c

            alpha_new = alpha + delta_alpha
            beta_new = beta + delta_beta
            gamma_new = gamma + delta_gamma

            v1_new = [lx_new, 0.0, 0.0]
            v2_new = [ly_new * np.cos(gamma_new), ly_new * np.sin(gamma_new), 0.0]
            v3_new = [
                lz_new * np.cos(beta_new),
                lz_new * (np.cos(alpha_new) - np.cos(beta_new) * np.cos(gamma_new)) / np.sin(gamma_new),
                lz_new * np.sqrt(
                    1.0
                    - np.cos(beta_new) ** 2
                    - ((np.cos(alpha_new) - np.cos(beta_new) * np.cos(gamma_new)) / np.sin(gamma_new)) ** 2
                ),
            ]
            # Construct new lattice vectors matrix
            axis_matrix_new = np.array([v1_new, v2_new, v3_new])
            volume_new = np.linalg.det(axis_matrix_new)

            scaling_factor = np.cbrt(oldvol / volume_new)

            # Scale the lattice vectors
            axis_matrix_new_scaled = axis_matrix_new * scaling_factor

            # Update the atomic positions
            x_new = x * scaling_factor + j*(2*np.random.random(1)[0]-1.0)
            y_new = y * scaling_factor + j*(2*np.random.random(1)[0]-1.0)
            z_new = z * scaling_factor + j*(2*np.random.random(1)[0]-1.0)


            # Write the modified POSCAR file
            poscar_file = open(f"{dirname}/POSCAR", "w")
            poscar_file.write("Modified POSCAR\n")
            poscar_file.write("1.0\n")
            np.savetxt(poscar_file, axis_matrix_new_scaled, fmt="%.16f", delimiter="   ", newline="\n")
            poscar_file.write(f"{ntype}\n")
            poscar_file.write(f"{natoms}\n")
            poscar_file.write("Direct\n")
            np.savetxt(poscar_file, np.column_stack((x_new, y_new, z_new)), fmt="%.16f", delimiter="   ", newline="\n")
            poscar_file.write("\n")
            poscar_file.close()

            # Copy the POSCAR file to all-poscars directory
            shutil.copy(f"{dirname}/POSCAR", f"all-poscars/POSCAR_{i}")

        filestart += incr
        fileend += incr

    print("All POSCAR files generated successfully.")

# usage
filename = 'replicated_POSCAR'
lmax = 0.2  # EACH length change in between +-0.1
amax = 0.05  # EACH angle change in between 2 degrees from ideal
tempstart = 0.008
tempincr = 0.005
filestart = 1
fileend = 100
incr = 100
num_iterations = 5 # 5*100 = 500 files

generate_POSCARs(filename, lmax, amax, tempstart, tempincr, filestart, fileend, incr, num_iterations)    
