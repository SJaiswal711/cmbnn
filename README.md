# Beam Convolution Simulations

This repository contains Fortran and Python codes for simulating satellite scanning of HEALPix sky maps with different beam models:

- **Elliptical Gaussian Beam** (idealized, parametric)
- **Realistic Beam** (based on pre-computed response grid)

The outputs can be used for forward modeling of time-ordered data (TOD), response matrices, and validation of scanning strategies.

---

## 📂 Repository Structure

```

.
├── README.md
├── RealBeam
│   ├── main.f90                     # Main MPI program for realistic beam simulation           
│   ├── constant.f90                 # Parameter definitions (angles, frequencies, steps)           
│   ├── math_library.f90
│   ├── subroutines.f90              # Subroutines for vectors, HEALPix queries, and time-step processing 
│   ├── neighbours_mat.f90           # Extracts neighboring pixel for response matrix neighbouring pixels
├── EllipticalBeam
│   ├── main.f90                     # Main program for elliptical beam simulation           
│   ├── global_variables.f90         # Parameter definitions (angles, frequencies, steps)           
│   ├── math_library.f90
├── HFI_ScanBeam_143-1a_R2.00.fits   # Example Planck HFI beam (used to generate grid.txt for real-beam code)
├── map.dat                          # Example HEALPix input map (used in elliptical convolution demo)
├── yearly_scan_frequency.f90        # Program to compute scan hit-count map (how many times each pixel is scanned)
├── python/
│   ├── gen_cmb_maps.py            # Generate 1000 CMB maps from CAMB Cl
│   ├── convolve_maps.py           # Convolve foreground and CMB maps
│   └── ...

````

---

## ⚙️ Code Descriptions

### 🔹 Elliptical Beam Convolution (`EllipticalBeam/main.f90`)
- **Input:**
  - `map.dat` → HEALPix sky map to be scanned.
- **Process:**
  - Simulates scanning with an **elliptical Gaussian beam** (`fwhm_x`, `fwhm_y`).
  - At each time step, computes:
    - Pointing direction (`pix_ring`).
    - Convolved temperature (detector signal).
- **Output:**
  - `convolved_map.dat` containing:
    ```
    time_step   pixel   convolved_temperature
    ```

---

### 🔹 Real Beam Convolution (`RealBeam/main.f90`)
- **Input:**
  - `grid.txt` → Pre-computed real beam response grid.
  - (Optionally) `map.dat` for testing.
- **Process:**
  - Uses MPI to parallelize over sky pixels.
  - Computes a **response matrix**:
    - Each row corresponds to a pixel.
    - Each entry gives beam weights for neighboring HEALPix pixels.
  - Stores intermediate results per rank (`beam_response_mat_0.dat … beam_response_mat_47.dat`).
- **Output:**
  - `beam_response_mat_{rank}.dat` files (one per MPI rank), containing:
    ```
    node_id   pixel   count   weight1   weight2   ...
    ```

---

### 🔹 Response Matrix Neighbors (`RealBeam/neighbors_mat.f90`)
- **Process:**
  - Extracts the **neighboring HEALPix pixel indices** for each response matrix pixel.
- **Output:**
  - `neighbors_mat_{rank}.dat` files with neighbor pixel indice.

---

### 🔹 Yearly Scan Frequency
- Generated as part of both simulations.
- Stores **hit counts**: how many times each HEALPix pixel was visited during the one-year scan.
- Useful for coverage maps and validation.

---

### 🔹 Python Scripts
- **`gen_cmb_maps.py`**: Generate CMB realizations from CAMB power spectra (`Cl`).
- **`convolve_maps.py`**: Convolve CMB + foreground maps with the beam model.
- Additional scripts for pre/post-processing.

---

## 🚀 Typical Workflow

1. **Generate or provide input maps:**
   - Place sky maps in `map.fits`.
   - Place beam grid in `grid.txt` (for real beam mode).

2. **Run elliptical beam scan:**
   **`gfortran EllipticalBeam/main.f90 EllipticalBeam/global_variables.f90 EllipticalBeam/math_library.f90 -o elliptical_convolution
   ./elliptical_convolution`**

→ Produces `convolved_map.dat`.

3. **Run real beam scan with MPI:**

   **`
   mpif90 RealBeam/math_library.f90 RealBeam/constants.f90 RealBeam/subroutines.f90 RealBeam/main.f90 -O3 -o real_beam
   mpirun -np 48 ./real_beam
   `**

   → Produces `beam_response_mat_0.dat … beam_response_mat_{rank}.dat`.

4. **Extract neighbors for response matrix (optional):**

   **`
   gfortran RealBeam/neighbors_matrix.f90 -o extract_neighbors
   ./extract_neighbors
   `**

   → Produces `neighbors_mat_0.dat … neighbors_mat_47.dat`.

---

## ⚠️ File Size Warning

* The `.dat` outputs (`beam_response_mat_*.dat`, `neighbors_mat_*.dat`, `convolved_map.dat`) can be **hundreds of GB** depending on scan duration and resolution.
* To reduce size, use a lower NSIDE (e.g., NSIDE=256) or adjust the scan length in constants.f90.
---

## 📌 Summary

* **Elliptical Beam Code** → produces **time-ordered convolved map**.
* **Real Beam Code** → produces **response matrix** (weights for later convolution).
* **Neighbors Code** → links response matrix pixels back to HEALPix neighbors.
* **Python Tools** → generate input maps and analyze outputs.
