# 🧬 Protein Structure Prediction using Hill Climbing on 2D HP Lattice Model

This project implements a **Hill Climbing algorithm** for predicting simplified **protein folding** structures using the **2D Hydrophobic–Polar (HP) lattice model**.  
It simulates how proteins fold into low-energy conformations based on hydrophobic collapse and polar exposure principles.

---

## 📘 Overview

Protein folding determines the three-dimensional shape and function of a protein.  
In this project, we use a **discrete lattice-based model** to simplify the folding process:

- Each residue is either **Hydrophobic (H)** or **Polar (P)**.
- The sequence forms a **self-avoiding walk (SAW)** on a 2D grid.
- **Hydrophobic–Hydrophobic (H–H)** contacts reduce energy (favorable).
- The **Hill Climbing** algorithm iteratively improves the conformation by exploring neighboring lattice configurations defined by the move set **{L, F, R}** (Left, Forward, Right).

The algorithm searches for a conformation that minimizes the energy:

\[
E(C) = - \sum_{(i,j)} \mathbf{1}_{\text{adj}}(i,j) \cdot \delta(s_i, H) \cdot \delta(s_j, H)
\]

---

## ⚙️ Features

- **2D HP Lattice Model** with self-avoiding walk constraint  
- **Relative Direction Encoding (L, F, R)** for folding paths  
- **Energy Function** based on H–H contacts  
- **Hill Climbing Optimization** with adaptive restarts  
- **Visualization** of protein lattice structures (colored H/P nodes)  
- **Experimental Graphs**: Energy convergence and runtime analysis  


---

## 🧠 Algorithm Outline

1. **Input:** HP sequence (e.g., `"HPHPPHHPHPPHPHHPPHPH"`)
2. **Generate** a random valid fold path (self-avoiding walk)
3. **Iteratively modify** one move (L/F/R) at a random index
4. **Evaluate** the energy of the new conformation  
   - Accept if energy decreases  
   - Reject otherwise  
5. **Terminate** when no improvement occurs after fixed iterations  
6. **Visualize** the final conformation and record energy statistics

---

## 🧩 Example Usage

```bash
# Clone the repository
git clone https://github.com/Kamehamehaaaaa/Protein-2D-Structure-Prediction.git
cd Protein-2D-Structure-Prediction

# Install dependencies
pip install -r requirements.txt

# Run the Hill Climbing algorithm
python hillClimbing.py <HP Sequence>
```

--- 

## Sample Sequences

sequence 1: "HPHPPHHPHPPH"
sequence 2: "HHPPHPHPH"
sequence 3: "HPHPPHHPHPPHPHHPPHPH"
sequence of length 100: "HPPHHPHPHHPPHPHPPHHPPHHPHPHPHHPPPHHPPHPHPHHPPHPPHPHHPPHPPHHPPHPHPPHPPHPHHPPPHHPPHPPHPPHPHHPPHPH"

--- 

📊 Visualization
The program automatically plots the best lattice structures found:
Orange nodes → Hydrophobic (H)
Blue nodes → Polar (P)
Gray lines → Peptide bonds

---

📜 License
© 2025 Rohit Bogulla
This project is released under the MIT License.
Please cite this work if used for research or educational purposes.
