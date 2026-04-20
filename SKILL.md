# ViennaRNA Visualization: Native Dimer & Local Duplex Mapping

This skill combines ViennaRNA physics engine with custom Python rendering logic to visualize precise dsRNA hybrid structures where sensors bind within large mRNA sequences.

## Core Operating Principles (3-Stage Scientific Reasoning)

This skill doesn't just draw pictures - it identifies the physical reality of binding sites through the following 3-stage logical process.

### 1. Optimal Binding Site Search (`duplexfold`)
Find the optimal coordinates for sensor binding within large mRNA (e.g., 5,200+ nt). Use `RNA.duplexfold` algorithm to compare sensor sequence against the entire mRNA, physically calculating and extracting the **thermodynamically lowest energy and most stable (Minimum Free Energy) binding site**.

### 2. Dot-bracket Notation Parsing
Precisely parse the binding structure data computed by the physics engine.
- **Parentheses `(` and `)`**: Indicate **'Match/Binding'** state where two bases form hydrogen bonds.
- **Dot `.`**: Indicates **'Mismatch'** or **'Bulge'** state where bases lack pairing partners or energy doesn't allow binding.

### 3. Data-driven Custom Visualization and Rendering
Translate parsed data into visual language using Python code.
- **Classification and Coloring**: Match regions (parentheses) are connected with gray lines and treated with light gray beads, while Bulge regions (dots) are emphasized with **intense red beads** for immediate structural defect detection.
- **Coordinate Layout (`get_xy_coordinates`)**: Based on binding information (Match/Bulge), calculate 2D coordinates (X, Y) by computing physical forces that pull binding regions together and push bulge regions apart.
- **Readability Optimization**: Expand strand spacing by 2x and finely adjust bead (1/10) and font (1/5) sizes to expose structural details of dsRNA.

## Core Engine and Implementation (Core Engine)

All scientific reasoning and rendering logic in this skill is executed by the following core Python engine.

- **Execution Engine**: `scripts/visualize_binding.py`
- **Role**: 
    - Pass user-provided mRNA and sensor sequence to ViennaRNA physics engine.
    - Parse `duplexfold` and `cofold` results to generate match/bulge data.
    - Perform custom rendering with matplotlib including bead size adjustment and strand spacing expansion.
- **Modification Guide**: Modify this file to change visualization ratios (bead size, etc.) or physics calculation parameters.

## Standard Workflow
1. **Input**: Receive target mRNA (full) and sensor RNA sequence.
2. **Fragmenting**: Automatically extract approximately 300bp target fragment including 50bp flanks on both sides.
3. **Layout**: Generate physical coordinates using linear (Naview) layout to prevent strand ends from tangling.
4. **Output**: Generate high-resolution PNG map with bulge emphasis along with binding energy (MFE) in `outcome/` folder.

## Engineering Guidelines
- **Zero-Vision Guessing**: All red dots are 'actual non-binding sites' computed by the physics engine - pure data results without user bias.
- **Scalability**: This approach is designed to always find optimal binding sites and perform precise mapping regardless of mRNA length.
