# ADAR-sensor: AI-powered RNA Sensor Design Platform

An ADAR editing sensor design and optimization platform combining ViennaRNA physics engine with AI reasoning.

## Main Components

### 1. ViennaRNA Visualization Skill (`.gemini/skills/viennarna-visualization`)
- **Core Function**: Visualize dsRNA hybrid structures where sensors bind within large mRNA
- **Operating Principles**:
  1. Optimal binding site search with `duplexfold` (MFE-based)
  2. Dot-bracket notation parsing
  3. Custom matplotlib rendering

### 2. C3 Sensor Project (`Projects/C3_Sensor/`)
- ALU-based sensor design and validation
- dsRNA dimer visualization
- Candidate scoring and ranking

## Directory Structure

```
ADAR-sensor/
├── SKILL.md                 # ViennaRNA visualization skill description
├── manifest.json            # Skill metadata
├── scripts/
│   ├── visualize_binding.py # Core visualization engine
│   ├── viennarna_enhanced.py
│   └── plot_native_dimer.py
└── c3_sensor/
    ├── plot_dsRNA.py
    ├── plot_native_dimer.py
    ├── plot_native_dimer_styled.py
    ├── visualize_c3_elite.py
    ├── c3_300bp_candidates_scored.csv
    ├── top6_c3_sites.txt
    └── human_C3_mRNA.fasta
```

## Usage

### Sensor Binding Visualization
```bash
python scripts/visualize_binding.py --mRNA <target_mRNA.fasta> --sensor <sensor_sequence>
```

### C3 Sensor Analysis
```bash
python c3_sensor/visualize_c3_elite.py
```

## Scientific Background

- **ADAR (Adenosine Deaminase Acting on RNA)**: Post-transcriptional RNA editing enzyme
- **C3 ALU Mimetic Sensor**: Sensor designed based on ALU sequences
- **ViennaRNA**: Physicochemical model for RNA secondary structure prediction and analysis

## License

For research use only.