# Data-Driven Model Predictive Control (DDMPC) Framework

This repository provides a structured implementation and comparison of multiple **data-driven Model Predictive Control (MPC)** strategies for nonlinear systems.

The goal is to demonstrate how different data-driven paradigms can be used to design predictive controllers **without relying on explicit parametric models**, and to benchmark their performance under a unified simulation framework.

---

## 🔍 Overview

We consider a nonlinear discrete-time system of the form:

\[
x_{k+1} = f(x_k, u_k), \quad y_k = g(x_k, u_k)
\]

Instead of identifying an explicit model, we design controllers directly from data using different paradigms:

- Behavioral / regression-based
- Kalman-filter-based estimation
- Optimization-based formulations
- Neural-network-based approximation

This repository provides **four implementations**:

| Method     | Description |
|------------|-------------|
| **DDRMPC** | Direct Data Regression MPC |
| **KF-DDMPC** | Kalman Filter-based Data-Driven MPC |
| **DDMPCA** | Optimization-based Data-Driven MPC (QP formulation) |
| **NN-DDMPC** | Neural Network-based Data-Driven MPC |

---

## 📁 Repository Structure

```text
data_driven_mpc/
│
├── README.md
├── LICENSE
│
├── src/
│   ├── DDRMPC/
│   │   └── ddrmpc.m
│   │
│   ├── KFDDMPC/
│   │   └── kfddmpc.m
│   │
│   ├── DDMPCA/
│   │   └── ddmpca.m
│   │
│   └── NNDDMPC/
│       ├── nnddmpc.m
│       └── mpccon.m
│
├── tests/
│   └── test1.m   % Main simulation script
│
└── docs/
    ├── theoretical_foundations.pdf
    ├── simulation_study_ddmpc.pdf
    └── simulation_results_ddmpc_comparison.pdf
  

---


## ▶️ Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/data_driven_mpc.git
cd data_driven_mpc

How to run:
addpath(genpath(pwd));
tests/test1



---

### 🧠 Methods + Study + Notes + Author**
```markdown
---

## 🧠 Methodological Details

### 1. DDRMPC (Direct Data Regression MPC)

- Constructs predictive models directly from input-output data
- Uses regression-based lifting
- No explicit state estimation

### 2. KF-DDMPC

- Integrates a Kalman filter for state estimation
- Combines data-driven prediction with recursive filtering
- Improved robustness to noise

### 3. DDMPCA (Optimization-Based)

- Formulates the problem as a Quadratic Program (QP)
- Uses explicit cost minimization with constraints
- Solved using `quadprog`

### 4. NN-DDMPC

- Uses neural-network-inspired structure for prediction
- Solves nonlinear optimization via `fmincon`
- Captures nonlinear dynamics more flexibly

---

## 📊 Simulation Study

The simulation framework evaluates controllers based on:

- Tracking performance
- Constraint satisfaction
- Input smoothness
- Robustness to nonlinearities

See:

- `docs/simulation_study_ddmpc.pdf`
- `docs/simulation_results_ddmpc_comparison.pdf`

---

## 📌 Notes

- The current implementation is primarily designed for **SISO systems**
- Extensions to MIMO systems require careful dimension handling
- The framework is intended for **research and benchmarking purposes**

---

## 🚀 Research Perspective

This repository reflects a broader research direction:

- Data-driven control without explicit identification
- Integration of learning and optimization
- Bridging MPC, estimation, and machine learning

It can serve as a foundation for:

- Real-time control implementations
- Embedded MPC (PLC / real-time systems)
- Energy systems and smart grid applications

---

## 📜 License

This project is licensed under the MIT License.  
See the `LICENSE` file for details.

---

## 👤 Author

**Roozbeh Abolpour**  
Control & Optimization Researcher  
GitHub: https://github.com/Roozbeh-Abolpour

---

## ⭐ Contribution

Contributions, issues, and suggestions are welcome.

If you find this repository useful, consider giving it a star.