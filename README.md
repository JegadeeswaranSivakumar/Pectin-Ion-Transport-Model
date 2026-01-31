
# Generalized Transport Modeling of Monovalent and Divalent Ion Conduction

This repository contains the source code and digitized datasets for the **Generalized Transport Framework** developed to analyze ion conduction in pectin-based polymer electrolytes.

The model utilizes a modified continuum percolation theory to quantify transport suppression and introduces the **Valency Penalty Index (VPI)** as a comparative metric for monovalent and divalent systems.

## 📄 Associated Publication

[Waiting...](https://doi.org/10.21203/rs.3.rs-8651772/v1)

## 🚀 Key Features

* **Modified Percolation Model:** Implements the transport equation: .
* **Valency Penalty Index (VPI):** Calculates the dimensionless parameter  to quantify charge-driven transport suppression.
* **Automated Fitting:** Uses `scipy.optimize` (Bounded Non-Linear Least Squares) to extract transport parameters from raw conductivity data.
* **Visualization:** Generates salt-resolved conductivity plots and comparisons.

## 📂 Repository Structure

* `data/`: Contains digitized conductivity datasets (`.csv`) for Li, NH4, Mg, and Zn systems.
* `src/`: Main Python scripts for model fitting and plotting.
* `results/`: Output plots and parameter logs.

## 💻 Usage

### Prerequisites

* Python 3.8+
* NumPy, SciPy, Matplotlib, Pandas

### Installation

```bash
git clone https://github.com/YourUsername/Pectin-Ion-Transport-Model.git
cd Pectin-Ion-Transport-Model
pip install -r requirements.txt

```

### Running the Model

```bash
python src/model_fitting.py

```

*This script will load the datasets, perform the non-linear regression, and output the fitted parameters (, , VPI) to the console and `results/` folder.*

## 🔗 Citation

Waiting...

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.

---

