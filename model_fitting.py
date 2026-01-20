import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score


#  PERCOLATION MODEL (t FIXED = 4.5)

def percolation_model_fixed_t(x, A, B, t=4.5):
    """
    Modified percolation model with fixed exponent t.
    A : effective mobility prefactor
    B : aggregation / crosslinking penalty
    """
    x_c = 0.05
    return A * (x - x_c)**t * np.exp(-B * x)



# DATA

data_registry = {

    "LiCl (Li+)": {
        "valency": "Monovalent",
        "charge": 1,
        "x": np.array([0.10, 0.20, 0.30, 0.40, 0.50]),
        "sigma": np.array([4.98e-7, 3.55e-6, 3.61e-5, 6.61e-4, 2.08e-3]),
        "fit": True
    },

    "NaN3 (Na+)": {
        "valency": "Monovalent",
        "charge": 1,
        "x": np.array([0.05, 0.10, 0.15, 0.20]),
        "sigma": np.array([1.85e-7, 2.70e-6, 2.40e-6, 2.20e-6]),
        "fit": False  # intentionally excluded, check paper for info...
    },

    "NH4I (NH4+)": {
        "valency": "Monovalent",
        "charge": 1,
        "x": np.array([0.30, 0.40, 0.50, 0.60, 0.70]),
        "sigma": np.array([5.1e-6, 9.7e-5, 3.7e-4, 1.2e-3, 4.5e-3]),
        "fit": True
    },

    "Mg(NO3)2 (Mg2+)": {
        "valency": "Divalent",
        "charge": 2,
        "x": np.array([0.30, 0.40, 0.50, 0.60]),
        "sigma": np.array([4.86e-6, 6.86e-5, 7.70e-4, 7.53e-4]),
        "fit": True
    },

    "ZnCl2 (Zn2+)": {
        "valency": "Divalent",
        "charge": 2,
        "x": np.array([0.50, 0.60, 0.70]),
        "sigma": np.array([6.72e-4, 4.49e-3, 4.21e-3]),
        "fit": True
    }
}



# FITTING + VPI...

fit_table = []
raw_data_table = []

print("\n" + "="*110)
print(f"{'System':<18} | {'A':<8} | {'B':<6} | {'VPI=B/A':<10} | {'RMSE':<10} | {'R²':<6}")
print("-"*110)

for system, d in data_registry.items():

    
    for xi, si in zip(d["x"], d["sigma"]):
        raw_data_table.append([system, d["charge"], xi, si])

    if not d["fit"]:
        print(f"{system:<18} | {'--':<8} | {'--':<6} | {'--':<10} | {'--':<10} | {'--':<6}")
        continue

    x = d["x"]
    y = d["sigma"]

    popt, _ = curve_fit(
        lambda x, A, B: percolation_model_fixed_t(x, A, B),
        x, y,
        p0=[0.05, 0.5],
        bounds=([1e-4, 0.01], [1.0, 10.0]),
        maxfev=30000
    )

    A, B = popt
    y_fit = percolation_model_fixed_t(x, A, B)

    rmse = np.sqrt(mean_squared_error(y, y_fit))
    r2 = r2_score(y, y_fit)
    vpi = B / A

    fit_table.append([system, d["charge"], A, B, vpi, rmse, r2])

    print(f"{system:<18} | {A:<8.3f} | {B:<6.2f} | {vpi:<10.2f} | {rmse:<10.1e} | {r2:<6.3f}")

print("="*110)



# FIGURE

plt.figure(figsize=(10, 7))
x_smooth = np.linspace(0.05, 0.80, 400)

colors = {
    "LiCl (Li+)": "red",
    "NaN3 (Na+)": "orange",
    "NH4I (NH4+)": "gold",
    "Mg(NO3)2 (Mg2+)": "purple",
    "ZnCl2 (Zn2+)": "brown"
}

markers = {"Monovalent": "o", "Divalent": "^"}

for system, d in data_registry.items():

    plt.scatter(
        d["x"], d["sigma"],
        s=120,
        edgecolor="black",
        marker=markers[d["valency"]],
        color=colors[system],
        label=system
    )

    if d["fit"]:
        row = next(r for r in fit_table if r[0] == system)
        A, B = row[2], row[3]
        plt.plot(
            x_smooth,
            percolation_model_fixed_t(x_smooth, A, B),
            color=colors[system],
            linewidth=2.5
        )

plt.yscale("log")
plt.xlabel("Salt Weight Fraction, x", fontsize=14)
plt.ylabel(r"Ionic Conductivity $\sigma$ (S cm$^{-1}$)", fontsize=14)
plt.title("Salt-Resolved Ion Transport in Pectin Biopolymer Electrolytes", fontsize=16, fontweight="bold")
plt.grid(True, which="both", alpha=0.3)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("Fig1_Salt_Resolved_Conductivity_FINAL.svg", dpi=300)
plt.show()



#SI EXPORT 



with open("SI_Raw_Digitized_Data.csv", "w") as f:
    f.write("System,Charge,Salt_Weight_Fraction,Conductivity_S_cm\n")
    for system, d in data_registry.items():
        for xi, si in zip(d["x"], d["sigma"]):
            f.write(f"{system},{d['charge']},{xi:.2f},{si:.3e}\n")

# FITTED PARAMETERS 
with open("SI_Fitted_Parameters_and_VPI.csv", "w") as f:
    f.write("System,Charge,A,B,VPI,RMSE_S_cm,R2\n")
    for row in fit_table:
        system, charge, A, B, VPI, rmse, r2 = row
        f.write(f"{system},{charge},{A:.4f},{B:.4f},{VPI:.2f},{rmse:.3e},{r2:.3f}\n")

print("\nSI FILES GENERATED (NO PANDAS REQUIRED):")
print("• SI_Raw_Digitized_Data.csv")
print("• SI_Fitted_Parameters_and_VPI.csv")
