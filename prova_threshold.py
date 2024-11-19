from scipy.stats import friedmanchisquare
import pandas as pd

# Input data
residue_data = {
    'PPI flex': [0.746, 0.740, 0.734, 0.728, 0.721, 0.715, 0.706, 0.697, 0.684],
    'PPI rec': [0.725, 0.718, 0.711, 0.704, 0.696, 0.688, 0.679, 0.669, 0.661],
    'PPI flex (Ab)': [0.64, 0.641, 0.635, 0.629, 0.640, 0.627, 0.617, 0.606, 0.595],
    'PPI rec (Ab)': [0.672, 0.667, 0.660, 0.651, 0.624, 0.618, 0.613, 0.612, 0.618],
    'Ab flex fine': [0.789, 0.789, 0.788, 0.788, 0.788, 0.787, 0.785, 0.783, 0.774],
    'Ab rec fine': [0.787, 0.784, 0.781, 0.780, 0.776, 0.774, 0.766, 0.759, 0.744],
    'Ab flex': [0.849, 0.849, 0.850, 0.851, 0.851, 0.850, 0.845, 0.841, 0.831],
    'Ab rec': [0.844, 0.843, 0.842, 0.838, 0.835, 0.832, 0.824, 0.814, 0.798]
}

atomic_data = {
    'PPI flex': [0.761, 0.758, 0.755, 0.752, 0.750, 0.747, 0.744, 0.739, 0.735],
    'PPI rec': [0.735, 0.732, 0.725, 0.719, 0.714, 0.706, 0.696, 0.685, 0.673],
    'PPI flex (Ab)': [0.718, 0.715, 0.716, 0.716, 0.715, 0.713, 0.712, 0.714, 0.717],
    'PPI rec (Ab)': [0.733, 0.734, 0.734, 0.735, 0.734, 0.734, 0.732, 0.727, 0.719],
    'Ab flex fine': [0.851, 0.850, 0.848, 0.847, 0.846, 0.841, 0.838, 0.834, 0.827],
    'Ab rec fine': [0.788, 0.786, 0.786, 0.786, 0.784, 0.782, 0.776, 0.769, 0.762],
    'Ab flex': [0.828, 0.827, 0.826, 0.823, 0.818, 0.811, 0.798, 0.777, 0.742],
    'Ab rec': [0.854, 0.854, 0.853, 0.850, 0.847, 0.842, 0.831, 0.815, 0.795]
}

thresholds = [60, 62.5, 65, 67.5, 70, 72.5, 75, 77.5, 80]

# Perform Friedman Test for each model
def perform_friedman_test(data, measure_type):
    print(f"\nFriedman Test Results for {measure_type} Data:")
    for model, values in data.items():
        stat, p = friedmanchisquare(*[[v] for v in values])  # Use *[[v] for v in values] to handle multiple thresholds
        print(f"Model: {model}, Statistic = {stat:.4f}, p-value = {p:.4f}")

# Residue
perform_friedman_test(residue_data, "Residue")

# Atomic
perform_friedman_test(atomic_data, "Atomic")
