import numpy as np
from matplotlib import pyplot as plt

model = np.array(list(map(float, open("log10_data_model.txt", "r").read().split())))
print("# data:", len(model), "\tmin value:", min(model), "\tmax value", max(model))
print("RMS:", np.sqrt(np.mean(model**2)))
print("Q1:", np.quantile(model, 0.25), "Q2:", np.quantile(model, 0.5), "Q3:", np.quantile(model, 0.75))
print("# val between -0.5 and 0.5:", np.count_nonzero(np.abs(model) <= 0.5), "\tthe percentage:", np.count_nonzero(np.abs(model) <= 0.5) / len(model))

#generate histogram
f, (ax_hist, ax_box) = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (.85, .15)})
ax_box.boxplot(model, vert=False)
ax_hist.hist(model, bins = np.linspace(-1.6, 1, 52))
plt.show()