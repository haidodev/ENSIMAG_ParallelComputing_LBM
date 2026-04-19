import re
import matplotlib
matplotlib.use("Agg")  # no display required

import matplotlib.pyplot as plt

file_path = "OAR.2811076.stdout"
pattern = re.compile(r"Total time:\s*([\d.]+)\s*seconds")

data = {}
current_process = 1
current_example = 0

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()

        if "Running reference case" in line:
            current_process = 1
            current_example = 0

        elif "Running with" in line:
            current_process = int(re.search(r"(\d+)", line).group(1))

        elif "Running example" in line:
            current_example = int(re.search(r"(\d+)", line).group(1))

        m = pattern.search(line)
        if m:
            t = float(m.group(1))
            key = (current_process, current_example)
            data.setdefault(key, []).append(t)

# ---- average ----
avg = {k: sum(v)/len(v) for k, v in data.items()}

examples = sorted(set(k[1] for k in avg if k[1] != 0))
processes = sorted(set(k[0] for k in avg))

# ---- plot ----
plt.figure()

for ex in examples:
    y = [avg.get((p, ex), None) for p in processes]
    plt.plot(processes, y, marker='o', label=f"ex {ex}")

plt.xlabel("Number of processes")
plt.ylabel("Average time (seconds)")
plt.title("LBM Performance")
plt.legend()
plt.grid()

# ---- save instead of show ----
plt.savefig("lbm_performance.png", dpi=300, bbox_inches="tight")