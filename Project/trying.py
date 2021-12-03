import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2)

for ax in fig.axes:
    ax.plot([0, 10], [0, 10], label='linear')

lines, labels = fig.axes[-1].get_legend_handles_labels()

fig.legend(lines, labels, loc='upper center')

print(fig.axes)
plt.show()
