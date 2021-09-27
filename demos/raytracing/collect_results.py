import numpy as np
import subprocess as sub
from lxml import etree
import pandas as pd

# import seaborn as sns
# import matplotlib
# # import matplotlib.patches as mpatches
# matplotlib.use('agg')
# import matplotlib.pyplot as plt
# plt.switch_backend('agg')

# matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['legend.frameon'] = 'True'
# matplotlib.rcParams["legend.framealpha"] = 0.75
# matplotlib.rcParams["legend.fancybox"] = True

# sns.set(color_codes=True, context="poster")
# sns.set_style("ticks")
# sns.set_palette("Set2", 8)

# colors = ["ocean green", "dark pink", "grey", "tan"]
# sns.set_palette(sns.xkcd_palette(colors), desat=.9)

# sns.set_palette("Set2", desat=.9)


SEEDS = range(0, 500)

fit_dict = {}
fit_list = []
for s in SEEDS:  
    try:
        root = etree.parse("output{0}.xml".format(s)).getroot()
        fitness = np.abs(float(root.findall("detail/bot_0/fitness_score")[0].text))
        fit_dict[s] = fitness
        fit_list += [fitness]

    except (IOError, IndexError):
        pass


# print results

# df = pd.DataFrame.from_dict(fit_dict)

sorted_ids = [k for k, v in sorted(fit_dict.items(), key=lambda item: item[1])]
sorted_fits = [v for k, v in sorted(fit_dict.items(), key=lambda item: item[1])]


worst = sorted_ids[0]
best = sorted_ids[-1]

print(best, ": ", fit_dict[best], "; ", worst, fit_dict[worst])

# fig, ax = plt.subplots(1, 1, figsize=(4, 3))

# plt.axhline(SPHERE_FIT, 
#             # ls=":", 
#             lw=1, color=sns.color_palette()[2])
# ax.annotate("default sphere", xy=(4500, SPHERE_FIT-0.2), va="top", ha="center", color=sns.color_palette()[2], fontsize=6.5)  # was: SPHERE_FIT-0.2

# plt.axhline(np.percentile(sphere_fits, 95, interpolation='lower'), ls=":", lw=1, color=sns.color_palette()[2])
# # plt.axhline(np.percentile(sphere_fits, 5), ls="--", lw=1, color=sns.color_palette()[2])

# ax.annotate("evolved shape", xy=(4500, 3.25), va="bottom", ha="center", color=sns.color_palette()[0], fontsize=6.5)

# if X_EVALS_NOT_GEN:
#     a = sns.tsplot(data=df, value="fitness", condition="trial", unit="rank", time="evals", alpha=0.6, lw=0.5, ci=0)
# else:
#     a = sns.tsplot(data=df, value="fitness", condition="rank", unit="trial", time="gen", ci=99)

# # for trial in range(1, 5):
# #     idx = np.array(results["trial"]) == trial
# #     ax.plot(np.array(results["gen"])[idx], np.array(results["fitness"])[idx])


# if X_EVALS_NOT_GEN:
#     ax.set_xlim([-2, 2000])
#     ax.set_ylim([-0.001, 5])
# else:
#     ax.set_xlim([-10, 5010])
#     ax.set_ylim([-0.01, 4.01])
#     ax.set_yticks(range(5))

# ax.set_ylabel('Filial generations')
# ax.set_xlabel('Designs evaluated')

# # handles, labels = ax.get_legend_handles_labels()
# # ax.legend(handles=handles[1:], labels=labels[1:])
# # legend = plt.legend(frameon=True, loc=2)
# # frame = legend.get_frame()
# # frame.set_facecolor('white')

# ax.get_legend().remove()

# sns.despine()
# plt.tight_layout()

# plt.savefig("plots/Evolutionary_Improvement.png".format(EXP_NAME), bbox_inches='tight', dpi=600) # , transparent=True)

# exit()

