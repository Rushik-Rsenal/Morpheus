import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import pyplot as plt
import sys
import numpy as np

plt.rcParams.update(
    {
        "font.size": 11,
        "font.family": 'serif',
        "font.serif": ['Liberation Serif'],
        "mathtext.fontset": 'stix',
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.figsize": (6,4),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.linewidth": 1.2,
        "grid.linewidth": 0.8,
        "lines.linewidth": 2,
    }
)

df = pd.read_excel(f"{sys.argv[2]}")

df = df.set_index("Celltype")

stressor = "DMSO"
feature = "mito_tmrm intensity within mtg"

layout = """
    000000000000
    111122223333
    444455556666
"""

fig, axes = plt.subplot_mosaic(layout, figsize=(30, 20), gridspec_kw={"wspace": 1, "hspace": 0.5})

# # 2. Add text title to the top left (Area '0')
axes['0'].axis('off') # Hide grid lines
axes['0'].text(0.5, 0.5, f"{feature} ({stressor})", 
               fontsize=50, ha='center', va='center')


def well_level_data_analysis_log_fold_change(df, col_name, timepoint, stressor, letter):
  global axes

  # Extract rows from a particular timepoint
  df_init = df
  df_extracted = df_init[df_init["Timepoint"] == timepoint]
  df_extracted = df_extracted[df_extracted["Stressor"] == stressor]
  
  # Extract the required columns
  df_extracted = df_extracted[["Celltype", "Cellline ID","Batch", col_name]]
  # Calculate the average of Control for each Batch
  batches = df_extracted["Batch"].unique()
  df_extracted = df_extracted.set_index("Batch")

  # Control_AVG
  control_avg_dict = {}
  cellline_df_list = []

  for batch in batches:
    selected_batch = df_extracted.loc[batch]

    # Get the control rows
    control_rows = selected_batch[selected_batch["Celltype"] == "Control"]

    # Get the average
    avg_df = control_rows.groupby(["Celltype"])[col_name].mean()

    # Store this in the dictionary created
    control_avg_dict[batch] = avg_df["Control"]

    # Divide each row in the selected_batch by the mean of controls
    selected_batch_normalised = selected_batch
    selected_batch_normalised[col_name] = np.log(selected_batch[col_name]/control_avg_dict[batch])


    # # Now get the average for cellines
    cellline_df = selected_batch_normalised.groupby(["Cellline ID", "Celltype"])[col_name].mean()

    # # Append to the list
    cellline_df_list.append(cellline_df)

  # Combine the cellline_df_list
  combined_cellline_df = pd.concat(cellline_df_list, axis=0)
  combined_cellline_df = combined_cellline_df.reset_index()
  combined_cellline_df = combined_cellline_df.sort_values(by="Celltype")

  # Visualise the stuff 
  combined_cellline_df_excluded = combined_cellline_df[combined_cellline_df["Celltype"] != "Control"]
  sns.barplot(data=combined_cellline_df_excluded, y="Cellline ID", x=col_name, hue="Celltype", ax=axes[str(timepoint)]).set_title(f"{letter}) log-fold change of {col_name} {timepoint}")
  plt.figure(figsize=(8, 10),)
  plt.show()



# Plot your graphs
def return_img():
  well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 1, "DMSO", "a")
  well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 2, "DMSO", "b")
  well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 3, "DMSO", "c")
  well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 4, "DMSO", "d")
  well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 5, "DMSO", "e")
  well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 6, "DMSO", "f")
  return fig.savefig(f"{feature}_{stressor}.png")

return_img()

