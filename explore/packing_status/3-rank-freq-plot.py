import pandas as pd
import matplotlib.pyplot as plt

pkgs_to_highlight = [
    "conda-forge::numpy",
    "conda-forge::pandas",
    "bioconda::samtools",
    "bioconda::art",
    "bioconda::art_modern",
    "bioconda::art_modern-openmpi",
    "bioconda::bwa",
    "bioconda::gatk4",
]

if __name__ == "__main__":
    for name in ["median", "total", "conda_velocity"]:
        df = pd.read_parquet(f"packing_status_{name}.parquet")
        for channel in ["conda-forge", "bioconda"]:
            pkgs_to_highlight_channel = {
                pkg.replace(f"{channel}::", "") for pkg in pkgs_to_highlight if pkg.startswith(f"{channel}::")
            }
            this_df = df.query(f"`data_source` == '{channel}'").query("`counts` >= 10")
            this_df.sort_values("counts", ascending=False, inplace=True)
            top_and_low = this_df.iloc[[0, -1]]
            print(f"\n{name.capitalize()} packing status for channel {channel}:")
            print(top_and_low[["pkg_name", "counts"]])
            pkgs_to_highlight_channel.update(top_and_low["pkg_name"].tolist())

            this_df.loc[:, "RANK"] = range(1, len(this_df) + 1)
            plt.figure(figsize=(10, 6))
            plt.loglog(this_df["RANK"], this_df["counts"], linestyle="-", alpha=0.5)
            this_df_highlight = this_df.query(f"`pkg_name` in @pkgs_to_highlight_channel")
            plt.loglog(
                this_df_highlight["RANK"],
                this_df_highlight["counts"],
                marker="o",
                color="red",
                linestyle="none",
            )

            for _, row in this_df_highlight.iterrows():
                pkg_name = row["pkg_name"]
                rank = row["RANK"]
                count = row["counts"]
                plt.text(rank, count, pkg_name, fontsize=9, ha="right", va="bottom")

            plt.xlabel("Rank")
            plt.ylabel(f"{name.capitalize()} Packing Status Count")
            plt.title(f"Rank-Frequency Plot of {name.capitalize()} Packing Status Counts")
            plt.grid(True, which="both", ls="--", lw=0.5)
            plt.savefig(f"packing_status_{name}_{channel}_rank_freq_plot.pdf")
            plt.cla()
