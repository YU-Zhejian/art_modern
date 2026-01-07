import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ABSENT = -1
MAX_RLEN = 300


if __name__ == "__main__":
    df_existing = pd.read_parquet("filtered_samples_readlengths.parquet")

    # Convert using following rule:
    # None -> [ABSENT, ABSENT]
    # [X] -> [X, ABSENT]
    # [X, Y] -> [X, Y]
    def convert_readlens(readlens):
        if readlens is None:
            return [ABSENT, ABSENT]
        elif len(readlens) == 1:
            return [readlens[0], ABSENT]
        elif len(readlens) >= 2:
            return [readlens[0], readlens[1]]
        else:
            return [ABSENT, ABSENT]

    df_existing["read_lengths"] = df_existing["read_lengths"].apply(convert_readlens)
    # Split into two columns
    df_existing[["read_length_1", "read_length_2"]] = pd.DataFrame(
        df_existing["read_lengths"].tolist(),
        index=df_existing.index,
    )
    x_values = [0, 50, 100, 150, 250, 300]
    colors = sns.color_palette("husl", len(x_values))
    # Plot using matplotlib with instrument_model as x-axis
    # Set fig size

    models = list(filter(lambda x: x is not None, set(df_existing["instrument_model"])))
    fig, axes = plt.subplots(len(models), 2, sharey=False, figsize=(16, 2 * len(models)), sharex=True)

    for i, model in enumerate(models):
        ax = axes[i, 0]
        df = df_existing.query(f"`instrument_model` == '{model}'")
        sns.histplot(data=df, x="read_length_1", ax=ax, bins=30, binrange=(ABSENT, MAX_RLEN))
        ax.set_title(f"RLen1 Dist. for {model} (n={sum(df['read_length_1'] > ABSENT)}/{df.shape[0]})")
        ax = axes[i, 1]
        sns.histplot(data=df, x="read_length_2", ax=ax, bins=30, binrange=(ABSENT, MAX_RLEN))
        ax.set_title(f"RLen2 Dist. for {model} (n={sum(df['read_length_2'] > ABSENT)}/{df.shape[0]})")
    for ax in axes.ravel():
        for x, c in zip(x_values, colors):
            ax.axvline(x=x, color=c, linestyle="--", linewidth=1)
    plt.savefig("read_length_1_by_instrument_model.pdf")
    plt.clf()
