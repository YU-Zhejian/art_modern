import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Load the data from the JSON file
    df = pd.read_json("ena_model.json")

    df.query(
        'instrument_model.str.contains("nextseq") |'
        'instrument_model.str.contains("illumina") |'
        'instrument_model.str.contains("hiseq") |'
        'instrument_model.str.contains("dnbseq") |'
        'instrument_model.str.contains("bgiseq") |'
        'instrument_model.str.contains("onso") |'
        'instrument_model.str.contains("mgiseq")',
        inplace=True,
        engine="python",
    )

    # Convert count to numeric (it's stored as string)
    df["count"] = pd.to_numeric(df["count"])

    # Sort by count in descending order
    df_sorted = df.sort_values("count", ascending=False)
    df_sorted.to_markdown("instrument_model_counts.md", index=False)
    df_sorted.to_csv("instrument_model_counts.csv", index=False)

    # Create the plot
    plt.figure(figsize=(12, 8))
    plt.barh(df_sorted["instrument_model"], df_sorted["count"])
    plt.xlabel("Count")
    plt.ylabel("Instrument Model")
    plt.title("Instrument Model Counts (High to Low)")
    plt.gca().invert_yaxis()  # Invert y-axis to have highest count at the top
    plt.tight_layout()

    # Show the plot
    plt.savefig("instrument_model_counts.pdf")
