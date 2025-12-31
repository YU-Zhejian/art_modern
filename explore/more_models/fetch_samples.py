import os
from urllib.parse import urlencode

import requests

FIELDS = ["run_accession", "fastq_bytes", "study_accession"]
MODELS = [
    "illumina novaseq 6000",
    "illumina hiseq 4000",
    "nextseq 550",
    "nextseq 2000",
    "illumina hiseq 3000",
    "illumina novaseq x",
    "nextseq 1000",
    "bgiseq-500",
    "illumina novaseq x plus",
    "dnbseq-g400",
    "dnbseq-t7",
    "illumina hiseq 1500",
    "illumina genome analyzer iix",
    "mgiseq-2000rs",
    "illumina iseq 100",
    "dnbseq-g50",
    "bgiseq-50",
    "dnbseq-g400 fast",
    "dnbseq-t10x4rs",
    "onso",
]

def fetch_ena_data(instrument_model="illumina novaseq 6000"):
    respt = "\t".join(FIELDS) + "\n"
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    query_params = {
            "result": "read_run",
            "query": f'instrument_model="{instrument_model}" AND (library_strategy="wgs" OR library_strategy="wxs")',
            "fields": ",".join(FIELDS),
            "limit": "0",
            "format": "tsv",
        }
    data = urlencode(query_params)
    response = requests.post(url, headers=headers, data=data)
    if response.status_code == 200:
        # Remove header line
        lines = response.text.split("\n")
        print(f"MODEL={instrument_model}: Success with {len(lines) - 2} records")
        respt += "\n".join(lines[1:])
    else:
        print(f"Failed: {response.status_code}: {response.text}")
    return respt


if __name__ == "__main__":
    os.makedirs("fetch_samples.d", exist_ok=True)
    for model_name in MODELS:
        print(f"Fetching data for model: {model_name}")
        result = fetch_ena_data(model_name)
        if result:
            with open(
                os.path.join("fetch_samples.d", f"{model_name.replace(' ', '_')}.tsv"),
                "w",
                encoding="utf-8",
            ) as f:
                f.write(result)
