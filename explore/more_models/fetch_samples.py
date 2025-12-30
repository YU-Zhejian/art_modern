import os

import pandas as pd
import requests
from urllib.parse import urlencode

E_COLI_TAX_ID = 562

def fetch_ena_data(instrument_model="illumina novaseq 6000", tax_id=E_COLI_TAX_ID):
    url = "https://www.ebi.ac.uk/ena/portal/api/search"

    headers = {
        "Content-Type": "application/x-www-form-urlencoded"
    }

    query_params = {
        'result': 'read_run',
        'query': f'tax_eq({tax_id}) AND instrument_model="{instrument_model}" AND library_strategy="wgs"',
        'fields': 'run_accession,experiment_title,tax_id,fastq_bytes,fastq_ftp,fastq_md5,study_accession',
        'limit': '1000',
        'format': 'tsv'
    }

    data = urlencode(query_params)

    response = requests.post(url, headers=headers, data=data)

    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed: {response.status_code}")
        return None

if __name__ == "__main__":
    os.makedirs("fetch_samples.d")
    model_counts_df = pd.read_csv("instrument_model_counts.csv")
    for model_name in model_counts_df["instrument_model"]:
        print( f"Fetching data for model: {model_name}" )
        result = fetch_ena_data(model_name)
        if result:
            with open(os.path.join("fetch_samples.d", f"{model_name.replace(' ', '_')}.tsv"), "w", encoding="utf-8") as f:
                f.write(result)
