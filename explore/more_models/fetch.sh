#!/usr/bin/env bash
curl -L 'https://www.ebi.ac.uk/ena/portal/api/count?format=json&result=read_run&field=instrument_model&_source=ENABrowser' | jq > ena_model.json
