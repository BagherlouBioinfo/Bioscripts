#Step_7
import pandas as pd
import requests
import json

# Load the gene list from Excel file
df = pd.read_excel("Processed/genes_lower_noise_in_PDLSC.xlsx")

# Extract gene names from the 'Gene' column
genes = df["Gene"].dropna().astype(str).tolist()

# Prepare the gene list string (separated by newlines)
gene_list_str = "\n".join(genes)

# Enrichr API endpoint for adding a gene list
add_list_url = 'https://maayanlab.cloud/Enrichr/addList'

# Define the payload for the POST request
payload = {
    'list': gene_list_str,
    'description': 'Genes with lower expression noise in PDLSC'
}

# Send POST request to Enrichr
response = requests.post(add_list_url, files=payload)

# Check if request was successful
if response.status_code == 200:
    data = response.json()
    user_list_id = data['userListId']
    print(f"Gene list successfully submitted to Enrichr.")
    print(f"Enrichr userListId: {user_list_id}")
    print(f"Open in Enrichr: https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}")
else:
    print("❌ Failed to submit gene list to Enrichr.")
    print("Response code:", response.status_code)
    print(response.text)

