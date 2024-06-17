import os
import logging
import pandas as pd
import requests
from bs4 import BeautifulSoup

def get_biosample_accession(assembly_id):
    # NCBI Datasets API endpoint for genome assembly summary
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{assembly_id}?format=json"

    # Make the API request
    response = requests.get(url)

    # Check if the request was successful
    if response.status_code != 200:
        return('Missing')

    # Parse the JSON response
    data = response.json()

    biosample_accession = data['assemblies'][0]['assembly']['biosample']['accession']
    return(biosample_accession)


def get_isolation_source(biosample_id):
    url = f'https://www.ncbi.nlm.nih.gov/biosample/{biosample_id}/'
    response = requests.get(url)
    
    if response.status_code != 200:
        return "Failed to retrieve the data"

    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Debug: Print raw HTML to verify the content
    # print(soup.prettify())

    data_elements = soup.find_all('tr')
    data = {}
    for element in data_elements:
        labels = element.find_all('th')
        values = element.find_all('td')

        for label, value in zip(labels, values):
            key = label.get_text(strip=True)
            val = value.get_text(strip=True)
            data[key] = val

    return data.get('isolation source', 'Missing')

def get_geo_loc_name(biosample_id):
    url = f'https://www.ncbi.nlm.nih.gov/biosample/{biosample_id}/'
    response = requests.get(url)
    
    if response.status_code != 200:
        return "Failed to retrieve the data"

    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Debug: Print raw HTML to verify the content
    # print(soup.prettify())

    data_elements = soup.find_all('tr')
    data = {}
    for element in data_elements:
        labels = element.find_all('th')
        values = element.find_all('td')

        for label, value in zip(labels, values):
            key = label.get_text(strip=True)
            val = value.get_text(strip=True)
            data[key] = val

    return data.get('geographic location', 'Missing')

def get_country(s):
    if s == 'Missing':
        return(s)
    else:
        return(s.split(':')[0])

def find_isolation_source(genome_ids, isolation_source_path):
    samples = pd.DataFrame({"genome_id": genome_ids})
    samples['biosample_accession'] = samples['genome_id'].apply(get_biosample_accession)
    samples['isolation_source'] = samples['biosample_accession'].apply(get_isolation_source)
    samples['geo_loc_name'] = samples['biosample_accession'].apply(get_geo_loc_name)
    samples['Country'] = samples['geo_loc_name'].apply(get_country)

    samples.to_csv(isolation_source_path, index=False)
    # samples.to_csv('../data/source_info/' + genus + '/df_ncbi_isolation_src.csv', index = False)