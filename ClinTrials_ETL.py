#!/usr/bin/env python3

import pandas as pd
import requests
import bs4
from bs4 import BeautifulSoup
import re
import collections
import os
import json
import numpy as np
import pickle
from functools import reduce
import time
# import concurrent
import concurrent.futures
import multiprocessing
import datetime as dt
from datetime import date
import pathlib
import configparser
import sys
import urllib
import zipfile
from thefuzz import fuzz # fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python

global sublist_length 
sublist_length = 985 # Name Resolver takes batches of 1000

global CAS_SERVERURL 
global II_SKR_SERVERURL 
global METAMAP_INTERACTIVE_URL 
global stserverurl 
global tgtserverurl
global apikey 
global serviceurl 
global ksource 

CAS_SERVERURL = "https://utslogin.nlm.nih.gov/cas/v1"
II_SKR_SERVERURL = 'https://ii.nlm.nih.gov/cgi-bin/II/UTS_Required'
METAMAP_INTERACTIVE_URL = II_SKR_SERVERURL + "/API_MM_interactive.pl"
stserverurl = "https://utslogin.nlm.nih.gov/cas/v1/tickets"
tgtserverurl = "https://utslogin.nlm.nih.gov/cas/v1/api-key"
serviceurl = METAMAP_INTERACTIVE_URL
ksource = '2020AB'
cfg = configparser.ConfigParser()
cfg.read('config.ini')
apikey = cfg['METAMAP']['apikey']
apikey = apikey.strip("\''")


def get_nr_response(chunk):
    
    """Runs Name Resolver"""
    nr_url = 'https://name-resolution-sri.renci.org/lookup'
    nr_terms = []
    for term in chunk:
        nr_term = {}
        params = {'string':term, 'limit':1} # limit -1 makes this return all available equivalent CURIEs name resolver can give            
        r = requests.post(nr_url, params=params)
        try:
            res = r.json()
            if res:
                for key, val in res.items():
                    nr_term[term] = [key, val[0]]
                    nr_terms.append(nr_term)
            else:
#                 print(term + " unable to be mapped by Name Resolver")
                pass
        except Exception as e:
#             print(e)
#             print(term + " unable to be mapped by Name Resolver")
            pass      
    time.sleep(5)
    return nr_terms

def run_parallel_threads_nr(unmapped_chunked):
    # multithread implementation for retrieving Name Resolver responses
    # Create a ThreadPoolExecutor with the desired number of threads
    with concurrent.futures.ThreadPoolExecutor(multiprocessing.cpu_count() - 1) as executor:
        # Submit the get_response() function for each item in the list
        futures = [executor.submit(get_nr_response, chunk) for chunk in unmapped_chunked]
        # Retrieve the results as they become available
        output = [future.result() for future in concurrent.futures.as_completed(futures)]
    return output

def split_list(lst, sublist_length):
    return [lst[i:i+sublist_length] for i in range(0, len(lst), sublist_length)]


def get_token_sort_ratio(str1, str2):
    """ fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python """
    try:
        return fuzz.token_sort_ratio(str1, str2)
    except:
        return None    

def get_token_set_ratio(str1, str2):
    """ fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python """
    try:
        return fuzz.token_set_ratio(str1, str2)
    except:
        return None  

def get_similarity_score(str1, str2):
    """ fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python """
    try:
        return fuzz.ratio(str1, str2)
    except:
        return None

def get_raw_ct_data():
    term_program_flag = True
    global data_dir
    global data_extracted
    
    # get all the links and associated dates of upload into a dict called date_link
    url_all = "https://aact.ctti-clinicaltrials.org/pipe_files"
    response = requests.get(url_all)
    soup = BeautifulSoup(response.text)
    body = soup.find_all('option') #Find all
    date_link = {}
    for el in body:
        tags = el.find('a')
        try:
            zip_name = tags.contents[0].split()[0]
            date = zip_name.split("_")[0]
            date = dt.datetime.strptime(date, '%Y%m%d').date()
            date_link[date] = tags.get('href')
        except:
            pass
    latest_file_date = max(date_link.keys())   # get the date of the latest upload
    url = date_link[latest_file_date]   # get the corresponding download link of the latest upload so we can download the raw data
    date_string = latest_file_date.strftime("%m_%d_%Y")
    data_dir = "{}/data".format(pathlib.Path.cwd())
    data_extracted = data_dir + "/{}_extracted".format(date_string)
    data_path = "{}/{}_pipe-delimited-export.zip".format(data_dir, date_string)
    
    if not os.path.exists(data_path):   # if folder containing most recent data doesn't exist, download and extract it into data folder
        
        term_program_flag = False   # flag below for terminating program if latest download exists (KG is assumed up to date)
        print("Downloading Clinical Trial data as of {}".format(date_string))
        response = requests.get(url)
        if response.status_code == 200:
            with open(data_path, 'wb') as file:
                file.write(response.content)
            print("Finished download of zip")
            with zipfile.ZipFile(data_path, 'r') as download:
                print("Unzipping data")
                download.extractall(data_extracted)
        else:
            print("KG is already up to date.")
    return {"term_program_flag": term_program_flag, "data_extracted_path": data_extracted}

def read_raw_ct_data(flag_and_path):
    if flag_and_path["term_program_flag"]:
        print("Exiting program. Assuming KG has already been constructed from most recent data dump from AACT.")
#         exit()
#         pass
    else:
        data_extracted = flag_and_path["data_extracted_path"]
        # read in pipe-delimited files 
        conditions_df = pd.read_csv(data_extracted + '/conditions.txt', sep='|', index_col=False, header=0)
        interventions_df = pd.read_csv(data_extracted + '/interventions.txt', sep='|', index_col=False, header=0)
        browse_conditions_df = pd.read_csv(data_extracted + '/browse_conditions.txt', sep='|', index_col=False, header=0)
        browse_interventions_df = pd.read_csv(data_extracted + '/browse_interventions.txt', sep='|', index_col=False, header=0)
        
#     ### GET RID OF....CHEAT LINE FOR TESTING
        conditions_df = conditions_df.iloc[:5000]
        interventions_df = interventions_df.iloc[:5000]

    return {"conditions": conditions_df, "interventions": interventions_df, "browse_conditions": browse_conditions_df, "browse_interventions": browse_interventions_df}

def exact_match_mesh(df_dict):
    
    # -------    CONDITIONS    ------- #
    conditions = df_dict["conditions"]
    browse_conditions = df_dict["browse_conditions"] 

    tomap_conditions = conditions["downcase_name"].values.tolist()
    tomap_conditions = list(set(tomap_conditions))
    print("Number of unique conditions in this Clinical Trials data dump: {}".format(len(tomap_conditions)))
    mesh_exact_mapped = list(set(tomap_conditions).intersection(browse_conditions.downcase_mesh_term.unique()))
    print("Number of unique conditions that have an exact MeSH term match given in this dump: {}".format(len(mesh_exact_mapped)))

    print("since these MeSH terms don't come with identifers, we retrieve them from Name Resolver")
    mesh_exact_mapped_chunked = split_list(mesh_exact_mapped, sublist_length)
    mesh_exact_mapped_curied = run_parallel_threads_nr(mesh_exact_mapped_chunked)

    mesh_exact_mapped_curied = [element for sublist in mesh_exact_mapped_curied for element in sublist] # flatten the list of lists
    mesh_exact_mapped_curied = {key: value for dictionary in mesh_exact_mapped_curied for key, value in dictionary.items()}

    mapped_conditions = pd.DataFrame({"condition_input": list(mesh_exact_mapped_curied.keys()), # get dataframe of exact MeSH mapped conditions
                                      "condition_CURIE_id": [value[0] for value in mesh_exact_mapped_curied.values()],
                                      "condition_CURIE_name": [value[-1] for value in mesh_exact_mapped_curied.values()],
                                      "source": "MeSH term exact mapped, Name Resolver CURIE"})
    unmapped_conditions = list(set(tomap_conditions)-set(mapped_conditions.condition_input))
    print("Number of unique conditions that are unmapped after finding exact MeSH mappings and using Name Resolver to get CURIES: {}".format(len(unmapped_conditions)))

    # -------    INTERVENTIONS    ------- #
    interventions = df_dict["interventions"]
    browse_interventions = df_dict["browse_interventions"] 

    tomap_interventions = interventions["name"].values.tolist()
    tomap_interventions = reduce(lambda a, b: a+[str(b)], tomap_interventions, [])
    tomap_interventions = [string.lower() for string in tomap_interventions] # lowercase the strings
    tomap_interventions = list(set(tomap_interventions))
    print("Number of unique interventions in this Clinical Trials data dump: {}".format(len(tomap_interventions)))
    mesh_exact_mapped = list(set(tomap_interventions).intersection(browse_interventions.downcase_mesh_term.unique()))
    print("Number of unique interventions that have an exact MeSH term match given in this dump: {}".format(len(mesh_exact_mapped)))

    print("since these MeSH terms don't come with identifers, we retrieve them from Name Resolver")
    mesh_exact_mapped_chunked = split_list(mesh_exact_mapped, sublist_length)
    mesh_exact_mapped_curied = run_parallel_threads_nr(mesh_exact_mapped_chunked)
    mesh_exact_mapped_curied = [element for sublist in mesh_exact_mapped_curied for element in sublist] # flatten the list of lists
    mesh_exact_mapped_curied = {key: value for dictionary in mesh_exact_mapped_curied for key, value in dictionary.items()}
    mapped_interventions = pd.DataFrame({"intervention_input": list(mesh_exact_mapped_curied.keys()),    # get dataframe of exact MeSH mapped interventions
                                         "intervention_CURIE_id": [value[0] for value in mesh_exact_mapped_curied.values()],
                                         "intervention_CURIE_name": [value[-1] for value in mesh_exact_mapped_curied.values()],
                                         "source": "MeSH term exact mapped, Name Resolver CURIE"})

    unmapped_interventions = list(set(tomap_interventions)-set(mapped_interventions.intervention_input))
    print("Number of unique interventions that are unmapped after finding exact MeSH mappings and using Name Resolver to get CURIES: {}".format(len(unmapped_interventions)))

    ct_terms = {'mapped_conditions': mapped_conditions,
                'unmapped_conditions': unmapped_conditions,
                'mapped_interventions': mapped_interventions,
                'unmapped_interventions': unmapped_interventions}
    return ct_terms


def inexact_match_mesh(df_dict, ct_terms):
    
    # get dataframes bc I'm going to compute fuzzy scores and dump into columns
    # find unmapped terms AND THEIR CORRESPONDING NCITS!
    # get the conditions that have exact MESH term matches, and conditions that don't have exact MESH term matches. We want to filter for rows that don't have exact MESH term matches bc we already captured those and don't want to run scoring on it

    # -------    CONDITIONS    ------- #

    print("Use fuzzy matching on MeSH terms from Clinical Trials dump to find more potential matches.")
    conditions = df_dict["conditions"] 
    conditions = conditions[["nct_id", "downcase_name"]]

    browse_conditions = df_dict["browse_conditions"] 
    all_mesh_conditions = browse_conditions.downcase_mesh_term.unique()
    mask = np.isin(conditions['downcase_name'], all_mesh_conditions)
    conditions = conditions.assign(mesh_conditions_exact_mapped = np.where(mask, conditions['downcase_name'], np.nan))# interventions_unmapped = interventions[interventions['mesh_interventions_exact_mapped'].isnull()] # get the rows where mesh_term is empty bc there was no match there, we will run fuzzy scoring on these rows (I'm getting unmapped interventions along with their NCT IDs, not just the condition, bc I need this to find possible MeSH terms by study)
    conditions_unmapped = conditions[conditions['mesh_conditions_exact_mapped'].isnull()] # get the rows where mesh_term is empty bc there was no match there, we will run fuzzy scoring on these rows (I'm getting unmapped conditions along with their NCT IDs, not just the condition, bc I need this to find possible MeSH terms by study)
    conditions_unmapped = conditions_unmapped.drop('mesh_conditions_exact_mapped', axis=1) # drop the empty column now

    mesh_conditions_per_study = pd.DataFrame(browse_conditions[["nct_id", "downcase_mesh_term", "mesh_type"]].groupby("nct_id")["downcase_mesh_term"].apply(list)) # get all MeSH terms available for each study

    conditions_unmapped_all_mesh_terms = pd.merge(conditions_unmapped, 
                                                  mesh_conditions_per_study,
                                                  how='left',
                                                  left_on=['nct_id'],
                                                  right_on = ['nct_id'])

    # some clinical trials are missing from browse_conditions (those nct_ids are not present in the browse_conditions text) They have NaN in the downcase_mesh_term column
    conditions_unmapped_all_mesh_terms = conditions_unmapped_all_mesh_terms[~conditions_unmapped_all_mesh_terms['downcase_mesh_term'].isnull()] # subset or delete rows where either column is empty/Nonetype bc fuzzymatching functions will throw error if handling
    conditions_unmapped_all_mesh_terms = conditions_unmapped_all_mesh_terms[~conditions_unmapped_all_mesh_terms['downcase_name'].isnull()] # subset or delete rows where either column is empty/Nonetype bc fuzzymatching functions will throw error if handling

    conditions_unmapped_all_mesh_terms = conditions_unmapped_all_mesh_terms.explode('downcase_mesh_term')

    sort_ratio = np.vectorize(get_token_sort_ratio)
    set_ratio = np.vectorize(get_token_set_ratio)
    sim_score = np.vectorize(get_similarity_score)

    conditions_unmapped_all_mesh_terms["sort_ratio"] = sort_ratio(conditions_unmapped_all_mesh_terms[["downcase_mesh_term"]].values, conditions_unmapped_all_mesh_terms[["downcase_name"]].values) # generate fuzzy scores based between original and MeSH term
    conditions_unmapped_all_mesh_terms["sim_score"] = sim_score(conditions_unmapped_all_mesh_terms[["downcase_mesh_term"]].values, conditions_unmapped_all_mesh_terms[["downcase_name"]].values)
    conditions_mesh_fuzz_scored = conditions_unmapped_all_mesh_terms[(conditions_unmapped_all_mesh_terms['sim_score'] > 88) | (conditions_unmapped_all_mesh_terms['sort_ratio'] > 88)]
    conditions_mesh_fuzz_scored = conditions_mesh_fuzz_scored.sort_values(by = ['nct_id', 'downcase_name'], ascending = [True, True], na_position = 'first')

    conditions_mesh_fuzz_scored = conditions_mesh_fuzz_scored.sort_values(by = ['sim_score', 'sort_ratio'], ascending=[False,False], na_position='last').drop_duplicates(['nct_id', 'downcase_name']).sort_index() # there may be many mesh terms that passed the ratio and score filter; this causes duplicates bc I exploded the df...this line of code picks one and throws away other potential matches for one disease-nct_id pair
    conditions_mesh_fuzz_scored = conditions_mesh_fuzz_scored.sort_values(['nct_id'], ascending=False)

    keys = list(conditions_mesh_fuzz_scored[["nct_id", "downcase_name"]].columns.values)
    i1 = conditions_unmapped.set_index(keys).index
    i2 = conditions_mesh_fuzz_scored.set_index(keys).index

    tomap_conditions = conditions_mesh_fuzz_scored["downcase_mesh_term"].values.tolist()
    tomap_conditions = list(set(tomap_conditions))
    
    print("Since these MeSH terms don't come with identifers, we retrieve them from Name Resolver")
    mesh_fuzz_mapped_chunked = split_list(tomap_conditions, sublist_length)
    mesh_fuzz_mapped_curied = run_parallel_threads_nr(mesh_fuzz_mapped_chunked)
    mesh_fuzz_mapped_curied = [element for sublist in mesh_fuzz_mapped_curied for element in sublist] # flatten the list of lists
    print("Number of unique conditions mapped using fuzzy matching to MeSH terms: {}".format(len(mesh_fuzz_mapped_curied)))
    mesh_fuzz_mapped_curied = {key: value for dictionary in mesh_fuzz_mapped_curied for key, value in dictionary.items()}
    fuzz_mapped_conditions = pd.DataFrame({"condition_input": list(mesh_fuzz_mapped_curied.keys()),
                                           "condition_CURIE_id": [value[0] for value in mesh_fuzz_mapped_curied.values()],
                                           "condition_CURIE_name": [value[-1] for value in mesh_fuzz_mapped_curied.values()],
                                           "source": "MeSH term fuzzy mapped, Name Resolver CURIE"})
    previously_mapped = ct_terms["mapped_conditions"]
    combined_mapped_conditions = pd.concat([previously_mapped, fuzz_mapped_conditions], ignore_index=True) # get dataframe of combined previously mapped conditions and additional fuzzy MeSH mapped conditions

    all_conditions_list = conditions["downcase_name"].values.tolist()
    all_conditions_list = list(set(all_conditions_list))
    unmapped_conditions = list(set(all_conditions_list)-set(list(combined_mapped_conditions.condition_input.values)))
    print("Number of unique conditions that are unmapped after finding fuzzy MeSH mappings and using Name Resolver to get CURIES: {}".format(len(unmapped_conditions)))

    # -------    INTERVENTIONS    ------- #

    interventions = df_dict["interventions"] 
    interventions['downcase_name'] = interventions['name'].str.lower()
    interventions = interventions[["nct_id", "downcase_name"]]
    browse_interventions = df_dict["browse_interventions"] 
    all_mesh_interventions = browse_interventions.downcase_mesh_term.unique()
    mask = np.isin(interventions['downcase_name'], all_mesh_interventions)
    interventions = interventions.assign(mesh_interventions_exact_mapped = np.where(mask, interventions['downcase_name'], np.nan))# interventions_unmapped = interventions[interventions['mesh_interventions_exact_mapped'].isnull()] # get the rows where mesh_term is empty bc there was no match there, we will run fuzzy scoring on these rows (I'm getting unmapped interventions along with their NCT IDs, not just the condition, bc I need this to find possible MeSH terms by study)
    interventions_unmapped = interventions[interventions['mesh_interventions_exact_mapped'].isnull()] # get the rows where mesh_term is empty bc there was no match there, we will run fuzzy scoring on these rows (I'm getting unmapped interventions along with their NCT IDs, not just the condition, bc I need this to find possible MeSH terms by study)
    interventions_unmapped = interventions_unmapped.drop('mesh_interventions_exact_mapped', axis=1) # drop the empty column now
    mesh_interventions_per_study = pd.DataFrame(browse_interventions[["nct_id", "downcase_mesh_term", "mesh_type"]].groupby("nct_id")["downcase_mesh_term"].apply(list)) # get all MeSH terms available for each study
    interventions_unmapped_all_mesh_terms = pd.merge(interventions_unmapped, 
                                                     mesh_interventions_per_study,
                                                     how='left',
                                                     left_on=['nct_id'],
                                                     right_on = ['nct_id'])

    # # some clinical trials are missing from browse_interventions (those nct_ids are not present in the browse_interventions text) They have NaN in the downcase_mesh_term column
    interventions_unmapped_all_mesh_terms = interventions_unmapped_all_mesh_terms[~interventions_unmapped_all_mesh_terms['downcase_mesh_term'].isnull()] # subset or delete rows where either column is empty/Nonetype bc fuzzymatching functions will throw error if handling
    interventions_unmapped_all_mesh_terms = interventions_unmapped_all_mesh_terms[~interventions_unmapped_all_mesh_terms['downcase_name'].isnull()] # subset or delete rows where either column is empty/Nonetype bc fuzzymatching functions will throw error if handling

    interventions_unmapped_all_mesh_terms = interventions_unmapped_all_mesh_terms.explode('downcase_mesh_term')

    sort_ratio = np.vectorize(get_token_sort_ratio)
    set_ratio = np.vectorize(get_token_set_ratio)
    sim_score = np.vectorize(get_similarity_score)

    interventions_unmapped_all_mesh_terms["sort_ratio"] = sort_ratio(interventions_unmapped_all_mesh_terms[["downcase_mesh_term"]].values, interventions_unmapped_all_mesh_terms[["downcase_name"]].values) # generate fuzzy scores based between original and MeSH term
    interventions_unmapped_all_mesh_terms["sim_score"] = sim_score(interventions_unmapped_all_mesh_terms[["downcase_mesh_term"]].values, interventions_unmapped_all_mesh_terms[["downcase_name"]].values)
    interventions_mesh_fuzz_scored = interventions_unmapped_all_mesh_terms[(interventions_unmapped_all_mesh_terms['sim_score'] > 88) | (interventions_unmapped_all_mesh_terms['sort_ratio'] > 88)]
    interventiaons_mesh_fuzz_scored = interventions_mesh_fuzz_scored.sort_values(by = ['nct_id', 'downcase_name'], ascending = [True, True], na_position = 'first')

    interventions_mesh_fuzz_scored = interventions_mesh_fuzz_scored.sort_values(by = ['sim_score', 'sort_ratio'], ascending=[False,False], na_position='last').drop_duplicates(['nct_id', 'downcase_name']).sort_index() # there may be many mesh terms that passed the ratio and score filter; this causes duplicates bc I exploded the df...this line of code picks one and throws away other potential matches for one disease-nct_id pair
    interventions_mesh_fuzz_scored = interventions_mesh_fuzz_scored.sort_values(['nct_id'], ascending=False)

    keys = list(interventions_mesh_fuzz_scored[["nct_id", "downcase_name"]].columns.values)
    i1 = interventions_unmapped.set_index(keys).index
    i2 = interventions_mesh_fuzz_scored.set_index(keys).index

    tomap_interventions = interventions_mesh_fuzz_scored["downcase_mesh_term"].values.tolist()
    tomap_interventions = list(set(tomap_interventions))

    print("Since these MeSH terms don't come with identifers, we retrieve them from Name Resolver")
    mesh_fuzz_mapped_chunked = split_list(tomap_interventions, sublist_length)
    mesh_fuzz_mapped_curied = run_parallel_threads_nr(mesh_fuzz_mapped_chunked)
    mesh_fuzz_mapped_curied = [element for sublist in mesh_fuzz_mapped_curied for element in sublist] # flatten the list of lists
    print("Number of unique interventions mapped using fuzzy matching to MeSH terms: {}".format(len(mesh_fuzz_mapped_curied)))
    mesh_fuzz_mapped_curied = {key: value for dictionary in mesh_fuzz_mapped_curied for key, value in dictionary.items()}
    fuzz_mapped_interventions = pd.DataFrame({"intervention_input": list(mesh_fuzz_mapped_curied.keys()),
                                              "intervention_CURIE_id": [value[0] for value in mesh_fuzz_mapped_curied.values()],
                                              "intervention_CURIE_name": [value[-1] for value in mesh_fuzz_mapped_curied.values()],
                                              "source": "MeSH term fuzzy mapped, Name Resolver CURIE"})
    previously_mapped = ct_terms["mapped_interventions"]
    combined_mapped_interventions = pd.concat([previously_mapped, fuzz_mapped_interventions], ignore_index=True) # get dataframe of combined previously mapped interventions and additional fuzzy MeSH mapped interventions

    all_interventions_list = interventions["downcase_name"].values.tolist()
    all_interventions_list = list(set(all_interventions_list))
    unmapped_interventions = list(set(all_interventions_list)-set(list(combined_mapped_interventions.intervention_input.values)))
    print("Number of unique interventions that are unmapped after finding fuzzy MeSH mappings and using Name Resolver to get CURIES: {}".format(len(unmapped_interventions)))


    ct_terms = {'mapped_conditions': combined_mapped_conditions,
                'unmapped_conditions': unmapped_conditions,
                'mapped_interventions': combined_mapped_interventions,
                'unmapped_interventions': unmapped_interventions,
                "mesh_conditions_per_study": mesh_conditions_per_study,
                "mesh_interventions_per_study": mesh_interventions_per_study}
    return ct_terms

def term_list_to_nr(df_dict, ct_terms):
    
    # -------    CONDITIONS    ------- #

    unmapped_conditions = ct_terms["unmapped_conditions"]
    conditions_unmapped_chunked = split_list(unmapped_conditions, sublist_length)
    nr_conditions = run_parallel_threads_nr(conditions_unmapped_chunked)
    nr_conditions = [element for sublist in nr_conditions for element in sublist] # flatten the list of lists
    print("Number of unique conditions mapped using Name Resolver: {}".format(len(nr_conditions)))
    nr_conditions = {key: value for dictionary in nr_conditions for key, value in dictionary.items()}
    nr_conditions_df = pd.DataFrame({"condition_input": list(nr_conditions.keys()),
                                     "condition_CURIE_id": [value[0] for value in nr_conditions.values()],
                                     "condition_CURIE_name": [value[-1] for value in nr_conditions.values()],
                                     "source": "Name Resolver, no further curation"})
    
    previously_mapped = ct_terms["mapped_conditions"]
    combined_mapped_conditions = pd.concat([previously_mapped, nr_conditions_df], ignore_index=True) # get dataframe of combined previously mapped conditions and additional fuzzy MeSH mapped conditions
    
    conditions = df_dict["conditions"]
    all_conditions_list = conditions["downcase_name"].values.tolist()
    all_conditions_list = list(set(all_conditions_list))
    unmapped_conditions = list(set(all_conditions_list)-set(list(combined_mapped_conditions.condition_input.values)))
    print("Number of unique conditions that are unmapped after using Name Resolver: {}".format(len(unmapped_conditions)))
    
    # -------    INTERVENTIONS    ------- #
    
    unmapped_interventions = ct_terms["unmapped_interventions"]
    interventions_unmapped_chunked = split_list(unmapped_interventions, sublist_length)
    nr_interventions = run_parallel_threads_nr(interventions_unmapped_chunked)
    nr_interventions = [element for sublist in nr_interventions for element in sublist] # flatten the list of lists
    print("Number of unique interventions mapped using Name Resolver: {}".format(len(nr_interventions)))
    nr_interventions = {key: value for dictionary in nr_interventions for key, value in dictionary.items()}
    nr_interventions_df = pd.DataFrame({"intervention_input": list(nr_interventions.keys()),
                                     "intervention_CURIE_id": [value[0] for value in nr_interventions.values()],
                                     "intervention_CURIE_name": [value[-1] for value in nr_interventions.values()],
                                     "source": "Name Resolver, no further curation"})
    
    previously_mapped = ct_terms["mapped_interventions"]
    combined_mapped_interventions = pd.concat([previously_mapped, nr_interventions_df], ignore_index=True) # get dataframe of combined previously mapped interventions and additional fuzzy MeSH mapped interventions
    
    interventions = df_dict["interventions"]
    all_interventions_list = interventions["downcase_name"].values.tolist()
    all_interventions_list = list(set(all_interventions_list))
    unmapped_interventions = list(set(all_interventions_list)-set(list(combined_mapped_interventions.intervention_input.values)))
    print("Number of unique interventions that are unmapped after using Name Resolver: {}".format(len(unmapped_interventions)))
    
    ct_terms = {'mapped_conditions': combined_mapped_conditions, 'unmapped_conditions': unmapped_conditions, 'mapped_interventions': combined_mapped_interventions, 'unmapped_interventions': unmapped_interventions}
    return ct_terms



def get_service_ticket(serverurl, ticket_granting_ticket, serviceurl):
    """ Obtain a Single-Use Proxy Ticket (also known as service ticket).
    Request for a Service Ticket:
        POST /cas/v1/tickets/{TGT id} HTTP/1.0
    data:
           service={form encoded parameter for the service url}
    Sucessful Response:
        200 OK
        ST-1-FFDFHDSJKHSDFJKSDHFJKRUEYREWUIFSD2132
    @param serverurl authentication server
    @param ticketGrantingTicket a Proxy Granting Ticket.
    @param serviceurl url of service with protected resources
    @return authentication ticket for service. """
    resp = requests.post("{}/{}".format(serverurl, ticket_granting_ticket),
                         {"service": serviceurl})
    if resp.status_code == 200:
        return resp.content
    return 'Error: status: {}'.format(resp.content)

def get_ticket(cas_serverurl, apikey, serviceurl):
    # set ticket granting ticket server url
    tgtserverurl = cas_serverurl + "/api-key"
    # set service ticket server url
    stserverurl = cas_serverurl + "/tickets"
    tgt = get_ticket_granting_ticket(tgtserverurl, apikey)
    return get_service_ticket(stserverurl, tgt, serviceurl)

def get_ticket_granting_ticket(tgtserverurl, apikey):
    # http://serviceurl/cas/v1/tickets/{TGT id}
    response = requests.post(tgtserverurl, {'apikey': apikey},
                             headers={'Accept': 'test/plain'})
    return extract_tgt_ticket(response.content)

def extract_tgt_ticket(htmlcontent):
#     print(htmlcontent)
    "Extract ticket granting ticket from HTML."    
    soup = BeautifulSoup(htmlcontent)
#     print(soup.find('form').get("action"))
    cas_url = soup.find("form").get("action")
    "Extract ticket granting ticket out of 'action' attribute"
    return cas_url.rsplit('/')[-1]

def get_redirect_target(resp):
        """Receives a Response. Returns a redirect URI or ``None``"""
        # Due to the nature of how requests processes redirects this method will
        # be called at least once upon the original response and at least twice
        # on each subsequent redirect response (if any).
        # If a custom mixin is used to handle this logic, it may be advantageous
        # to cache the redirect location onto the response object as a private
        # attribute.
        if resp.is_redirect:
            location = resp.headers["location"]
            # Currently the underlying http module on py3 decode headers
            # in latin1, but empirical evidence suggests that latin1 is very
            # rarely used with non-ASCII characters in HTTP headers.
            # It is more likely to get UTF8 header rather than latin1.
            # This causes incorrect handling of UTF8 encoded location headers.
            # To solve this, we re-encode the location in latin1.
#             print(location)
            location = location.encode("latin1")
#             print(location)
#             print(to_native_string(location, "utf8"))
            return to_native_string(location, "utf8")
        return None


def get_metamap_mappings(chunk, args):
    
    form = {}
    form['KSOURCE'] = ksource
    form['COMMAND_ARGS'] = args
    headers = {'Accept': 'application/json'}

    mm_terms = {}
    cui_pattern = r"C\d+(?=:)"
    name_pattern = r"(?<=:)[^[]+"
    semtype_pattern = r"\[(.*?)\]"
    
    form['APIText'] = chunk
    service_ticket = get_ticket(CAS_SERVERURL, apikey, serviceurl)
    params = {'ticket': service_ticket}

    s = requests.Session()
    trycnt = 5  # max try count to receive response from MetaMap Interactive API
    while trycnt > 0:
        try:
            response = s.post(serviceurl, form, headers=headers, params=params, allow_redirects=False)
            if response.status_code == 302:
                newurl = s.get_redirect_target(response)
                response = s.post(newurl, form, headers=headers, params=params, allow_redirects=False)
            trycnt = 0 # success, recieved response from MetaMap Interactive Server
        except (ConnectionResetError,OSError) as ex:  # Catch ProtocolError or socket.error in requests that raises a ConnectionError as "OSError" ....https://stackoverflow.com/questions/74253820/cannot-catch-requests-exceptions-connectionerror-with-try-except
            if trycnt <= 0: print("Failed to retrieve MetaMap response\n" + str(ex))  # done retrying
            else: trycnt -= 1  # retry
            time.sleep(5)  # wait 5 seconds, then retry
            
    for line in response.text.splitlines():
        if not any(s in line for s in ["Meta Mapping", "Processing", "/dmzfiler/"]):
            if "Phrase:" in line:
                cuis_per_input = []
                mm_input = line.split(":")[1].strip()
                cui_match_count = 0
            else:
                cui_match = re.findall(cui_pattern, line)
                if cui_match: 
                    cui_match_count +=1
                    if cui_match_count > 1: # get only 1st CUI/CURIE per Phrase; continue to next loop iteration to skip lines with more available CUIs
                        continue
                    cui_info = []
                    name_match = re.findall(name_pattern, line)
                    semtype_match = re.findall(semtype_pattern, line)
                    try: cui_info.append(cui_match[0].strip())
                    except: cui_info.append(None)
                    try: cui_info.append(name_match[0].strip())
                    except: cui_info.append(None)
                    try: cui_info.append(semtype_match[0].strip())
                    except: cui_info.append(None)
                    cuis_per_input.append(cui_info)
                    mm_terms[mm_input] = cuis_per_input
    return mm_terms


def run_parallel_threads_mm(terms_chunked, args):
    # multithread implementation for retrieving MetaMap API responses
    # Create a ThreadPoolExecutor with the desired number of threads
    with concurrent.futures.ThreadPoolExecutor(multiprocessing.cpu_count() - 1) as executor:
        # Submit the get_response() function for each item in the list
        futures = [executor.submit(get_metamap_mappings, term, args) for term in terms_chunked]

        # Retrieve the results as they become available
        output = [future.result() for future in concurrent.futures.as_completed(futures)]
    # mm_dict = reduce(lambda d1, d2: {**d1, **d2}, output) # merge the list of dicts of MetaMap responses in output into 1 dict
    mm_dict = reduce(lambda d1, d2: {**d1, **d2}, output, {})
    return mm_dict

def split_list_by_char_lim(lst):
    result = []
    current_sublist = []
    current_length = 0
    for item in lst:
        item_length = len(item)
        if current_length + item_length > 9500: # max is 10,000 char allowed by MetaMap
            result.append(current_sublist)
            current_sublist = []
            current_length = 0
        item = item + "\n"  # add a "\n" for term processing option by MetaMap, the terms in the input file must be separated by blank lines (https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/TermProcessing.pdf)
        current_sublist.append(item)
        current_length += item_length
    result.append(current_sublist)
    return result


def term_list_to_mm(df_dict, ct_terms):
    
    # -------    CONDITIONS    ------- #
    print("Using UMLS MetaMap to get more mappings for conditions. MetaMap returns mappings, CUIs, and semantic type of mapping.")
    unmapped_conditions = ct_terms["unmapped_conditions"]
    conditions_unmapped_chunked = split_list_by_char_lim(unmapped_conditions)
    # see MetaMap Usage instructions: https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/MM_2016_Usage.pdf
    # removing sosy semantic type (sign or symptom) - often get MetaMap matches to the sign or symptom instead of the full disease...for example, will get back "exercise-induced" instead of "immune dysfunction" for "exercise-induced immune dysfunction" bc it matches the descriptive quality "exercise-induced" is matched on 
    condition_args = ['--sldi -I -C -J acab,anab,cgab,comd,dsyn,inpo,mobd,neop,patf,sosy -z -i -f']  # see https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/SemanticTypes_2018AB.txt for semantic types ("acab,anab,etc.")
    mm_conditions = run_parallel_threads_mm(conditions_unmapped_chunked, condition_args)
    flattened_mm_conditions = {key: [item for sublist in value for item in sublist] for key, value in mm_conditions.items()}
    mm_conditions_df = pd.DataFrame({"condition_input": list(flattened_mm_conditions.keys()),
                                     "condition_CURIE_id": [value[0] for value in flattened_mm_conditions.values()],
                                     "condition_CURIE_name": [value[1] for value in flattened_mm_conditions.values()],
                                     "condition_semantic_type": [value[-1] for value in flattened_mm_conditions.values()],
                                     "source": "MetaMap via UMLS, term and CURIE"})
    
    mm_conditions_df[['condition_CURIE_name_1', 'condition_CURIE_name_2']] = mm_conditions_df['condition_CURIE_name'].str.extract(r'^(.*?)\s*\((.*?)\)$').fillna('NA') # 

    sort_ratio = np.vectorize(get_token_sort_ratio)
    set_ratio = np.vectorize(get_token_set_ratio)
    sim_score = np.vectorize(get_similarity_score)

    # many MetaMap terms are returned as "term (term)". For example, "Nonessential Amino Acid (Nonessential amino acid)". This repetition messes up the sort ratio and sim score, so we extract the substrings out of the parenthesis to conduct scoring on those
    mm_conditions_scored = mm_conditions_df.copy()
    mm_conditions_scored["sort_ratio"] = sort_ratio(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name"]].values) # generate fuzzy scores based between original and MeSH term
    mm_conditions_scored["sim_score"] = sim_score(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name"]].values)

    mm_conditions_scored["sort_ratio_1"] = sort_ratio(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name_1"]].values) # generate fuzzy scores based between original and MetaMap term
    mm_conditions_scored["sim_score_1"] = sim_score(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name_1"]].values)

    mm_conditions_scored["sort_ratio_2"] = sort_ratio(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name_2"]].values) # generate fuzzy scores based between original and MetaMap term
    mm_conditions_scored["sim_score_2"] = sim_score(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name_2"]].values)

    mm_conditions_scored_thresholded = mm_conditions_scored.copy() 
    
    mm_conditions_scored["sort_ratio"] = sort_ratio(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name"]].values) # generate fuzzy scores based between original and MetaMap term
    mm_conditions_scored["sim_score"] = sim_score(mm_conditions_scored[["condition_input"]].values, mm_conditions_scored[["condition_CURIE_name"]].values)
    mm_conditions_scored_thresholded = mm_conditions_scored.copy() 
    mm_conditions_scored_thresholded = mm_conditions_scored_thresholded[(mm_conditions_scored_thresholded['sim_score'] > 88) |
                                                                        (mm_conditions_scored_thresholded['sort_ratio'] > 88) |
                                                                        (mm_conditions_scored_thresholded['sim_score_1'] > 88) |
                                                                        (mm_conditions_scored_thresholded['sort_ratio_1'] > 88) |
                                                                        (mm_conditions_scored_thresholded['sim_score_2'] > 88) |
                                                                        (mm_conditions_scored_thresholded['sort_ratio_2'] > 88)]
    
    print("Number of unique conditions that are mapped after using MetaMap and similarity and ratio score thresholds of 88: {}".format(mm_conditions_scored_thresholded.shape[0]))
    
    mm_conditions_scored_thresholded = mm_conditions_scored_thresholded.drop(['condition_CURIE_name_1',
                                                                              'condition_CURIE_name_2',
                                                                              'sort_ratio',
                                                                              'sim_score',
                                                                              'sort_ratio_1',
                                                                              'sim_score_1',
                                                                              'sort_ratio_2',
                                                                              'sim_score_2'], axis=1)
    previously_mapped = ct_terms["mapped_conditions"]
    combined_mapped_conditions = pd.concat([previously_mapped, mm_conditions_scored_thresholded], ignore_index=True) # get dataframe of combined previously mapped conditions and additional MetaMapped interventions that passed threshold scoring

    conditions = df_dict["conditions"]
    all_conditions_list = conditions["downcase_name"].values.tolist()
    all_conditions_list = list(set(all_conditions_list))
    unmapped_conditions = list(set(all_conditions_list)-set(list(combined_mapped_conditions.condition_input.values)))
    print("Number of unique conditions that are unmapped after using MetaMap and similarity and ratio score thresholds of 88: {}".format(len(unmapped_conditions)))
          
    # -------    INTERVENTIONS    ------- #
    print("Using UMLS MetaMap to get more mappings for interventions. MetaMap returns mappings, CUIs, and semantic type of mapping.")
    unmapped_interventions = ct_terms["unmapped_interventions"]
    interventions_unmapped_chunked = split_list_by_char_lim(unmapped_interventions)
    # see MetaMap Usage instructions: https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/MM_2016_Usage.pdf
    # removing sosy semantic type (sign or symptom) - often get MetaMap matches to the sign or symptom instead of the full disease...for example, will get back "exercise-induced" instead of "immune dysfunction" for "exercise-induced immune dysfunction" bc it matches the descriptive quality "exercise-induced" is matched on 
    intervention_args = ['--sldi -I -C -k acab,anab,cgab,comd,dsyn,inpo,mobd,neop,patf,sosy -z -i -f']  # see https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/SemanticTypes_2018AB.txt for semantic types ("acab,anab,etc.") (I used inverse of semantic terms picked for conditions here)
    mm_interventions = run_parallel_threads_mm(interventions_unmapped_chunked, intervention_args)
    flattened_mm_interventions = {key: [item for sublist in value for item in sublist] for key, value in mm_interventions.items()}
    mm_interventions_df = pd.DataFrame({"intervention_input": list(flattened_mm_interventions.keys()),
                                        "intervention_CURIE_id": [value[0] for value in flattened_mm_interventions.values()],
                                        "intervention_CURIE_name": [value[1] for value in flattened_mm_interventions.values()],
                                        "intervention_semantic_type": [value[-1] for value in flattened_mm_interventions.values()],
                                        "source": "MetaMap via UMLS, term and CURIE"})

    mm_interventions_df[['intervention_CURIE_name_1', 'intervention_CURIE_name_2']] = mm_interventions_df['intervention_CURIE_name'].str.extract(r'^(.*?)\s*\((.*?)\)$').fillna('NA') # 

    sort_ratio = np.vectorize(get_token_sort_ratio)
    set_ratio = np.vectorize(get_token_set_ratio)
    sim_score = np.vectorize(get_similarity_score)

    # many MetaMap terms are returned as "term (term)". For example, "Nonessential Amino Acid (Nonessential amino acid)". This repetition messes up the sort ratio and sim score, so we extract the substrings out of the parenthesis to conduct scoring on those
    mm_interventions_scored = mm_interventions_df.copy()
    mm_interventions_scored["sort_ratio"] = sort_ratio(mm_interventions_scored[["intervention_input"]].values, mm_interventions_scored[["intervention_CURIE_name"]].values) # generate fuzzy scores based between original and MeSH term
    mm_interventions_scored["sim_score"] = sim_score(mm_interventions_scored[["intervention_input"]].values, mm_interventions_scored[["intervention_CURIE_name"]].values)

    mm_interventions_scored["sort_ratio_1"] = sort_ratio(mm_interventions_scored[["intervention_input"]].values, mm_interventions_scored[["intervention_CURIE_name_1"]].values) # generate fuzzy scores based between original and MetaMap term
    mm_interventions_scored["sim_score_1"] = sim_score(mm_interventions_scored[["intervention_input"]].values, mm_interventions_scored[["intervention_CURIE_name_1"]].values)

    mm_interventions_scored["sort_ratio_2"] = sort_ratio(mm_interventions_scored[["intervention_input"]].values, mm_interventions_scored[["intervention_CURIE_name_2"]].values) # generate fuzzy scores based between original and MetaMap term
    mm_interventions_scored["sim_score_2"] = sim_score(mm_interventions_scored[["intervention_input"]].values, mm_interventions_scored[["intervention_CURIE_name_2"]].values)

    mm_interventions_scored_thresholded = mm_interventions_scored.copy() 
    mm_interventions_scored_thresholded = mm_interventions_scored_thresholded[(mm_interventions_scored_thresholded['sim_score'] > 88) |
                                                                              (mm_interventions_scored_thresholded['sort_ratio'] > 88) |
                                                                              (mm_interventions_scored_thresholded['sim_score_1'] > 88) |
                                                                              (mm_interventions_scored_thresholded['sort_ratio_1'] > 88) |
                                                                              (mm_interventions_scored_thresholded['sim_score_2'] > 88) |
                                                                              (mm_interventions_scored_thresholded['sort_ratio_2'] > 88)]
    
    print("Number of unique interventions that are mapped after using MetaMap and similarity and ratio score thresholds of 88: {}".format(mm_interventions_scored_thresholded.shape[0]))

    mm_interventions_scored_thresholded = mm_interventions_scored_thresholded.drop(['intervention_CURIE_name_1',
                                                                                    'intervention_CURIE_name_2',
                                                                                    'sort_ratio',
                                                                                    'sim_score',
                                                                                    'sort_ratio_1',
                                                                                    'sim_score_1',
                                                                                    'sort_ratio_2',
                                                                                    'sim_score_2'], axis=1)
    previously_mapped = ct_terms["mapped_interventions"]
    combined_mapped_interventions = pd.concat([previously_mapped, mm_interventions_scored_thresholded], ignore_index=True) # get dataframe of combined previously mapped interventions and additional MetaMapped interventions that passed threshold scoring
    interventions = df_dict["interventions"]
    all_interventions_list = interventions["downcase_name"].values.tolist()
    all_interventions_list = list(set(all_interventions_list))
    unmapped_interventions = list(set(all_interventions_list)-set(list(combined_mapped_interventions.intervention_input.values)))
    print("Number of unique interventions that are unmapped after using MetaMap and similarity and ratio score thresholds of 88: {}".format(len(unmapped_interventions)))
    ct_terms = {'mapped_conditions': combined_mapped_conditions,
                'unmapped_conditions': unmapped_conditions,
                'mapped_interventions': combined_mapped_interventions,
                'unmapped_interventions': unmapped_interventions,
                'all_metamapped_conditions': mm_conditions_df,
                'all_metamapped_interventions': mm_interventions_df}


    return ct_terms

# output all results to TSVs
def compile_and_output(df_dict, ct_terms, remaining_unmapped_possible):
    print("\n")
    print("#   -------- -------- -------- --------  ")
    print("Final Tallies:")
    print("Total # of conditions mapped: {}".format(ct_terms["mapped_conditions"].shape[0]))
    print("Total # of interventions mapped: {}".format(ct_terms["mapped_interventions"].shape[0]))
    print("Total # of conditions unmapped or not mapped: {}".format(len(ct_terms["unmapped_conditions"])))
    print("Total # of interventions unmapped or not mapped: {}".format(len(ct_terms["unmapped_interventions"])))    
    # How many Clinical Trials are there? Well, it's different depending on the Conditions or Interventions dataframes...
    conditions_nctids = len(df_dict["conditions"].nct_id.unique())
    interventions_nctids = len(df_dict["interventions"].nct_id.unique())
    print("Number of Clinical Trials NCITs in Conditions table: {}".format(conditions_nctids))      
    print("Number of Clinical Trials NCITs in Interventions table: {}".format(interventions_nctids))
    print("#   -------- -------- -------- --------  ")

    """ create tables of unused MeSH and MetaMap CURIEs that could be used for unmapped Conditions and Interventions """
    # -------    CONDITIONS    ------- #
    all_conditions = df_dict["conditions"][["nct_id", "downcase_name"]]
    conditions_mesh = pd.merge(all_conditions, 
                               mesh_conditions_per_study,
                               how='left',
                               left_on=['nct_id'],
                               right_on = ['nct_id'])
    
    metamap_possibilities = remaining_unmapped_possible["all_metamapped_conditions"][["condition_input", "condition_CURIE_id", "condition_CURIE_name", "condition_semantic_type"]]
    conditions_mesh_metamap = pd.merge(conditions_mesh, 
                                       metamap_possibilities,
                                       how='left',
                                       left_on=['downcase_name'],
                                       right_on = ['condition_input'])
    
    unmapped_conditions_possible_terms = conditions_mesh_metamap[conditions_mesh_metamap['downcase_name'].isin(ct_terms["unmapped_conditions"])]
    unmapped_conditions_possible_terms = unmapped_conditions_possible_terms.drop('condition_input', axis=1) # drop the redundant column now
    
    # -------    INTERVENTIONS    ------- #
    all_interventions = df_dict["interventions"][["nct_id", "downcase_name"]]
    interventions_mesh = pd.merge(all_interventions, 
                               mesh_interventions_per_study,
                               how='left',
                               left_on=['nct_id'],
                               right_on = ['nct_id'])
    
    metamap_possibilities = remaining_unmapped_possible["all_metamapped_interventions"][["intervention_input", "intervention_CURIE_id", "intervention_CURIE_name", "intervention_semantic_type"]]
    interventions_mesh_metamap = pd.merge(interventions_mesh, 
                                       metamap_possibilities,
                                       how='left',
                                       left_on=['downcase_name'],
                                       right_on = ['intervention_input'])
    
    unmapped_interventions_possible_terms = interventions_mesh_metamap[interventions_mesh_metamap['downcase_name'].isin(ct_terms["unmapped_interventions"])]
    unmapped_interventions_possible_terms = unmapped_interventions_possible_terms.drop('intervention_input', axis=1) # drop the redundant column now
          
        
    """   Output all to TSVs   """    
    pd.Series(ct_terms["unmapped_conditions"]).to_csv('unmapped_conditions.tsv', sep="\t", index=False, header=False) # convert the list to a pandas series, then output to TSV
    pd.Series(ct_terms["unmapped_interventions"]).to_csv('unmapped_interventions.tsv', sep="\t", index=False, header=False) # convert the list to a pandas series, then output to TSV
    ct_terms["mapped_conditions"].to_csv('mapped_conditions.tsv', sep="\t", index=False)
    ct_terms["mapped_interventions"].to_csv('mapped_interventions.tsv', sep="\t", index=False)
    unmapped_conditions_possible_terms.to_csv('unmapped_conditions_possible_mappings.tsv', sep="\t", index=False)
    unmapped_interventions_possible_terms.to_csv('unmapped_interventions_possible_mappings.tsv', sep="\t", index=False)
    

def main():
    # flag_and_path = get_raw_ct_data() # uncomment for production
    # flag_and_path = {'term_program_flag': False, 'data_extracted_path': '/Users/Kamileh/Work/ISB/NCATS_BiomedicalTranslator/Projects/ClinicalTrials/ETL_Python/data/07_25_2023_extracted'} # comment for production
    flag_and_path = {'term_program_flag': False, 'data_extracted_path': './data/07_25_2023_extracted'} # comment for production
    df_dict = read_raw_ct_data(flag_and_path)
    ct_terms = exact_match_mesh(df_dict)
    ct_terms = inexact_match_mesh(df_dict, ct_terms)

    # pull the available MeSH terms per study out of the returned ct_terms dict 
    mesh_conditions_per_study = ct_terms["mesh_conditions_per_study"]
    mesh_interventions_per_study = ct_terms["mesh_interventions_per_study"]

    ct_terms = term_list_to_nr(df_dict, ct_terms)
    ct_terms = term_list_to_mm(df_dict, ct_terms)

    # pull the available UMLS terms per study out of the returned ct_terms dict 
    all_metamapped_conditions = ct_terms["all_metamapped_conditions"]
    all_metamapped_interventions = ct_terms["all_metamapped_interventions"]

    remaining_unmapped_possible = {"mesh_conditions_per_study": mesh_conditions_per_study,
                                   "mesh_interventions_per_study": mesh_interventions_per_study,
                                   "all_metamapped_conditions": all_metamapped_conditions,
                                   "all_metamapped_interventions": all_metamapped_interventions}

    compile_and_output(df_dict, ct_terms, remaining_unmapped_possible)


if __name__ == "__main__":
    main()





