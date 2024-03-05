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
from time import sleep
import concurrent.futures
import multiprocessing
import datetime as dt
from datetime import date
import pathlib
import configparser
import sys
import urllib
import zipfile
import csv
import gc
# sys.path.insert(0, '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/pymetamap-master') # for local
sys.path.insert(0, '/users/knarsinh/projects/clinical_trials/metamap/pymetamap') # for hypatia
from pymetamap import MetaMap  # https://github.com/AnthonyMRios/pymetamap/blob/master/pymetamap/SubprocessBackend.py
from pandas import ExcelWriter
import ast
import glob
from tqdm import tqdm
import subprocess
import shlex
from collections import Counter
from ratelimit import limits, sleep_and_retry
import threading
from threading import Thread
from joblib import Parallel, delayed
csv_writer_lock = threading.Lock()

# %pip install thefuzz
# %pip install levenshtein
# %pip install xlsxwriter
# %pip install ratelimit
# %pip install timeout_decorator

from thefuzz import fuzz # fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python

# 40 calls per minute
CALLS = 40
RATE_LIMIT = 60


def get_token_sort_ratio(str1, str2):
    """ fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python """
    try:
        return fuzz.token_sort_ratio(str1, str2)
    except:
        return None

def get_similarity_score(str1, str2):
    """ fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python """
    try:
        return fuzz.ratio(str1, str2)
    except:
        return None
    
def convert_seconds_to_hms(seconds):
    """ converts the elapsed time or runtime to hours, min, sec """
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return hours, minutes, seconds

def de_ascii_er(text):
    non_ascii = "[^\x00-\x7F]"
    pattern = re.compile(r"[^\x00-\x7F]")
    non_ascii_text = re.sub(pattern, ' ', text)
    return non_ascii_text

# def start_metamap_servers(metamap_base_dir, metamap_bin_dir):
#     global metamap_pos_server_dir
#     global metamap_wsd_server_dir
#     metamap_pos_server_dir = 'bin/skrmedpostctl' # Part of speech tagger
#     metamap_wsd_server_dir = 'bin/wsdserverctl' # Word sense disambiguation 
    
#     metamap_executable_path_pos = os.path.join(metamap_base_dir, metamap_pos_server_dir)
#     metamap_executable_path_wsd = os.path.join(metamap_base_dir, metamap_wsd_server_dir)
#     command_pos = [metamap_executable_path_pos, 'start']
#     command_wsd = [metamap_executable_path_wsd, 'start']

#     # Start servers, with open portion redirects output of metamap server printing output to NULL
#     with open(os.devnull, "w") as fnull:
#         result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
#         result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
#     sleep(5)

def start_metamap_servers(metamap_dirs):
    global metamap_pos_server_dir
    global metamap_wsd_server_dir
    metamap_pos_server_dir = 'bin/skrmedpostctl' # Part of speech tagger
    metamap_wsd_server_dir = 'bin/wsdserverctl' # Word sense disambiguation 
    
    metamap_base_dir = metamap_dirs["metamap_base_dir"]

    metamap_executable_path_pos = os.path.join(metamap_base_dir, metamap_pos_server_dir)
    metamap_executable_path_wsd = os.path.join(metamap_base_dir, metamap_wsd_server_dir)
    command_pos = [metamap_executable_path_pos, 'start']
    command_wsd = [metamap_executable_path_wsd, 'start']

    # Start servers, with open portion redirects output of metamap server printing output to NULL
    with open(os.devnull, "w") as fnull:
        result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
        result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
    sleep(5)    

def stop_metamap_servers(metamap_dirs):

    metamap_base_dir = metamap_dirs["metamap_base_dir"]

    metamap_executable_path_pos = os.path.join(metamap_base_dir, metamap_pos_server_dir)
    metamap_executable_path_wsd = os.path.join(metamap_base_dir, metamap_wsd_server_dir)
    command_pos = [metamap_executable_path_pos, 'stop']
    command_wsd = [metamap_executable_path_wsd, 'stop']
    
    # Stop servers, with open portion redirects output of metamap server printing output to NULL
    with open(os.devnull, "w") as fnull:
        result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
        result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
    sleep(5)  
    
def check_os():
    if "linux" in sys.platform:
        print("Linux platform detected")
        metamap_base_dir = "/users/knarsinh/projects/clinical_trials/metamap/public_mm/"    # /users/knarsinh/projects/clinical_trials/metamap/public_mm 
        metamap_bin_dir = 'bin/metamap20'
        metamap_version = '2020'
    else:
        metamap_base_dir = '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/public_mm/' # for running on local
        metamap_bin_dir = 'bin/metamap18'
        metamap_version = '2018'
        
    return {"metamap_base_dir":metamap_base_dir, "metamap_bin_dir":metamap_bin_dir, "metamap_version":metamap_version} 

@sleep_and_retry
@limits(calls=CALLS, period=RATE_LIMIT)
def check_limit():
    ''' Empty function just to check for calls to Name Resolver API '''
    return

def wrap(x): # use this to convert string objects to dicts 
    try:
        a = ast.literal_eval(x)
        return(a)
    except:
        pass


def get_raw_ct_data():
    term_program_flag = True
    global data_dir
    global data_extracted
    
    try:
        # get all the links and associated dates of upload into a dict called date_link
        url_all = "https://aact.ctti-clinicaltrials.org/download"
        response = requests.get(url_all)
        soup = BeautifulSoup(response.text, features="lxml")
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
    except:
        print("continue")

    if not os.path.exists(data_path):   # if folder containing most recent data doesn't exist, download and extract it into data folder

        term_program_flag = False   # flag below for terminating program if latest download exists (KG is assumed up to date)
        print("Attempting download of Clinical Trial data as of {}\n".format(date_string))
        try:
            response = requests.get(url)
            if response.status_code == 200:
                with open(data_path, 'wb') as file:
                    file.write(response.content)
                print("Finished download of zip")
                with zipfile.ZipFile(data_path, 'r') as download:
                    print("Unzipping data")
                    download.extractall(data_extracted)
        except:
            print("Failed to scrape AACT for download. Please navigate to https://aact.ctti-clinicaltrials.org/download and manually download zip file.")
            print("Please store the downloaded zip in the /data directory. This should be the only item besides the cache file, condition manual review file, and intervention manual review file, in the directory at this time.")
            done = input("Type Done when done: ")
            if done == "Done":
                data_dir = "{}/data".format(pathlib.Path.cwd())
                # list_of_files = glob.glob(data_dir + "/*") # get all files in directory
                try:
                    # latest_file = max(list_of_files, key=os.path.getctime) # get the most recent file in the directory
                    pattern = os.path.join(data_dir, "*.zip")
                    zip_file = glob.glob(pattern) # look for file in directory that ends in ".zip"
                    zip_file = zip_file[0]
                    print("File found at: ")
                    print(zip_file)
                    # print(latest_file)
                    print("Please make sure this the correct zip file from AACT")
                    if not os.path.exists(data_extracted):   # if folder of unzipped data does not exist, unzip
                        try:
                            with zipfile.ZipFile(zip_file, 'r') as download:
                                print("Unzipping data into")
                                cttime = os.path.getctime(zip_file)
                                date_string = dt.datetime.fromtimestamp(cttime).strftime('%m_%d_%Y')
                                data_extracted = data_dir + "/{}_extracted".format(date_string)
                                print(data_extracted)
                                download.extractall(data_extracted)
                        except:
                            pattern = os.path.join(data_dir, "*_extracted")
                            extracted_file = glob.glob(pattern) # look for file in directory that ends in "_extracted"
                            data_extracted = extracted_file[0]
                            extracted_name = os.path.basename(os.path.normpath(extracted_file[0]))
                            date_string = extracted_name.replace('_extracted', '')
                            print("Assuming data is already unzipped")
                        
                except:
                    print("Unable to download and extract Clincal Trial data.")
                    print("Cannot find pipe-delimited zip in /data folder.")
    else:
        print("KG is already up to date.")

    return {"term_program_flag": term_program_flag, "data_extracted_path": data_extracted, "date_string": date_string}

def read_raw_ct_data(flag_and_path, subset_size):
    if flag_and_path["term_program_flag"]:
        print("Exiting program. Assuming KG has already been constructed from most recent data dump from AACT.")
        exit()
    else:
        data_extracted = flag_and_path["data_extracted_path"]
        # read in pipe-delimited files 
        conditions_df = pd.read_csv(data_extracted + '/conditions.txt.gz', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
        interventions_df = pd.read_csv(data_extracted + '/interventions.txt.gz', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
        interventions_alts_df = pd.read_csv(data_extracted + '/intervention_other_names.txt.gz', sep='|', index_col=False, header=0, on_bad_lines = 'warn')

        if subset_size:   # if a subset size is given, we are running this script on a small subset of the dataset
            conditions_df = conditions_df.sample(n=subset_size)
            interventions_df = interventions_df.sample(n=subset_size)
            interventions_alts_df = interventions_alts_df.sample(n=subset_size)
    
    df_dict = {"conditions": conditions_df, "interventions": interventions_df, "interventions_alts": interventions_alts_df}
    return df_dict



def cache_manually_selected_terms():    

    def return_curie_dict(curie_info_delimited):
        keys = ["mapped_name", "mapped_curie", "mapped_score", "mapped_semtypes"]
        curie_list = curie_info_delimited.split(" | ")
        curie_dict = dict(zip(keys, curie_list))
        return curie_dict

    files = glob.glob("*.xlsx")
    manually_selected_file = [i for i in files if "manual_review" in i if not i.startswith("~")][0] # find the file of manual selections
    manually_selected = pd.read_excel(manually_selected_file)
    cols_to_fill = ["mapping_tool", "term_type", "clintrial_term", "input_term"]
    manually_selected.loc[:,cols_to_fill] = manually_selected.loc[:,cols_to_fill].ffill()

    manually_selected = manually_selected[~manually_selected['manually_selected_CURIE'].isnull()] # get rows where terms were manually chosen
    manually_selected.drop(["mapping_tool_response"], axis = 1, inplace = True)

    manually_selected["manually_selected_CURIE"] = manually_selected["manually_selected_CURIE"].apply(lambda x: return_curie_dict(x)) # convert | delimited strings to CURIE dict
    manually_selected["score"] = 1000  # human curated score = 1000
    manually_selected.rename(columns = {'manually_selected_CURIE':'mapping_tool_response'}, inplace = True)
    manually_selected = manually_selected[["mapping_tool", "term_type", "clintrial_term", "input_term", "mapping_tool_response", "score"]] # reorder columns to be same as the cache files we're appending to 
    manually_selected.to_csv("mapping_cache.tsv", mode='a', header=False, sep ="\t", index=False)


def check_against_cache(df_dict):

    print("Finding new terms to map, compare to cache")
    
    conditions_list = df_dict['conditions'].name.unique().tolist()
    conditions_list = [str(i) for i in conditions_list]
    conditions_list = list(set([i.lower() for i in conditions_list]))
    
    interventions_list = df_dict['interventions'].name.unique().tolist()
    interventions_list = [str(i) for i in interventions_list]
    interventions_list = list(set([i.lower() for i in interventions_list]))
    
    interventions_alts_list = df_dict['interventions_alts'].name.unique().tolist()
    interventions_alts_list = [str(i) for i in interventions_alts_list]
    interventions_alts_list = list(set([i.lower() for i in interventions_alts_list]))
    
    try:
        cache_manually_selected_terms()
    except:
        print("No manually selected terms file found")
    
    try:        
        cache_df = pd.read_csv("mapping_cache.tsv", sep ="\t", usecols = ['term_type', 'clintrial_term'], index_col=False, header=0, on_bad_lines = 'skip', encoding="utf-8")
        
        conditions_cache = cache_df[cache_df["term_type"] == "condition"]
        conditions_cache = conditions_cache['clintrial_term'].unique().tolist()
        conditions_cache = list(set([i.lower() for i in conditions_cache]))

        
        conditions_new = [x for x in conditions_list if x not in conditions_cache] # find conditions not in the cache (i.g. new conditions to map)
        conditions_new = list(filter(None, conditions_new))
        conditions_new = [str(i) for i in conditions_new]
        
        interventions_cache = cache_df[cache_df["term_type"] == "intervention"]
        interventions_cache = interventions_cache['clintrial_term'].unique().tolist()
        interventions_cache = list(set([i.lower() for i in interventions_cache]))
        
        interventions_new = [x for x in interventions_list if x not in interventions_cache] # find interventions not in the cache (i.g. new interventions to map)
        interventions_new = list(filter(None, interventions_new))
        interventions_new = [str(i) for i in interventions_new]
        
        interventions_alts_cache = cache_df[cache_df["term_type"] == "intervention_alternate"]
        interventions_alts_cache = interventions_alts_cache['clintrial_term'].unique().tolist()
        interventions_alts_cache = list(set([i.lower() for i in interventions_alts_cache]))
        
        interventions_alts_new = [x for x in interventions_alts_list if x not in interventions_alts_cache] # find interventions_alts not in the cache (i.g. new interventions_alts to map)
        interventions_alts_new = list(filter(None, interventions_alts_new))
        interventions_alts_new = [str(i) for i in interventions_alts_new]
        
    except:
        print("No cache of terms found. Proceeding to map entire KG from scratch")
        conditions_new = conditions_list
        interventions_new = interventions_list
        interventions_alts_new = interventions_alts_list
        
    dict_new_terms = {"conditions": conditions_new, "interventions": interventions_new, "interventions_alts": interventions_alts_new}

    return dict_new_terms


def get_nr_response(orig_term):
    def create_session():
        s = requests.Session()
        return s
 
    sess = create_session()
 
    """   Runs Name Resolver   """
    nr_url = 'https://name-resolution-sri.renci.org/lookup'
    max_retries = 3 
    
    input_term = orig_term # in MetaMap, we have to potentially deascii the term and lower case it...for Name Resolver, we don't need to do that. To keep columns consist with MetaMap output, we just keep it and say the original term and the input term are the same. For MetaMap, they might be different
    retries = 0
    params = {'string':orig_term, 'limit':1} # limit -1 makes this return all available equivalent CURIEs name resolver can give (deprecated)
    while retries <= max_retries:
        try:
            r = sess.post(nr_url, params=params)
            check_limit() # counts how many requests have been sent to NR. If limit of 40 have been sent, sleeps for 1 min
            if r.status_code == 200:
                mapping_tool_response = r.json()  # process Name Resolver response
                return mapping_tool_response
            else:
                return None
        except (requests.RequestException, ConnectionResetError, OSError) as ex:
            print(f"\nName Resolver request failed for term: {term}. Error: {ex}")
            retries += 1
            if retries < max_retries:
                print(f"Retrying ({retries}/{max_retries}) after a delay.")
                time.sleep(2 ** retries)  # Increase the delay between retries exponentially
            else:
                print(f"Max retries (Name Resolver) reached for term: {term}.")
                return None

# I'm only getting 1 concept from Name Resolver. 
# Both MetaMap and Name Resolver return several, 
# but I only take 1 from Name Resolver bc they have a preferred concept.
# MetaMap's 2nd or 3rd result is often the best one, so I collect all of them and try to score"

def process_metamap_concept(concept):
    concept = concept._asdict()
    concept_dict  = {"mapped_name": concept.get("preferred_name"),
                     "mapped_curie": concept.get("cui"),
                     "mapped_score": concept.get("score"),
                     "mapped_semtypes": concept.get("semtypes")}
    if not concept.get("preferred_name"): # if condition triggered if the concept dict looks like following, where AA is for Abbreviation.... {'index': 'USER', 'aa': 'AA', 'short_form': 'copd', 'long_form': 'chronic obstructive pulmonary disease', 'num_tokens_short_form': '1', 'num_chars_short_form': '4', 'num_tokens_long_form': '7', 'num_chars_long_form': '37', 'pos_info': '43:4'}
        concept_dict = None
    return concept_dict

def process_nameresolver_response(nr_response):              
    nr_curie = nr_response[0]["curie"]
    nr_name = nr_response[0]["label"]
    nr_type = nr_response[0]["types"][0]
    nr_score = nr_response[0]["score"]
    concept_dict = {"mapped_name": nr_name,
                    "mapped_curie": nr_curie,
                    "mapped_score": nr_score,
                    "mapped_semtypes": nr_type}
    return concept_dict


def run_mappers(term_pair, params, term_type, mapping_filename):
    # check_count()
    # while not terminate_flag.is_set():

        output = open(mapping_filename, 'a', newline='', encoding="utf-8") 
        csv_writer = csv.writer(output, delimiter='\t')

        orig_term = term_pair[0]
        input_term = term_pair[1]
        from_mapper = []
        mm = MetaMap.get_instance(metamap_dirs["metamap_base_dir"] + metamap_dirs["metamap_bin_dir"])

        # Format of output TSV: header = ['mapping_tool', 'term_type', 'clintrial_term', 'input_term', 'mapping_tool_response', 'score']

        if params.get("exclude_sts") is None: # exclude_sts is used for Interventions. restrict_to_sts is used for Conditions. So, the logic is, if we're mapping Conditions, execute "if" part of code. If we're mapping Interventions, execute "else" part of code
            try:
                concepts,error = mm.extract_concepts([input_term],
                                                     restrict_to_sts = params["restrict_to_sts"],
                                                     term_processing = params["term_processing"],
                                                     ignore_word_order = params["ignore_word_order"],
                                                     strict_model = params["strict_model"],)
                                                        
                if concepts:   # if MetaMap gives response, process response
                    mapping_tool = "metamap"
                    for concept in concepts:
                        concept_info = []
                        new_concept_dict = process_metamap_concept(concept)
                        concept_info.extend([mapping_tool, term_type, orig_term, input_term, new_concept_dict]) # score column is empty, Format of output TSV: header = ['mapping_tool', 'term_type', 'clintrial_term', 'input_term', 'mapping_tool_response', 'score']
                        from_mapper.append(concept_info)
                else:   # if MetaMap fails, try using Name Resolver and process response
                    nr_response = get_nr_response(orig_term)
                    if nr_response: # if Name Resolver gives response, process repsonse
                        input_term = orig_term # no preprocessing (lowercasing or deascii-ing) necessary to submit terms to Name Resolver (unlike MetaMap)
                        mapping_tool = "nameresolver"
                        concept_info = []
                        new_concept_dict = process_nameresolver_response(nr_response)
                        concept_info.extend([mapping_tool, term_type, orig_term, input_term, new_concept_dict]) # Add None for score column, empty bc not scored yet
                        from_mapper.append(concept_info)
                    else:
                        concept_info = []
                        # print("Nothing returned from NR or Metamap")
                        concept_info.extend(["mapping_tools_failed", term_type, orig_term, input_term, "mapping_tools_failed"])
                        from_mapper.append(concept_info)
            except:
                concept_info = []
                # print("Nothing returned from NR or Metamap")
                concept_info.extend(["mapping_tools_failed", term_type, orig_term, input_term, "mapping_tools_failed"])
                from_mapper.append(concept_info)
                
        else:   # Else block triggered if mapping Interventions
            try:
                concepts,error = mm.extract_concepts([input_term],
                                                     exclude_sts = params["exclude_sts"],
                                                     term_processing = params["term_processing"],
                                                     ignore_word_order = params["ignore_word_order"],
                                                     strict_model = params["strict_model"],) 
                                                       
                if concepts:   # if MetaMap gives response, process response
                    mapping_tool = "metamap"
                    for concept in concepts:
                        concept_info = []
                        new_concept_dict = process_metamap_concept(concept)
                        concept_info.extend([mapping_tool, term_type, orig_term, input_term, new_concept_dict]) # score column is empty, Format of output TSV: header = ['mapping_tool', 'term_type', 'clintrial_term', 'input_term', 'mapping_tool_response', 'score']
                        from_mapper.append(concept_info)
                else:   # if MetaMap fails, try using Name Resolver and process response
                    nr_response = get_nr_response(orig_term) 
                    if nr_response: # if Name Resolver gives response, process repsonse
                        input_term = orig_term # no preprocessing (lowercasing or deascii-ing) necessary to submit terms to Name Resolver (unlike MetaMap)
                        mapping_tool = "nameresolver"
                        concept_info = []
                        new_concept_dict = process_nameresolver_response(nr_response)
                        concept_info.extend([mapping_tool, term_type, orig_term, input_term, new_concept_dict])
                        from_mapper.append(concept_info)
                    else:
                        concept_info = []
                        # print("Nothing returned from NR or Metamap")
                        concept_info.extend(["mapping_tools_failed", term_type, orig_term, input_term, "mapping_tools_failed"])
                        from_mapper.append(concept_info)
            except:
                concept_info = []
                # print("Nothing returned from NR or Metamap")
                concept_info.extend(["mapping_tools_failed", term_type, orig_term, input_term, "mapping_tools_failed"])
                from_mapper.append(concept_info)
          
        for result in from_mapper:
            # print(result)
            if result[0] == "mapping_tools_failed":
                result.append(-1)
            else:
                result.append("unscored")
            # print(result)
        with csv_writer_lock:
            csv_writer.writerow(result)
    

def parallelize_mappers(term_pair_list, params, term_type, mapping_filename):
    # n_workers = 2 * multiprocessing.cpu_count() - 1
    n_workers = 6

    Parallel(n_jobs=n_workers,backend="multiprocessing")(
        delayed(run_mappers)
        (term_pair, params, term_type, mapping_filename) 
  for term_pair in term_pair_list
  )
          


def term_list_to_mappers(dict_new_terms):   
    metamap_version = [int(s) for s in re.findall(r'\d+', metamap_dirs.get('metamap_bin_dir'))] # get MetaMap version being run 
    deasciier = np.vectorize(de_ascii_er) # vectorize function
    

    # open mapping cache to add mapped terms
    mapping_filename = "mapping_cache.tsv"
    if os.path.exists(mapping_filename):
        output = open(mapping_filename, 'a', newline='', encoding="utf-8") 
        output.close()
    else:
        output = open(mapping_filename, 'w+', newline='', encoding='utf-8')
        col_names = ['mapping_tool', 'term_type', 'clintrial_term', 'input_term', 'mapping_tool_response', 'score']
        csv_writer = csv.writer(output, delimiter='\t')
        csv_writer.writerow(col_names)
        output.close()

    #  - Conditions
    condition_semantic_type_restriction = ['acab,anab,cgab,comd,dsyn,inpo,mobd,neop,patf,clna,fndg']  # see https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/SemanticTypes_2018AB.txt for semantic types ("acab,anab,etc.")
    conditions = dict_new_terms.get("conditions")
    condition_params = {"restrict_to_sts":condition_semantic_type_restriction, "term_processing":True, "ignore_word_order":True, "strict_model":False, "prune":30} # strict_model and relaxed_model are presumably opposites? relaxed_model = True is what I want, but that option appears to be broken in Pymetamap (returns no results when used). Using strict_model = False instead...

    #  - Interventions
    condition_semantic_type_restriction = ['acab,anab,cgab,comd,dsyn,inpo,mobd,neop,patf,clna,fndg']  # see https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/SemanticTypes_2018AB.txt for semantic types ("acab,anab,etc.")
    interventions = dict_new_terms.get("interventions")
    intervention_params = {"exclude_sts":condition_semantic_type_restriction, "term_processing":True, "ignore_word_order":True, "strict_model":False, "prune":30} # strict_model and relaxed_model are presumably opposites? relaxed_model = True is what I want, but that option appears to be broken in Pymetamap (returns no results when used). Using strict_model = False instead...we are also excluding all semantic types of condition bc interventions can be anything and moreover, prune=30 for memory issue (https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/OutOfMemory.pdf)

    #  - Alternate Interventions
    condition_semantic_type_restriction = ['acab,anab,cgab,comd,dsyn,inpo,mobd,neop,patf,clna,fndg']  # see https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/SemanticTypes_2018AB.txt for semantic types ("acab,anab,etc.")
    interventions_alts = dict_new_terms.get("interventions_alts")
    intervention_alts_params = intervention_params # same params as interventions
    
    chunksize = 10
    
    if metamap_version[0] >= 20:
        
        cons_processed = list(zip(conditions, conditions))  # these are lists of the same term repeated twice, bc MetaMap 2020 does not require deasciing, so the 2nd term remains unchanged and is a repeat of the first term
        ints_processed = list(zip(interventions, interventions))
        ints_alts_processed = list(zip(interventions_alts, interventions_alts))
    
        conditions_chunked = [cons_processed[i:i + chunksize] for i in range(0, len(cons_processed), chunksize)]  
        interventions_chunked = [ints_processed[i:i + chunksize] for i in range(0, len(ints_processed), chunksize)]  
        interventions_alts_chunked = [ints_alts_processed[i:i + chunksize] for i in range(0, len(ints_alts_processed), chunksize)] 
        
        print("MetaMap version >= 2020, conduct mapping on original terms")
        
        start_metamap_servers(metamap_dirs) # start the MetaMap servers

        LENGTH = len(cons_processed)  # Number of iterations required to fill progress bar (pbar)
        # pbar = tqdm(total=LENGTH, desc="% conditions mapped", position=0, leave=True, mininterval = LENGTH/20, bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}')  # Init progress bar
        pbar = tqdm(total=LENGTH, desc="% conditions mapped", position=0, leave=True, mininterval = LENGTH/40, bar_format='{l_bar}{bar:40}{r_bar}{bar:-10b}')  # Init progress bar
        for chunk in conditions_chunked:
            parallelize_mappers(chunk, condition_params, "condition", mapping_filename)
            pbar.update(n=len(chunk))


        LENGTH = len(ints_processed)  # Number of iterations required to fill progress bar (pbar)
        pbar = tqdm(total=LENGTH, desc="% interventions mapped", position=0, leave=True, mininterval = LENGTH/40, bar_format='{l_bar}{bar:40}{r_bar}{bar:-10b}')  # Init progress bar
        for chunk in interventions_chunked:
            parallelize_mappers(chunk, intervention_params, "intervention", mapping_filename)
            pbar.update(n=len(chunk))


        LENGTH = len(ints_alts_processed)  # Number of iterations required to fill progress bar (pbar)
        pbar = tqdm(total=LENGTH, desc="% alternate interventions mapped", position=0, leave=True, mininterval = LENGTH/40, bar_format='{l_bar}{bar:40}{r_bar}{bar:-10b}')  # Init progress bar
        for chunk in interventions_alts_chunked:
            parallelize_mappers(chunk, intervention_alts_params, "alternate_intervention", mapping_filename)
            pbar.update(n=len(chunk))

        stop_metamap_servers(metamap_dirs) # stop the MetaMap servers

        
    else:
        print("MetaMap version < 2020, conduct mapping on terms after removing ascii characters")
        
        deascii_cons = deasciier(conditions)
        deascii_ints = deasciier(interventions)
        deascii_int_alts = deasciier(interventions_alts)
                
        cons_processed = list(zip(conditions, deascii_cons)) # these are lists of the original term, and the deasciied term, bc MetaMap 2018 does not process ascii characters
        ints_processed = list(zip(interventions, deascii_ints))
        ints_alts_processed = list(zip(interventions_alts, deascii_int_alts))
        
        conditions_chunked = [cons_processed[i:i + chunksize] for i in range(0, len(cons_processed), chunksize)]  
        interventions_chunked = [ints_processed[i:i + chunksize] for i in range(0, len(ints_processed), chunksize)]  
        interventions_alts_chunked = [ints_alts_processed[i:i + chunksize] for i in range(0, len(ints_alts_processed), chunksize)] 
        
        start_metamap_servers(metamap_dirs) # start the MetaMap servers

        LENGTH = len(cons_processed)  # Number of iterations required to fill progress bar (pbar)
        pbar = tqdm(total=LENGTH, desc="% conditions mapped", position=0, leave=True, mininterval = LENGTH/40, bar_format='{l_bar}{bar:40}{r_bar}{bar:-10b}')  # Init progress bar
        for chunk in conditions_chunked:
            # parallelize_mappers(chunk, condition_params, "condition", csv_writer)
            parallelize_mappers(chunk, condition_params, "condition", mapping_filename)
            pbar.update(n=len(chunk))


        LENGTH = len(ints_processed)  # Number of iterations required to fill progress bar (pbar)
        pbar = tqdm(total=LENGTH, desc="% interventions mapped", position=0, leave=True, mininterval = LENGTH/40, bar_format='{l_bar}{bar:40}{r_bar}{bar:-10b}')  # Init progress bar
        for chunk in interventions_chunked:
            parallelize_mappers(chunk, intervention_params, "intervention", mapping_filename)

            pbar.update(n=len(chunk))


        LENGTH = len(ints_alts_processed)  # Number of iterations required to fill progress bar (pbar)
        pbar = tqdm(total=LENGTH, desc="% alternate interventions mapped", position=0, leave=True, mininterval = LENGTH/40, bar_format='{l_bar}{bar:40}{r_bar}{bar:-10b}')  # Init progress bar
        for chunk in interventions_alts_chunked:
            # parallelize_mappers(chunk, intervention_alts_params, "alternate_intervention", csv_writer)
            parallelize_mappers(chunk, intervention_alts_params, "alternate_intervention", mapping_filename)

            pbar.update(n=len(chunk))

        stop_metamap_servers(metamap_dirs) # stop the MetaMap servers
    
    # """ Remove duplicate rows """
    mapping_filename = "mapping_cache.tsv"
    cache = pd.read_csv(mapping_filename, sep='\t', index_col=False, header=0, encoding_errors='ignore', on_bad_lines='skip')
    cache = cache.drop_duplicates()
    cache.to_csv(mapping_filename, sep="\t", index=False, header=True) # output deduplicated cache terms to TSV
    

def score_mappings():
    print("Scoring cache")

    def get_max_score(str1, str2, old_score):
        
        try:
            if old_score == "unscored":
                sortratio_score = get_token_sort_ratio(str1, str2)
                similarity_score = get_similarity_score(str1, str2)
                max_score = max(sortratio_score, similarity_score)
                score = max_score
            else:
                score = old_score   
        except:
            score = old_score
        return score

    def wrap(x): # use this to convert string objects to dicts 
        try:
            a = ast.literal_eval(x)
            return(a)
        except:
            pass

    with pd.read_csv("mapping_cache.tsv", sep='\t', index_col=False, header=0, on_bad_lines = 'warn', usecols=lambda c: not c.startswith('Unnamed:'), chunksize=1000) as reader:
        write_header = True
        for chunk in reader:
            chunk["mapping_tool_response"] = chunk["mapping_tool_response"].apply(lambda x: wrap(x))
            mapping_info = chunk["mapping_tool_response"].apply(pd.Series, dtype='object')
            chunk["mapped_name"] = mapping_info["mapped_name"]
            chunk["score"] = chunk.apply(lambda x: get_max_score(x['input_term'], x['mapped_name'], x['score']), axis=1) # get score for score rows that are empty/not scored yet
            chunk.drop(["mapped_name"], axis = 1, inplace = True)
            chunk.to_csv(f'mapping_cache_scored_temp.tsv', sep="\t", index=False, header=write_header, mode = 'a', encoding="utf-8") # output to TSV
            write_header = False

    os.rename('mapping_cache.tsv','mapping_cache_backup.tsv')       
    os.rename('mapping_cache_scored_temp.tsv','mapping_cache.tsv')


def output_terms_files():
    print("Generating output files")

    """   Get high scorers   """
    cache = pd.read_csv("mapping_cache.tsv", sep='\t', index_col=False, header=0)
    cache['score'] = pd.to_numeric(cache['score'], errors='coerce')
    highscorers = cache[cache['score'] >= 80] 
    # test = highscorers.groupby('clintrial_term')
    idx = highscorers.groupby('clintrial_term')['score'].idxmax()  # group by the clinical trial term and get the highest scoring
    auto_selected = highscorers.loc[idx]
    auto_selected.to_csv(f'autoselected_terms.tsv', sep="\t", index=False, header=True) # output to TSV

    """   Get low scorers, aggregate for manual selections  """
    low_scorers = cache[cache['score'] < 80]
    manual_review = low_scorers[~low_scorers.clintrial_term.isin(highscorers['clintrial_term'].unique().tolist())] # there are terms autoselected that have mappings that didn't pass threshold too, but we want to consider that term mapped. So get rid of these rows too
    mapping_tool_response = manual_review['mapping_tool_response'].apply(lambda x: wrap(x))
    manual_review = manual_review.copy()
    mapping_tool_response = mapping_tool_response.apply(pd.Series, dtype='object')
    manual_review.loc[:, 'mapping_tool_response_lists'] = mapping_tool_response.values.tolist()
    manual_review.drop('mapping_tool_response', axis=1, inplace=True)
    manual_review = manual_review[["mapping_tool", "term_type", "clintrial_term", "mapping_tool_response_lists", "input_term", "score"]]
    # manual_review['mapping_tool_response_lists'] = manual_review['mapping_tool_response_lists'].apply(lambda x: ' | '.join(x) if isinstance(x, list) else None)  # Multindexing does not work on lists, so remove the CURIE information out of the list to enable this
    manual_review['mapping_tool_response'] = [' | '.join(map(str, l)) for l in manual_review['mapping_tool_response_lists']]
    manual_review.drop('mapping_tool_response_lists', axis=1, inplace=True)
    manual_review = manual_review.sort_values(by=["mapping_tool", "term_type", "clintrial_term", "input_term"], ascending=False)
    manual_review.set_index(["mapping_tool", "term_type", "clintrial_term", "input_term"], inplace=True)   # create index
    manual_review['manually_selected_CURIE'] = None # make a column 
    manual_review.to_excel('manual_review.xlsx', engine='xlsxwriter', index=True)

    sys.stdout.flush() 

    print("Done\n")


if __name__ == "__main__":
    # flag_and_path = get_raw_ct_data() # download raw data

    flag_and_path = {"term_program_flag": False, "data_extracted_path": "/15TB_2/gglusman/datasets/clinicaltrials/AACT-20240227", "date_string": "02_27_2024"}
    # flag_and_path = {"term_program_flag": False, "data_extracted_path": "/Users/Kamileh/Work/ISB/NCATS_BiomedicalTranslator/Projects/ClinicalTrials/ETL_Python/data/02_27_2024_extracted", "date_string": "02_27_2024"}
    global metamap_dirs
    metamap_dirs = check_os()
    subset_size = None
    df_dict = read_raw_ct_data(flag_and_path, subset_size) # read the clinical trial data
    dict_new_terms = check_against_cache(df_dict) # use the existing cache of MetaMapped terms so that only new terms are mapped
    term_list_to_mappers(dict_new_terms)
    score_mappings()
    output_terms_files()

    