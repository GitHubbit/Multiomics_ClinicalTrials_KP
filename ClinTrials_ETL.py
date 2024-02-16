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
sys.path.insert(0, '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/pymetamap-master')
from pymetamap import MetaMap  # https://github.com/AnthonyMRios/pymetamap/blob/master/pymetamap/SubprocessBackend.py
from pandas import ExcelWriter
import ast
import glob
from tqdm import tqdm
import subprocess
import shlex
from collections import Counter

# %pip install thefuzz
# %pip install levenshtein
# %pip install xlsxwriter

from thefuzz import fuzz # fuzzy matching explained: https://www.datacamp.com/tutorial/fuzzy-string-python


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

def start_metamap_servers(metamap_dirs):
    global metamap_pos_server_dir
    global metamap_wsd_server_dir
    metamap_pos_server_dir = 'bin/skrmedpostctl' # Part of speech tagger
    metamap_wsd_server_dir = 'bin/wsdserverctl' # Word sense disambiguation 
    
    metamap_executable_path_pos = os.path.join(metamap_dirs['metamap_base_dir'], metamap_pos_server_dir)
    metamap_executable_path_wsd = os.path.join(metamap_dirs['metamap_base_dir'], metamap_wsd_server_dir)
    command_pos = [metamap_executable_path_pos, 'start']
    command_wsd = [metamap_executable_path_wsd, 'start']

    # Start servers, with open portion redirects output of metamap server printing output to NULL
    with open(os.devnull, "w") as fnull:
        result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
        result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
    sleep(5)

def stop_metamap_servers(metamap_dirs):
    metamap_executable_path_pos = os.path.join(metamap_dirs['metamap_base_dir'], metamap_pos_server_dir)
    metamap_executable_path_wsd = os.path.join(metamap_dirs['metamap_base_dir'], metamap_wsd_server_dir)
    command_pos = [metamap_executable_path_pos, 'stop']
    command_wsd = [metamap_executable_path_wsd, 'stop']
    
    # Stop servers, with open portion redirects output of metamap server printing output to NULL
    with open(os.devnull, "w") as fnull:
        result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
        result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
    sleep(2)  
    
def add_mappings_to_cache(flag_and_path):
    relevant_date = flag_and_path["date_string"]   # get date of bulk download of clinical trial data
    with open("metamapped_terms_cache.tsv", 'a+', encoding="utf-8") as cache:
        with open(f"{relevant_date}_metamap_output.tsv", 'r', encoding="utf-8", errors='ignore') as new_metamapped_terms:
            # Read the first line from new_metamapped_terms to move the cursor
            line = new_metamapped_terms.readline()

            # Move the cursor to the position after the first line
            while line:
                line = new_metamapped_terms.readline()
                if line:
                    # Append the line to file_1
                    cache.write(line)
    """ Remove duplicate rows from cache """
    cache = pd.read_csv("metamapped_terms_cache.tsv", sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    cache = cache.drop_duplicates()
    cache.to_csv('metamapped_terms_cache.tsv', sep="\t", index=False, header=True) # output deduplicated cache terms to TSV

def add_manually_selected_terms_to_cache():
    # -----     ------     GENERATE MANUALLY SELECTED CACHE     -----     ------  #
    try:
        #  --- --- --   CONDITIONS     --- --- --   #
        files = glob.glob("*.xlsx")
        conditions_manselected_files = [i for i in files if "conditions_manual_review" in i if not i.startswith("~")][0]  
        conditions_manselected = pd.read_excel(conditions_manselected_files)
        conditions_manselected.name.ffill(inplace=True)
        conditions_manselected.orig_con.ffill(inplace=True)
        conditions_manselected = conditions_manselected[~conditions_manselected['manually_selected_CURIE'].isnull()] # check if the conditions got mapped to any CURIEs
        conditions_manselected.drop(["curie_info"], axis = 1, inplace = True)
        conditions_manselected.rename(columns = {'name':'original_clin_trial_term', 'orig_con':'modified_clin_trial_term'}, inplace = True)

        with open('conditions_manually_selected_cache.tsv', 'a') as output:
            conditions_manselected.to_csv(output, mode='a',sep="\t", index=False, header=output.tell()==0)
        """ Remove duplicate rows from cache """
        cache = pd.read_csv("conditions_manually_selected_cache.tsv", sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
        cache = cache.drop_duplicates()
        cache.to_csv('conditions_manually_selected_cache.tsv', sep="\t", index=False, header=True) # output deduplicated cache terms to TSV

        #  --- --- --   INTERVENTIONS and Alternate INTERVENTIONS   --- --- --   #
        files = glob.glob("*.xlsx")
        interventions_manselected_files = [i for i in files if "interventions_manual_review" in i if not i.startswith("~")][0]  
        interventions_manselected = pd.read_excel(interventions_manselected_files)
        interventions_manselected.name.ffill(inplace=True)
        interventions_manselected.orig_int.ffill(inplace=True)
        interventions_manselected = interventions_manselected[~interventions_manselected['manually_selected_CURIE'].isnull()] # check if the conditions got mapped to any CURIEs
        interventions_manselected.drop(["curie_info", "description"], axis = 1, inplace = True)
        interventions_manselected.rename(columns = {'name':'original_clin_trial_term', 'orig_int':'modified_clin_trial_term'}, inplace = True)

        with open('interventions_manually_selected_cache.tsv', 'a') as output:
            interventions_manselected.to_csv(output, mode='a',sep="\t", index=False, header=output.tell()==0)
        """ Remove duplicate rows from cache """
        cache = pd.read_csv("interventions_manually_selected_cache.tsv", sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
        cache = cache.drop_duplicates()
        cache.to_csv('interventions_manually_selected_cache.tsv', sep="\t", index=False, header=True) # output deduplicated cache terms to TSV
    except:
        print("No terms in manual select column; either column is empty or bug. Proceeding without them")
        
def check_os():
    if "linux" in sys.platform:
        print("Linux platform detected")
        metamap_base_dir = "{}/metamap/".format(pathlib.Path.cwd().parents[0])
        metamap_bin_dir = 'bin/metamap20'
    else:
        metamap_base_dir = '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/public_mm/' # for running on local
        metamap_bin_dir = 'bin/metamap18'
        
    return {"metamap_base_dir":metamap_base_dir, "metamap_bin_dir":metamap_bin_dir}  

def get_raw_ct_data():
    term_program_flag = True
    global data_dir
    global data_extracted

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
            print("\nFailed to scrape AACT for download. Please navigate to https://aact.ctti-clinicaltrials.org/download and manually download zip file.")
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
        conditions_df = pd.read_csv(data_extracted + '/conditions.txt', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
        interventions_df = pd.read_csv(data_extracted + '/interventions.txt', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
        interventions_alts_df = pd.read_csv(data_extracted + '/intervention_other_names.txt', sep='|', index_col=False, header=0, on_bad_lines = 'warn')

        if subset_size:   # if a subset size is given, we are running this script on a small subset of the dataset
            conditions_df = conditions_df.sample(n=subset_size)
            interventions_df = interventions_df.sample(n=subset_size)
            interventions_alts_df = interventions_alts_df.sample(n=subset_size)

    return {"conditions": conditions_df, "interventions": interventions_df, "interventions_alts": interventions_alts_df}


def run_metamap(term_pair, params, mm, cond_or_inter, csv_writer):
    
    orig_term = term_pair[0]
    input_term = term_pair[1]
    from_metamap = []
    
    if params.get("exclude_sts") is None: # exclude_sts is used for Interventions. restrict_to_sts is used for Conditions. So, the logic is, if we're mapping Conditions, execute "if" part of code. If we're mapping Interventions, execute "else" part of code
        try:
            concepts,error = mm.extract_concepts([input_term],
                                                 restrict_to_sts = params["restrict_to_sts"],
                                                 term_processing = params["term_processing"],
                                                 ignore_word_order = params["ignore_word_order"],
                                                 strict_model = params["strict_model"],
                                                )
            for concept in concepts:
                concept_info = []
                concept = concept._asdict()
                concept_info.extend([cond_or_inter, orig_term, input_term])
                concept_info.extend([concept.get(k) for k in ['preferred_name', 'cui', 'score', 'semtypes']])
                from_metamap.append(concept_info)
        except:
            from_metamap.extend([cond_or_inter, orig_term, input_term, None, None, None, None])
    else:
        try:
            concepts,error = mm.extract_concepts([input_term],
                                                 exclude_sts = params["exclude_sts"],
                                                 term_processing = params["term_processing"],
                                                 ignore_word_order = params["ignore_word_order"],
                                                 strict_model = params["strict_model"],
                                                )
            for concept in concepts:
                concept_info = []
                concept = concept._asdict()
                concept_info.extend([cond_or_inter, orig_term, input_term])
                concept_info.extend([concept.get(k) for k in ['preferred_name', 'cui', 'score', 'semtypes']])
                from_metamap.append(concept_info)
        except:
            from_metamap.extend([cond_or_inter, orig_term, input_term, None, None, None, None])
        
    for result in from_metamap:
        # print(result)
        csv_writer.writerow(result)
    return from_metamap

def parallelize_metamap(term_pair_list, params, cond_or_inter, flag_and_path, csv_writer):
    LENGTH = len(term_pair_list)  # Number of iterations required to fill progress bar (pbar)
    pbar = tqdm(total=LENGTH, desc="% {}s mapped".format(cond_or_inter), position=0, leave=True, mininterval = LENGTH/20, bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}')  # Init progress bar

    start_metamap_servers(metamap_dirs) # start the MetaMap servers
    mm = MetaMap.get_instance(metamap_dirs["metamap_base_dir"] + metamap_dirs["metamap_bin_dir"])
    with concurrent.futures.ThreadPoolExecutor((multiprocessing.cpu_count()*2) - 1) as executor:
        futures = [executor.submit(run_metamap, term_pair, params, mm, cond_or_inter, csv_writer) for term_pair in term_pair_list]
        for _ in concurrent.futures.as_completed(futures):
            pbar.update(n=1)  # Increments counter
    stop_metamap_servers(metamap_dirs) # stop the MetaMap servers
    

def term_list_to_cache(df_dict, flag_and_path):
    print("Sending terms to cache; picking terms that have not been mapped previously")

    #  --- --- --     CONDITIONS     --- --- --   #

    # retrieve conditions from the new data dump and from the cache, compare, get the diff so as to map just newly encountered terms (conditions)
    conditions_list = df_dict['conditions'].name.unique().tolist()
    conditions_list = [str(i) for i in conditions_list]
    conditions_list = set([i.lower() for i in conditions_list])

    # Get CONDITIONS already run through MetaMap (1st cache)
    mm_cache_file = "metamapped_terms_cache.tsv"
    mm_cache_df = pd.read_csv(mm_cache_file, sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    mm_conditions_cache = mm_cache_df[mm_cache_df["term_type"] == "condition"]
    mm_conditions_cache = mm_conditions_cache['original_clin_trial_term'].unique().tolist()
    mm_conditions_cache = list(set([i.lower() for i in mm_conditions_cache]))

    # Get CONDITIONS already manually mapped (2nd cache)
    manual_cache_file = "conditions_manually_selected_cache.tsv"
    manual_cache_df = pd.read_csv(manual_cache_file, sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    manual_conditions_cache = manual_cache_df['original_clin_trial_term'].unique().tolist()
    manual_conditions_cache = list(set([i.lower() for i in manual_conditions_cache]))

    # merge the 2 caches
    conditions_cache = mm_conditions_cache + manual_conditions_cache

    conditions_new = [x for x in conditions_list if x not in conditions_cache] # find conditions not in the cache (i.g. new conditions to map)
    conditions_new = list(filter(None, conditions_new))
    conditions_new = [str(i) for i in conditions_new]

    #  --- --- --     INTERVENTIONS     --- --- --   #

    # retrieve conditions from the new data dump and from the cache, compare, get the diff so as to map just newly encountered terms (interventions)
    interventions_list = df_dict['interventions'].name.unique().tolist()
    interventions_list = [str(i) for i in interventions_list]
    interventions_list = set([i.lower() for i in interventions_list])

    # Get INTERVENTIONS already run through MetaMap (1st cache)
    mm_cache_file = "metamapped_terms_cache.tsv"
    mm_cache_df = pd.read_csv(mm_cache_file, sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    mm_interventions_cache = mm_cache_df[mm_cache_df["term_type"] == "intervention"]
    mm_interventions_cache = mm_interventions_cache['original_clin_trial_term'].unique().tolist()
    mm_interventions_cache = list(set([i.lower() for i in mm_interventions_cache]))

    # Get INTERVENTIONS already manually mapped (2nd cache)
    manual_cache_file = "interventions_manually_selected_cache.tsv"
    manual_cache_df = pd.read_csv(manual_cache_file, sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    manual_interventions_cache = manual_cache_df['original_clin_trial_term'].unique().tolist()
    manual_interventions_cache = list(set([i.lower() for i in manual_interventions_cache]))

    # merge the 2 caches
    interventions_cache = mm_interventions_cache + manual_interventions_cache

    interventions_new = [x for x in interventions_list if x not in interventions_cache] # find conditions not in the cache (i.g. new conditions to map)
    interventions_new = list(filter(None, interventions_new))
    interventions_new = [str(i) for i in interventions_new]

    #  --- --- --     ALTERNATE INTERVENTIONS     --- --- --   #

    # retrieve conditions from the new data dump and from the cache, compare, get the diff so as to map just newly encountered terms (interventions)
    interventions_alts_list = df_dict['interventions_alts'].name.unique().tolist()
    interventions_alts_list = [str(i) for i in interventions_alts_list]
    interventions_alts_list = set([i.lower() for i in interventions_alts_list])

    # Get ALTERNATE INTERVENTIONS already run through MetaMap (1st cache)
    mm_cache_file = "metamapped_terms_cache.tsv"
    mm_cache_df = pd.read_csv(mm_cache_file, sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    mm_interventions_alts_cache = mm_cache_df[mm_cache_df["term_type"] == "alternate_intervention"]
    mm_interventions_alts_cache = mm_interventions_alts_cache['original_clin_trial_term'].unique().tolist()
    mm_interventions_alts_cache = list(set([i.lower() for i in mm_interventions_alts_cache]))

    # Get ALTERNATE INTERVENTIONS already manually mapped (2nd cache)
    manual_cache_file = "interventions_manually_selected_cache.tsv"
    manual_cache_df = pd.read_csv(manual_cache_file, sep='\t', index_col=False, header=0, on_bad_lines = 'warn')
    manual_interventions_alts_cache = manual_cache_df['original_clin_trial_term'].unique().tolist()
    manual_interventions_alts_cache = list(set([i.lower() for i in manual_interventions_alts_cache]))

    # merge the 2 caches
    interventions_alts_cache = mm_interventions_alts_cache + manual_interventions_alts_cache

    interventions_alts_new = [x for x in interventions_alts_list if x not in interventions_alts_cache] # find conditions not in the cache (i.g. new conditions to map)
    interventions_alts_new = list(filter(None, interventions_alts_new))
    interventions_alts_new = [str(i) for i in interventions_alts_new]

    dict_new_terms = {"conditions": conditions_new, "interventions": interventions_new, "interventions_alts": interventions_alts_new}
    return dict_new_terms


def term_list_to_mm(dict_new_terms, flag_and_path):   
    metamap_version = [int(s) for s in re.findall(r'\d+', metamap_dirs.get('metamap_bin_dir'))] # get MetaMap version being run 
    relevant_date = flag_and_path["date_string"]   # get date of bulk download of clinical trial data
    deasciier = np.vectorize(de_ascii_er) # vectorize function

    # prep output file of Metamap results
    filename = f"{relevant_date}_metamap_output.tsv"
    metamap_output = open(filename, 'w+', newline='')
    col_names = ['term_type', 'original_clin_trial_term', 'modified_clin_trial_term', 'metamap_preferred_name', 'metamap_cui', 'metamap_score', 'metamap_semantic_type']
    csv_writer = csv.writer(metamap_output, delimiter='\t')
    csv_writer.writerow(col_names)

    condition_semantic_type_restriction = ['acab,anab,cgab,comd,dsyn,inpo,mobd,neop,patf,clna,fndg']  # see https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/SemanticTypes_2018AB.txt for semantic types ("acab,anab,etc.")

    orig_cons = dict_new_terms.get("conditions")
    condition_params = {"restrict_to_sts":condition_semantic_type_restriction, "term_processing":True, "ignore_word_order":True, "strict_model":False} # strict_model and relaxed_model are presumably opposites? relaxed_model = True is what I want, but that option appears to be broken in Pymetamap (returns no results when used). Using strict_model = False instead...

    orig_ints = dict_new_terms.get("interventions")
    intervention_params = {"exclude_sts":condition_semantic_type_restriction, "term_processing":True, "ignore_word_order":True, "strict_model":False} # strict_model and relaxed_model are presumably opposites? relaxed_model = True is what I want, but that option appears to be broken in Pymetamap (returns no results when used). Using strict_model = False instead...

    orig_int_alts = dict_new_terms.get("interventions_alts")
    intervention_alts_params = intervention_params

    if metamap_version[0] >= 20:
        print("MetaMap version >= 2020, conduct mapping on original terms")
        # parallelize_metamap(orig_cons, condition_params, "condition", flag_and_path, csv_writer)
        parallelize_metamap(list(zip(orig_cons, orig_cons)), condition_params, "condition", flag_and_path, csv_writer)
        parallelize_metamap(list(zip(orig_ints, orig_ints)), intervention_params, "intervention", flag_and_path, csv_writer)
        parallelize_metamap(list(zip(orig_int_alts, orig_int_alts)), intervention_alts_params, "alternate_intervention", flag_and_path, csv_writer)
    else:
        print("MetaMap version < 2020, conduct mapping on terms after removing ascii characters")
        deascii_cons = deasciier(orig_cons)
        deascii_ints = deasciier(orig_ints)
        deascii_int_alts = deasciier(orig_int_alts)
        parallelize_metamap(list(zip(orig_cons, deascii_cons)), condition_params, "condition", flag_and_path, csv_writer)
        parallelize_metamap(list(zip(orig_ints, deascii_ints)), intervention_params, "intervention", flag_and_path, csv_writer)
        parallelize_metamap(list(zip(orig_int_alts, deascii_int_alts)), intervention_alts_params, "alternate_intervention", flag_and_path, csv_writer)

    metamap_output.close()
    add_mappings_to_cache(flag_and_path)
    add_manually_selected_terms_to_cache()


def map_to_trial(flag_and_path):
    print("\nMapping UMLS CURIEs and names back to clinical trials")
    relevant_date = flag_and_path["date_string"]   # get date of bulk download of clinical trial data
    metamap_version = [int(s) for s in re.findall(r'\d+', metamap_dirs.get('metamap_bin_dir'))] # get MetaMap version being run 

    metamap_input = "metamapped_terms_cache.tsv"
    metamapped = pd.read_csv(metamap_input, sep='\t', index_col=False, header=0, encoding="utf-8", on_bad_lines = 'warn')

    metamap_semantic_types = pd.read_csv("MetaMap_SemanticTypes_2018AB.txt", on_bad_lines = 'warn') # get the full names of the semantic types so we know what we're looking at
    metamapped['metamap_semantic_type'] = metamapped['metamap_semantic_type'].str.replace(r'\[|\]', '', regex=True)
    sem_type_col_names = ["abbv", "group", "semantic_type_full"]
    metamap_semantic_types = pd.read_csv("MetaMap_SemanticTypes_2018AB.txt", sep="|", index_col=False, header=None, names=sem_type_col_names, on_bad_lines = 'warn')
    sem_type_dict = dict(zip(metamap_semantic_types['abbv'], metamap_semantic_types['semantic_type_full'])) # make a dict of semantic type abbv and full name

    metamapped['metamap_semantic_type'] = metamapped['metamap_semantic_type'].apply(lambda x: x.split(',') if isinstance(x, str) else np.nan) # Handle NaN (None) values in metamap_semantic_type column
    metamapped['metamap_semantic_type'] = metamapped['metamap_semantic_type'].apply(lambda x: '|'.join([sem_type_dict[term] if term in sem_type_dict else term for term in x]) if isinstance(x, list) else x) # map semantic type abbreviations to the full name of the semantic type

    metamapped['metamap_preferred_name'] = metamapped['metamap_preferred_name'].str.lower()
    metamapped = metamapped.dropna(axis=0)
    metamapped = metamapped[["term_type", "original_clin_trial_term", "modified_clin_trial_term", "metamap_cui","metamap_preferred_name", "metamap_semantic_type"]]

    metamapped["metamap_term_info"] = metamapped[["metamap_cui", "metamap_preferred_name", "metamap_semantic_type"]].values.tolist() 
    metamapped.drop(["metamap_cui", "metamap_preferred_name", "metamap_semantic_type"], axis = 1, inplace = True)
    metamapped = metamapped.groupby(['term_type', 'original_clin_trial_term'])['metamap_term_info'].agg(list).reset_index()

    data_extracted = flag_and_path["data_extracted_path"]
    # read in pipe-delimited files 
    conditions_df = pd.read_csv(data_extracted + '/conditions.txt', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
    conditions_mapped = conditions_df.copy()

    metamapped_con = metamapped.loc[metamapped['term_type'] == "condition"]
    mapper_con = dict(zip(metamapped_con['original_clin_trial_term'], metamapped_con['metamap_term_info'])) # make a dict from the metamapped cache to map conditions
    conditions_mapped['curie_info'] = conditions_mapped['downcase_name'].map(mapper_con)

    data_extracted = flag_and_path["data_extracted_path"]
    # read in pipe-delimited files 
    interventions_df = pd.read_csv(data_extracted + '/interventions.txt', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
    interventions_mapped = interventions_df.copy()
    interventions_mapped["downcase_name"] = interventions_mapped['name'].str.lower()

    metamapped_int = metamapped.loc[(metamapped['term_type'] == "intervention") | (metamapped['term_type'] == "alternate_intervention")]
    mapper_int = dict(zip(metamapped_int['original_clin_trial_term'], metamapped_int['metamap_term_info'])) # make a dict from the metamapped cache to map interventions
    interventions_mapped['curie_info'] = interventions_mapped['downcase_name'].map(mapper_int)

    interventions_alts_df = pd.read_csv(data_extracted + '/intervention_other_names.txt', sep='|', index_col=False, header=0, on_bad_lines = 'warn')
    interventions_alts_mapped = interventions_alts_df.copy()
    interventions_alts_mapped.drop(["id"], axis=1, inplace=True)
    interventions_alts_mapped = interventions_alts_mapped.merge(interventions_df[["id", "intervention_type", "description"]], left_on='intervention_id', right_on='id', how='left') 
    interventions_alts_mapped.drop(["id"], axis=1, inplace=True)
    interventions_alts_mapped["downcase_name"] = interventions_alts_mapped['name'].str.lower()

    metamapped_int = metamapped.loc[(metamapped['term_type'] == "intervention") | (metamapped['term_type'] == "alternate_intervention")]
    mapper_int = dict(zip(metamapped_int['original_clin_trial_term'], metamapped_int['metamap_term_info'])) # make a dict from the metamapped cache to map interventions
    interventions_alts_mapped['curie_info'] = interventions_alts_mapped['downcase_name'].map(mapper_int)
    interventions_alts_mapped = interventions_alts_mapped[["intervention_id", "nct_id", "intervention_type", "name", "description", "downcase_name", "curie_info"]]

    # # # conditions_actually_mapped = conditions_mapped[~conditions_mapped['curie_info'].isnull()] # check if the conditions got mapped to any CURIEs
    # # # interventions_actually_mapped = interventions_mapped[~interventions_mapped['curie_info'].isnull()] # check if the interventions got mapped to any CURIEs
    # # # interventions_alts_actually_mapped = interventions_alts_mapped[~interventions_alts_mapped['curie_info'].isnull()] # check if the interventions got mapped to any CURIEs

    conditions_mapped.to_csv('{}_conditions_mapped.tsv'.format(relevant_date), sep="\t", index=False, header=True) # output conditions to TSV
    interventions_mapped.to_csv('{}_interventions_mapped.tsv'.format(relevant_date), sep="\t", index=False, header=True) # output interventions to TSV
    interventions_alts_mapped.to_csv('{}_interventions_alternates_mapped.tsv'.format(relevant_date), sep="\t", index=False, header=True) # output alternate interventions to TSV

def score_mappings(flag_and_path):

    relevant_date = flag_and_path["date_string"]   # get date of bulk download of clinical trial data
    pattern_outside = r'(?<=\().+?(?=\))|([^(]+)'
    pattern_inside = r'\(([^)]+)\)'

    # -----     ------     CONDITIONS     -----     ------  #

    print("Scoring mappings for CONDITIONS")
    
    # prep output file of condition scored results
    filename = f"{relevant_date}_conditions_scored.tsv"
    conditions_scored_output = open(filename, 'w+', newline='')
    col_names = ["condition_id", "nct_id", "name", "orig_con", "curie_info", "orig_con_outside", "orig_con_inside"]
    csv_writer = csv.writer(conditions_scored_output, delimiter='\t')
    csv_writer.writerow(col_names)

    # with pd.read_csv("{}_conditions_mapped.tsv".format(relevant_date), sep='\t', index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000, nrows=4000) as reader:
    with pd.read_csv("{}_conditions_mapped.tsv".format(relevant_date), sep='\t', index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000) as reader:

        for conditions_df_chunk in reader:
            conditions = conditions_df_chunk.copy()
            conditions.rename(columns = {'downcase_name':'orig_con','id': 'condition_id'}, inplace = True)

            matches_outside = conditions['orig_con'].str.extract(pattern_outside)
            conditions['orig_con_outside'] = matches_outside[0].fillna(np.nan).replace([np.nan], [None]) # if no matches inside ()/to regex pattern, then fill with None
            matches_inside = conditions['orig_con'].str.extract(pattern_inside)
            conditions['orig_con_inside'] = matches_inside[0].fillna(np.nan).replace([np.nan], [None]) # # if no matches outside ()/to regex pattern, then fill with None

            conditions = conditions.fillna(np.nan).replace([np.nan], [None]) # replace NaN values in dataframe with None (just for consistency)
            conditions_cols = conditions.columns.tolist() # get all column names as list so we can dynamically change column values as we iterate through rows 

            cols_to_check = [ele for ele in conditions.columns if any([substr in ele for substr in ['_con']])] # find columns with "con" in it 

            for i, row in conditions.iterrows():
                orig_row = row
                new_row = orig_row.mask(orig_row.duplicated(), None) # if term/condition inside or outside parentheses is duplicated, replace the duplicate term with None so we don't waste time scoring it
                conditions.loc[i, conditions_cols] = new_row 
                curies_sublists_scored = []
                for col_name in cols_to_check: # check only columns with terms to score
                    value = row[col_name]
                    curie_info = row["curie_info"]
                    if None not in [value, curie_info]: # if the column has a term to score in that row...
                        curie_sublists = ast.literal_eval(curie_info)
                        for sublist in curie_sublists:
                            sublist.append(f'sort_ratio: {get_token_sort_ratio(value, sublist[1])}') # get the sort ratio score for that term and the CURIE from MetaMap
                            sublist.append(f'similarity_score: {get_similarity_score(value, sublist[1])}') # get the similarity score for that term and the CURIE from MetaMap
                            curies_sublists_scored.append(sublist)
                if curies_sublists_scored:
                    curies_sublists_scored = [list(y) for y in set([tuple(x) for x in curies_sublists_scored])] # remove any duplicate sublists
                    conditions.at[i, "curie_info"] = curies_sublists_scored
                else:
                    conditions.at[i, "curie_info"] = None
                scored_row = conditions.loc[i, :].values.tolist()
                scored_row = list(map(lambda x: str(x) if x is not None else "", scored_row))
                csv_writer.writerow(scored_row)                
    conditions_scored_output.close() 

    # # -----     ------     INTERVENTIONS     -----     ------  #

    print("Scoring mappings for INTERVENTIONS")

    # prep output file of intervention scored results
    filename = f"{relevant_date}_interventions_scored.tsv"
    interventions_scored_output = open(filename, 'w+', newline='')
    col_names = ["intervention_id", "nct_id", "intervention_type", "name", "description", "orig_int", "curie_info", "orig_int_inside", "orig_int_outside"] 
    csv_writer = csv.writer(interventions_scored_output, delimiter='\t')
    csv_writer.writerow(col_names)

    # with pd.read_csv("{}_interventions_mapped.tsv".format(relevant_date), sep='\t', index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000, nrows=4000) as reader:
    with pd.read_csv("{}_interventions_mapped.tsv".format(relevant_date), sep='\t', index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000) as reader:

        for interventions_df_chunk in reader:
            interventions = interventions_df_chunk.copy()
            interventions.rename(columns = {'downcase_name':'orig_int','id': 'intervention_id'}, inplace = True)

            matches_outside = interventions['orig_int'].str.extract(pattern_outside)
            interventions['orig_int_outside'] = matches_outside[0].fillna(np.nan).replace([np.nan], [None]) # if no matches inside ()/to regex pattern, then fill with None
            matches_inside = interventions['orig_int'].str.extract(pattern_inside)
            interventions['orig_int_inside'] = matches_inside[0].fillna(np.nan).replace([np.nan], [None]) # # if no matches outside ()/to regex pattern, then fill with None

            interventions = interventions.fillna(np.nan).replace([np.nan], [None]) # replace NaN values in dataframe with None (just for consistency)
            interventions_cols = interventions.columns.tolist() # get all column names as list so we can dynamically change column values as we iterate through rows 

            cols_to_check = [ele for ele in interventions.columns if any([substr in ele for substr in ['_int']])] # find columns with "int" in it 

            for i, row in interventions.iterrows():
                orig_row = row
                new_row = orig_row.mask(orig_row.duplicated(), None) # if term/condition inside or outside parentheses is duplicated, replace the duplicate term with None so we don't waste time scoring it
                interventions.loc[i, interventions_cols] = new_row 
                curies_sublists_scored = []
                for col_name in cols_to_check: # check only columns with terms to score
                    value = row[col_name]
                    curie_info = row["curie_info"]
                    if None not in [value, curie_info]: # if the column has a term to score in that row...
                        curie_sublists = ast.literal_eval(curie_info)
                        for sublist in curie_sublists:
                            sublist.append(f'sort_ratio: {get_token_sort_ratio(value, sublist[1])}') # get the sort ratio score for that term and the CURIE from MetaMap
                            sublist.append(f'similarity_score: {get_similarity_score(value, sublist[1])}') # get the similarity score for that term and the CURIE from MetaMap
                            curies_sublists_scored.append(sublist)
                if curies_sublists_scored:
                    curies_sublists_scored = [list(y) for y in set([tuple(x) for x in curies_sublists_scored])] # remove any duplicate sublists
                    interventions.at[i, "curie_info"] = curies_sublists_scored
                else:
                    interventions.at[i, "curie_info"] = None
                scored_row = interventions.loc[i, :].values.tolist()
                scored_row = list(map(lambda x: str(x) if x is not None else "", scored_row)) 
                csv_writer.writerow(scored_row)
    interventions_scored_output.close()
    
    # # -----     ------     ALTERNATE INTERVENTIONS     -----     ------  #
    
    print("Scoring mappings for ALTERNATE INTERVENTIONS")

    # prep output file of intervention scored results
    filename = f"{relevant_date}_interventions_alternates_scored.tsv"
    interventions_alts_scored_output = open(filename, 'w+', newline='')
    col_names = ["intervention_id", "nct_id", "intervention_type", "name", "description", "orig_int_alt", "curie_info", "orig_int_alt_inside", "orig_int_alt_outside"] 

    csv_writer = csv.writer(interventions_alts_scored_output, delimiter='\t')
    csv_writer.writerow(col_names)

    # with pd.read_csv("{}_interventions_alternates_mapped.tsv".format(relevant_date), sep='\t', index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000, nrows=4000) as reader:
    with pd.read_csv("{}_interventions_alternates_mapped.tsv".format(relevant_date), sep='\t', index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000) as reader:

        for interventions_df_chunk in reader:
            interventions = interventions_df_chunk.copy()
            interventions.rename(columns = {'downcase_name':'orig_int_alt'}, inplace = True)

            matches_outside = interventions['orig_int_alt'].str.extract(pattern_outside)
            interventions['orig_int_alt_outside'] = matches_outside[0].fillna(np.nan).replace([np.nan], [None]) # if no matches inside ()/to regex pattern, then fill with None
            matches_inside = interventions['orig_int_alt'].str.extract(pattern_inside)
            interventions['orig_int_alt_inside'] = matches_inside[0].fillna(np.nan).replace([np.nan], [None]) # # if no matches outside ()/to regex pattern, then fill with None

            interventions = interventions.fillna(np.nan).replace([np.nan], [None]) # replace NaN values in dataframe with None (just for consistency)
            interventions_cols = interventions.columns.tolist() # get all column names as list so we can dynamically change column values as we iterate through rows 

            cols_to_check = [ele for ele in interventions.columns if any([substr in ele for substr in ['_int']])] # find columns with "int" in it 

            for i, row in interventions.iterrows():
                orig_row = row
                new_row = orig_row.mask(orig_row.duplicated(), None) # if term/condition inside or outside parentheses is duplicated, replace the duplicate term with None so we don't waste time scoring it
                interventions.loc[i, interventions_cols] = new_row 
                curies_sublists_scored = []
                for col_name in cols_to_check: # check only columns with terms to score
                    value = row[col_name]
                    curie_info = row["curie_info"]
                    if None not in [value, curie_info]: # if the column has a term to score in that row...
                        curie_sublists = ast.literal_eval(curie_info)
                        for sublist in curie_sublists:
                            sublist.append(f'sort_ratio: {get_token_sort_ratio(value, sublist[1])}') # get the sort ratio score for that term and the CURIE from MetaMap
                            sublist.append(f'similarity_score: {get_similarity_score(value, sublist[1])}') # get the similarity score for that term and the CURIE from MetaMap
                            curies_sublists_scored.append(sublist)
                if curies_sublists_scored:
                    curies_sublists_scored = [list(y) for y in set([tuple(x) for x in curies_sublists_scored])] # remove any duplicate sublists
                    interventions.at[i, "curie_info"] = curies_sublists_scored
                else:
                    interventions.at[i, "curie_info"] = None
                scored_row = interventions.loc[i, :].values.tolist()
                scored_row = list(map(lambda x: str(x) if x is not None else "", scored_row)) 
                csv_writer.writerow(scored_row)
    interventions_scored_output.close()   

def auto_select_curies(flag_and_path):
    print("Auto-selecting high scoring CURIEs")
    relevant_date = flag_and_path["date_string"]   # get date of bulk download of clinical trial data

    # select sublist with highest scoring term
    def filter_and_select_sublist(sublists):  # function to find the highest score of a CURIE, and pick that curie if it's greater than threshold of 88
        try:
            selected_sublist = None
            if pd.isnull(sublists):
                return selected_sublist
            else:
                high_score = -1
                sublists = ast.literal_eval(sublists)
                for sublist in sublists:
                    if len(sublist) >= 4:
                        sort_ratio = int(sublist[3].split(": ")[1])
                        sim_score = int(sublist[4].split(": ")[1])
                        max_score = max(sort_ratio, sim_score)
                        if max_score > high_score:
                            high_score  = max_score
                            selected_sublist = sublist
                        else:
                            continue
                        if max_score < 80:
                            selected_sublist = None
            return selected_sublist
        except:
            return None

    # convert string to list, return empty list
    def convert_to_list(x):
        try:
            listx = ast.literal_eval(x)
            return listx
        except Exception as e:
            # print(e)
            return []

    # # -----     ------     CONDITIONS     -----     ------  #

    # if the previous autoselected or dump of "no CURIE selected" files exist, delete them
    if os.path.exists(f'{relevant_date}_conditions_autoselected.tsv'):
        os.remove(f'{relevant_date}_conditions_autoselected.tsv')
    if os.path.exists(f'{relevant_date}_conditions_manual_review.tsv'):
        os.remove(f'{relevant_date}_conditions_manual_review.tsv')

    # with pd.read_csv(f"{relevant_date}_conditions_scored.tsv", sep='\t', usecols=["condition_id", "nct_id", "name", "orig_con", "curie_info"], index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000, nrows=4000) as reader:
    with pd.read_csv(f"{relevant_date}_conditions_scored.tsv", sep='\t', usecols=["condition_id", "nct_id", "name", "orig_con", "curie_info"], index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000) as reader:

        print("Auto-selecting CONDITIONS CURIEs")
        write_header = True

        for chunk_df in reader:
            """  Create an output TSV of CURIEs that are auto-selected based on passing the threshold of scoring > 80  """
            chunk_df['auto_selected_curie'] = chunk_df['curie_info'].apply(filter_and_select_sublist)  # select CURIE that scores highest using filter_and_select_sublist function = auto-select
            auto_selected = chunk_df.loc[chunk_df['auto_selected_curie'].notnull(),]  # get the rows where a CURIE has been auto-selected
            auto_selected.to_csv(f'{relevant_date}_conditions_autoselected.tsv', sep="\t", index=False, header=write_header, mode = 'a') # output to TSV

            manual_review = chunk_df.loc[~chunk_df['auto_selected_curie'].notnull(),] # get the rows where a CURIE has not been selected or MetaMap did not return a term
            manual_review = manual_review[["name", "orig_con", "curie_info", "auto_selected_curie"]]
            manual_review = manual_review.drop_duplicates() # run this on the chunked df just bc...we repeat this on the entire dataframe later
            manual_review.to_csv(f'{relevant_date}_conditions_manual_review.tsv', sep="\t", index=False, header=write_header, mode = 'a') # output to TSV

            write_header = False  

    """  Create an output TSV of possible CURIEs available for each term that was not auto-selected  """
    manual_review = pd.read_csv(f'{relevant_date}_conditions_manual_review.tsv', sep="\t", on_bad_lines = 'warn')
    manual_review.drop(["auto_selected_curie"], axis = 1, inplace = True)   # drop the autoselect column bc it's empty, these rows are specifically ones where nothing was selected
    manual_review = manual_review.drop_duplicates()
    manual_review.loc[:, "curie_info"] = manual_review.curie_info.apply(lambda x: convert_to_list(x)) # in order to multi-index, we have to group-by the original input term. To do this, first convert the curie_info column to list of lists
    manual_review = manual_review.explode('curie_info')  # explode that column so every sublist is on a separate row
    manual_review['curie_info'] = manual_review['curie_info'].apply(lambda x: x[:3] if isinstance(x, list) else None)   # remove the scores (sort_ratio and similarity score) from the list, don't need them and they compromise readability of manual outputs 
    manual_review['curie_info'] = manual_review['curie_info'].apply(lambda x: '|--|'.join(x) if isinstance(x, list) else None)  # Multindexing does not work on lists, so remove the CURIE information out of the list to enable this

    manual_review['temp'] = "temp"   # create a temp column to facilitate multi-indexing
    manual_review.set_index(["name", "orig_con", "curie_info"], inplace=True)   # create index
    manual_review.drop(["temp"], axis = 1, inplace = True)   # drop the temp column
    manual_review['manually_selected_CURIE'] = None # make a column 

    manual_review.to_excel('{}_conditions_manual_review.xlsx'.format(relevant_date), engine='xlsxwriter', index=True)

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
    #     display(manual_review[:2000])        

    # # # -----     ------     INTERVENTIONS     -----     ------  #

    # if the previous autoselected or dump of "no CURIE selected" files exist, delete them
    if os.path.exists(f'{relevant_date}_interventions_autoselected.tsv'):
        os.remove(f'{relevant_date}_interventions_autoselected.tsv')
    if os.path.exists(f'{relevant_date}_interventions_manual_review.tsv'):
        os.remove(f'{relevant_date}_interventions_manual_review.tsv')

    # with pd.read_csv(f"{relevant_date}_interventions_scored.tsv", sep='\t', usecols=["intervention_id", "description", "nct_id", "name", "orig_int", "curie_info"], index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000, nrows=4000) as reader:
    with pd.read_csv(f"{relevant_date}_interventions_scored.tsv", sep='\t', usecols=["intervention_id", "description", "nct_id", "name", "orig_int", "curie_info"], index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000) as reader:

        print("Auto-selecting INTERVENTIONS CURIEs")
        write_header = True

        for chunk_df in reader:
            """  Create an output TSV of CURIEs that are auto-selected based on passing the threshold of scoring > 80  """
            chunk_df['auto_selected_curie'] = chunk_df['curie_info'].apply(filter_and_select_sublist)  # select CURIE that scores highest using filter_and_select_sublist function = auto-select
            auto_selected = chunk_df.loc[chunk_df['auto_selected_curie'].notnull(), ]  # get the rows where a CURIE has been auto-selected
            auto_selected.to_csv(f'{relevant_date}_interventions_autoselected.tsv', sep="\t", index=False, header=write_header, mode = 'a') # output to TSV

            manual_review = chunk_df.loc[~chunk_df['auto_selected_curie'].notnull(),] # get the rows where a CURIE has not been selected or MetaMap did not return a term
            manual_review = manual_review[["name", "orig_int", "curie_info", "auto_selected_curie", "description"]]
            manual_review = manual_review.drop_duplicates(["orig_int", "description"]) # run this on the chunked df just bc...we repeat this on the entire dataframe later
            manual_review.to_csv(f'{relevant_date}_interventions_manual_review.tsv', sep="\t", index=False, header=write_header, mode = 'a') # output to TSV

            write_header = False  

    """  Create an output TSV of possible CURIEs available for each term that was not auto-selected  """
    manual_review = pd.read_csv(f'{relevant_date}_interventions_manual_review.tsv', sep="\t", on_bad_lines = 'warn')
    manual_review.drop(["auto_selected_curie"], axis = 1, inplace = True)   # drop the autoselect column bc it's empty, these rows are specifically ones where nothing was selected
    manual_review = manual_review.drop_duplicates()
    manual_review.loc[:, "curie_info"] = manual_review.curie_info.apply(lambda x: convert_to_list(x)) # in order to multi-index, we have to group-by the original input term. To do this, first convert the curie_info column to list of lists
    manual_review = manual_review.explode('curie_info')  # explode that column so every sublist is on a separate row
    manual_review['curie_info'] = manual_review['curie_info'].apply(lambda x: x[:3] if isinstance(x, list) else None)   # remove the scores (sort_ratio and similarity score) from the list, don't need them and they compromise readability of manual outputs 
    manual_review['curie_info'] = manual_review['curie_info'].apply(lambda x: '|--|'.join(x) if isinstance(x, list) else None)  # Multindexing does not work on lists, so remove the CURIE information out of the list to enable this

    manual_review['temp'] = "temp"   # create a temp column to facilitate multi-indexing
    manual_review.set_index(["name", "description", "orig_int", "curie_info"], inplace=True)   # create index
    manual_review.drop(["temp"], axis = 1, inplace = True)   # drop the temp column
    manual_review['manually_selected_CURIE'] = None # make a column 

    manual_review.to_excel('{}_interventions_manual_review.xlsx'.format(relevant_date), engine='xlsxwriter', index=True)  


    # # # -----     ------     ALTERNATE INTERVENTIONS     -----     ------  #

    # if the previous autoselected or dump of "no CURIE selected" files exist, delete them
    if os.path.exists(f'{relevant_date}_interventions_alternates_autoselected.tsv'):
        os.remove(f'{relevant_date}_interventions_alternates_autoselected.tsv')
    if os.path.exists(f'{relevant_date}_interventions_alternates_manual_review.tsv'):
        os.remove(f'{relevant_date}_interventions_alternates_manual_review.tsv')

    # with pd.read_csv(f"{relevant_date}_interventions_alternates_scored.tsv", sep='\t', usecols=["intervention_id", "description", "nct_id", "name", "orig_int_alt", "curie_info"], index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000, nrows=4000) as reader:
    with pd.read_csv(f"{relevant_date}_interventions_alternates_scored.tsv", sep='\t', usecols=["intervention_id", "description", "nct_id", "name", "orig_int_alt", "curie_info"], index_col=False, header=0, on_bad_lines = 'warn', chunksize=1000) as reader:

        print("Auto-selecting ALTERNATE INTERVENTIONS CURIEs")
        write_header = True

        for chunk_df in reader:
            """  Create an output TSV of CURIEs that are auto-selected based on passing the threshold of scoring > 80  """
            chunk_df['auto_selected_curie'] = chunk_df['curie_info'].apply(filter_and_select_sublist)  # select CURIE that scores highest using filter_and_select_sublist function = auto-select
            auto_selected = chunk_df.loc[chunk_df['auto_selected_curie'].notnull(), ]  # get the rows where a CURIE has been auto-selected
            auto_selected.to_csv(f'{relevant_date}_interventions_alternates_autoselected.tsv', sep="\t", index=False, header=write_header, mode = 'a') # output to TSV

            manual_review = chunk_df.loc[~chunk_df['auto_selected_curie'].notnull(),] # get the rows where a CURIE has not been selected or MetaMap did not return a term
            manual_review = manual_review[["name", "orig_int_alt", "curie_info", "auto_selected_curie", "description"]]
            manual_review = manual_review.drop_duplicates(["orig_int_alt", "description"]) # run this on the chunked df just bc...we repeat this on the entire dataframe later
            manual_review.to_csv(f'{relevant_date}_interventions_alternates_manual_review.tsv', sep="\t", index=False, header=write_header, mode = 'a') # output to TSV

            write_header = False  

    """  Create an output TSV of possible CURIEs available for each term that was not auto-selected  """
    manual_review = pd.read_csv(f'{relevant_date}_interventions_alternates_manual_review.tsv', sep="\t", on_bad_lines = 'warn')
    manual_review.drop(["auto_selected_curie"], axis = 1, inplace = True)   # drop the autoselect column bc it's empty, these rows are specifically ones where nothing was selected
    manual_review = manual_review.drop_duplicates()
    manual_review.loc[:, "curie_info"] = manual_review.curie_info.apply(lambda x: convert_to_list(x)) # in order to multi-index, we have to group-by the original input term. To do this, first convert the curie_info column to list of lists
    manual_review = manual_review.explode('curie_info')  # explode that column so every sublist is on a separate row
    manual_review['curie_info'] = manual_review['curie_info'].apply(lambda x: x[:3] if isinstance(x, list) else None)   # remove the scores (sort_ratio and similarity score) from the list, don't need them and they compromise readability of manual outputs 
    manual_review['curie_info'] = manual_review['curie_info'].apply(lambda x: '|--|'.join(x) if isinstance(x, list) else None)  # Multindexing does not work on lists, so remove the CURIE information out of the list to enable this

    manual_review['temp'] = "temp"   # create a temp column to facilitate multi-indexing
    manual_review.set_index(["name", "description", "orig_int_alt", "curie_info"], inplace=True)   # create index
    manual_review.drop(["temp"], axis = 1, inplace = True)   # drop the temp column
    manual_review['manually_selected_CURIE'] = None # make a column 

    manual_review.to_excel('{}_interventions_alternates_manual_review.xlsx'.format(relevant_date), engine='xlsxwriter', index=True)    


def run_ETL(subset_size):

    start_time_begin = time.time()
    flag_and_path = get_raw_ct_data() # download raw data
    end_time_download = time.time()
    elapsed_time = end_time_download - start_time_begin
    hours, minutes, seconds = convert_seconds_to_hms(elapsed_time)
    print(f"\nApproximate runtime for downloading or locating raw data: {hours} hours, {minutes} minutes, {seconds} seconds")

    global metamap_dirs
    metamap_dirs = check_os()
    df_dict = read_raw_ct_data(flag_and_path, subset_size) # read the clinical trial data
    dict_new_terms = term_list_to_cache(df_dict, flag_and_path) # use the existing cache of MetaMapped terms so that only new terms are mapped

    start_time_mm = time.time()
    term_list_to_mm(dict_new_terms, flag_and_path) # map new terms using MetaMap
    end_time_mm = time.time()
    elapsed_time = end_time_mm - start_time_mm
    hours, minutes, seconds = convert_seconds_to_hms(elapsed_time)
    print(f"Approximate runtime for mapping: {hours} hours, {minutes} minutes, {seconds} seconds")

    map_to_trial(flag_and_path) # map MetaMap terms back to trial 
    score_mappings(flag_and_path) # score the mappings
    auto_select_curies(flag_and_path) # select CURIEs automatically that pass score threshold
    
    # compile_curies_for_trials(flag_and_path) # select CURIEs automatically that pass score threshold

    end_time_end = time.time()
    elapsed_time = end_time_end - start_time_begin
    hours, minutes, seconds = convert_seconds_to_hms(elapsed_time)
    print(f"Approximate runtime for overall mapping: {hours} hours, {minutes} minutes, {seconds} seconds")


def test_or_prod():
    global subset_size
    subset_size = 100
    print(f"The test run of this code performs the construction of the KG on a random subset of {subset_size} Conditions, {subset_size} Interventions, and {subset_size} Alternate Interventions from Clinical Trials.")
    test_or_prod = input("Is this a test run or the production of a new version of the KG? Enter Test for test, or Prod for production: ")
    if test_or_prod == "Test":
        run_ETL(subset_size)
    elif test_or_prod == "Prod":
        subset_size = None
        run_ETL(subset_size)
    else:
        print("Bad input")
        sys.exit(0)


if __name__ == "__main__":
    test_or_prod()
