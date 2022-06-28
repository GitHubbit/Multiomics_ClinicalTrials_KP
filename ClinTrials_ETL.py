#!/usr/bin/env python3

import os
import pandas as pd
# shell command: brew install postgresql 
import psycopg2
import pandas.io.sql as sqlio
from config import API_KEY
from time import sleep
import re

# import requests
# from bs4 import BeautifulSoup


#####	----	****	----	----	****	----	----	****	----    #####

# ACCESSING DATA BY DOWNLOADING STATIC COPY OF DATABASE, AND INSTALLING POSTGRESQL SOFTWARE TO THEN POPULATE DATABASE ON MACHINE
# THEN CONNECTING TO DB

# connect to DB and get the column names of the table

# extract data from ClinicalTrials.gov
def get_trial_data():
    con = None
    con = psycopg2.connect(database="aact")
    con.rollback()
    cursor = con.cursor()

    con.autocommit = True # SQL statement is treated as a transaction and is automatically committed right after it is executed
    # grab the conditions
    sql = '''SELECT * FROM ctgov.conditions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    conditions_df = pd.DataFrame(tuples, columns=column_names)

    #grab the browse_conditions
    sql = '''SELECT * FROM ctgov.browse_conditions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    browse_conditions_df = pd.DataFrame(tuples, columns=column_names)

    #grab the interventions
    sql = '''SELECT * FROM ctgov.interventions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    interventions_df = pd.DataFrame(tuples, columns=column_names)

    #grab the browse_interventions
    sql = '''SELECT * FROM ctgov.browse_interventions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    browse_interventions_df = pd.DataFrame(tuples, columns=column_names)

    con.close()


    # rename and drop df relevant columns to prepare for merging
    interventions_df = interventions_df.rename(columns={'id': 'int_id',
                                                        'nct_id': 'int_nctid',
                                                        'intervention_type': 'int_type',
                                                        'name': 'int_name',
                                                        'description': 'int_description'})
    interventions_df = interventions_df.drop(columns=['int_id', 'int_description'])
    conditions_df = conditions_df.rename(columns={'id': 'con_id',
                                                  'nct_id': 'con_nctid',
                                                  'name': 'con_name',
                                                  'downcase_name': 'con_downcase_name'})
    conditions_df = conditions_df.drop(columns=['con_id', 'con_name'])
    browse_interventions_df = browse_interventions_df.rename(columns={'id': 'browseint_id',
                                                                      'nct_id': 'browseint_nctid',
                                                                      'mesh_term': 'browseint_meshterm',
                                                                      'downcase_mesh_term': 'browseint_meshterm_downcase',
                                                                      'mesh_type': 'browseint_meshtype'})

    browse_interventions_df = browse_interventions_df.drop(columns=['browseint_id', 'browseint_meshterm'])
    browse_conditions_df = browse_conditions_df.rename(columns={'id': 'browsecon_id',
                                                                'nct_id': 'browsecon_nctid',
                                                                'mesh_term': 'browsecon_meshterm',
                                                                'downcase_mesh_term': 'browsecon_meshterm_downcase',
                                                                'mesh_type': 'browsecon_meshtype'})
    browse_conditions_df = browse_conditions_df.drop(columns=['browsecon_id', 'browsecon_meshterm'])                                                                                                                          

    # merge conditions_df and interventions_df since they have relevant terms 
    df = pd.merge(conditions_df, interventions_df, left_on='con_nctid', right_on = 'int_nctid')
    df_dedup = df.drop_duplicates(subset = ['con_downcase_name', 'int_name'],
                                          keep = 'first').reset_index(drop = True)

    return df_dedup

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     display(df_dedup.head(100))


    # with pd.option_context('expand_frame_repr', False, 'display.max_rows', None):
    # 	print(df_dedup[['con_nctid', 'con_downcase_name', 'int_name']].head(100))


    # df_dedup[['con_downcase_name', 'int_name']].to_csv('ClinTrials_KG_nodes.csv', sep ='\t', index=False)

# the boss wants to see the properly formatted .tsv file output first, even if it's filled with dummy data...

# Massage data extracted from ClinicalTrials.gov
def preprocess_ct_data(ct_data):    
    # first get only relevant columns from DB
    ct_extract = pd.DataFrame(ct_data[['con_nctid', 'con_downcase_name', 'int_type', 'int_name']])
    ct_extract = ct_extract.rename(columns={'con_nctid': 'nctid'})
    # get CURIE column for nct_id column (https://bioregistry.io/registry/clinicaltrials)

    ct_extract['nctid_curie'] = ct_extract['nctid']
    ct_extract['nctid_curie'] = 'clinicaltrials:' + ct_extract['nctid'].astype(str)  # generate CURIEs for each clinical trials study

    ct_final = pd.DataFrame(columns=['subject','predicate','object', 'subject_name','object_name','category'])

    ct_final['subject'] = 'condition:' + ct_extract['con_downcase_name'].astype(str)
    ct_final['predicate'] = 'biolink:hypothesized_for'
    ct_final['object'] = 'intervention:' + ct_extract['int_name'].astype(str)   # this will not all be RxNorm CURIEs since some interventions are not drugs
    ct_final.subject_name = ct_extract.con_downcase_name
    ct_final.object_name = ct_extract.int_name
    ct_final.category = 'biolink:Association'
    ct_final['nctid'] = ct_extract['nctid']
    ct_final['nctid_curie'] = ct_extract['nctid_curie']

    return(ct_final)


def select_cols_to_preprocess(ct_preprocessed):
    # prep a list to hold dicts of de-asciied cols from the ClinTrials df
    processed_ct_cols = []
    conditions_metamap_compatible = {}
    interventions_metamap_compatible = {}

    # get unique lists of the concepts to MetaMap from the ClinTrials df
    conditions = list(set(list(ct_preprocessed['subject_name'])))
    interventions = list(filter(None, ct_preprocessed['object_name'].unique()))

    if not all(condition.isascii() for condition in conditions):
        print("Non-ASCII chars are detected in Conditions col from ClinTrial (not checking other cols), proceed with text processing")
        print("This step is unnecessary for MetaMap 2020+")
        # process lists to remove non-ascii chars (this is not required for MetaMap 2020!!!!)
        conditions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in conditions)
        interventions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in interventions)
        
        conditions_metamap_compatible["conditions"] = conditions_translated
        interventions_metamap_compatible["interventions"] = interventions_translated

        processed_ct_cols.append(conditions_metamap_compatible)
        processed_ct_cols.append(interventions_metamap_compatible)
    else:
        print("No non-ASCII chars detected, no text processing required")
        conditions_metamap_compatible["conditions"] = conditions
        interventions_metamap_compatible["interventions"] = interventions

        processed_ct_cols.append(conditions_metamap_compatible)
        processed_ct_cols.append(interventions_metamap_compatible)

    return(processed_ct_cols)
    



    # while ascii_flag:
    #     for condition in conditions:
    #         if condition.isascii():





        # if not condition.isascii():
        #     print("Non-ASCII chars are detected in Conditions col from ClinTrial (not checking other cols), proceed with text processing")
        #     # process lists to remove non-ascii chars (this is not required for MetaMap 2020!!!!)
        #     conditions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in conditions)
        #     interventions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in interventions)
        #     processed_ct_cols.append(("conditions", conditions_translated))
        #     processed_ct_cols.append(("interventions", interventions_translated))
        #     break
        # else:
        #     print("All ASCII char, no text processing required")
        #     processed_ct_cols.append(("conditions", conditions))
        #     processed_ct_cols.append(("interventions", interventions))
            # break  

    # print(processed_ct_cols)




    # if any(0 <= ord(char) <= 127 for char in condition for condition in conditions): # note we're only checking the Conditions col
    #     print("Non-ASCII chars are detected in Conditions col from ClinTrial (not checking other cols), proceed with text processing")
    #     # process lists to remove non-ascii chars (this is not required for MetaMap 2020!!!!)
    #     conditions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in conditions)
    #     interventions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in interventions)
    #     processed_ct_cols.append(("conditions", conditions_translated))
    #     processed_ct_cols.append(("interventions", interventions_translated))
    # else:
    #     print("All ASCII char, no text processing required")
    #     processed_ct_cols.append(("conditions", conditions))
    #     processed_ct_cols.append(("interventions", interventions))  

    # print(processed_ct_cols)
    # return(processed_ct_cols)
    #     # # process lists to remove non-ascii chars (this is not required for MetaMap 2020!!!!)
    # conditions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in conditions)
    # interventions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in interventions)


    # conditions_translated = ("conditions", conditions_translated)
    # print(conditions_translated[1])
    # interventions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in interventions)


    # de_ascii_lists.append({"conditions", conditions_translated})
    # de_ascii_lists.append({"interventions", interventions_translated})



def preprocess_cols_for_metamap(text):
    non_ascii = "[^\x00-\x7F]"
    pattern = re.compile(r"[^\x00-\x7F]")
    non_ascii_text = re.sub(pattern, ' ', text)
    return non_ascii_text

# # setup MetaMap
# def prep_for_metamap():
#     # Setup UMLS Server
#     metamap_base_dir = '/gwshare/umls_2021/metamap/public_mm/'
#     metamap_bin_dir = 'bin/metamap20'
#     metamap_pos_server_dir = 'bin/skrmedpostctl'
#     metamap_wsd_server_dir = 'bin/wsdserverctl'

#     # Start servers
#     os.system(metamap_base_dir + metamap_pos_server_dir + ' start') # Part of speech tagger
#     os.system(metamap_base_dir + metamap_wsd_server_dir + ' start') # Word sense disambiguation 

#     # Sleep a bit to give time for these servers to start up
#     sleep(20)


# def run_metamap():
#     pass






def driver():
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     display(df_dedup.head(100))
    ct_data = get_trial_data()
    ct_preprocessed = preprocess_ct_data(ct_data)
    select_cols_to_preprocess(ct_preprocessed)

    # translated_lists = select_cols_to_preprocess(ct_preprocessed)
    # print(translated_lists)

    # preprocess_cols_for_metamap(ct_preprocessed)
    # prep_for_metamap(ct_preprocessed)

    # # print(ct_preprocessed.head(20))
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     print(ct_preprocessed.head(20))
    # # mapper(ct_preprocessed)


driver()



#####	----	****	----	----	****	----	----	****	----    #####


