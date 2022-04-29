#!/usr/bin/env python3

import os
import pandas as pd
# shell command: brew install postgresql 
import psycopg2
import pandas.io.sql as sqlio

# import requests
# from bs4 import BeautifulSoup


#####	----	****	----	----	****	----	----	****	----    #####

# ACCESSING DATA BY DOWNLOADING STATIC COPY OF DATABASE, AND INSTALLING POSTGRESQL SOFTWARE TO THEN POPULATE DATABASE ON MACHINE
# THEN CONNECTING TO DB

# connect to DB and get the column names of the table
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

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#     display(df_dedup.head(100))


with pd.option_context('expand_frame_repr', False, 'display.max_rows', None):
	print(df_dedup[['con_nctid', 'con_downcase_name', 'int_name']].head(100))


df_dedup[['con_downcase_name', 'int_name']].to_csv('ClinTrials_KG_nodes.csv', sep ='\t', index=False)












#####	----	****	----	----	****	----	----	****	----    #####


