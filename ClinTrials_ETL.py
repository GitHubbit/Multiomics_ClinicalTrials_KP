#!/usr/bin/env python3

import os
import pandas as pd
# shell command: brew install postgresql 
import psycopg2
import pandas.io.sql as sqlio

import requests
from bs4 import BeautifulSoup


#####	----	****	----	----	****	----	----	****	----    #####

# ACCESSING DATA BY DOWNLOADING STATIC COPY OF DATABASE, AND INSTALLING POSTGRESQL SOFTWARE TO THEN POPULATE DATABASE ON MACHINE
# THEN CONNECTING TO DB

connect to DB and get the column names of the table
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

# testing if dataframes successfully inhaled database data
# print(conditions_df.head(5))
# print(browse_conditions_df.head(5))
# print(interventions_df.head(5).to_string())

# prob best to stick to MESH terms since they're standardized, these are in the Browse_Conditions and Browse_Interventions tables
# print(browse_conditions_df.head(5).to_string())
# print(browse_interventions_df.head(5).to_string())

# see list of unique MESH conditions
unique_cond = browse_conditions_df.downcase_mesh_term.unique()
# print(*(sorted(unique_cond)), sep = "\n")

# see list of unique MESH interventions
unique_int = browse_interventions_df.downcase_mesh_term.unique()
# print(*(sorted(unique_int)), sep = "\n")

# need to link condition MESH terms to intervention MESH terms via NCT ID
# get dfs ready

browse_interventions_df = browse_interventions_df.drop(columns=['id', 'mesh_term', 'mesh_type'])
browse_conditions_df = browse_conditions_df.drop(columns=['id', 'mesh_term', 'mesh_type'])
# rename columns called "MESH term"...otherwise when merging later we will get 2 columns in merged table with identical name
browse_interventions_df = browse_interventions_df.rename(columns={'downcase_mesh_term': 'Intervention_MESH'})
browse_conditions_df = browse_conditions_df.rename(columns={'downcase_mesh_term': 'Condition_MESH'})


print(browse_conditions_df.head(5).to_string())
print(browse_interventions_df.head(5).to_string())


# performing right outer join, where the left table contains conditions, and the right interventions
# I only want to preserve conditions in the table for which an intervention exists

df = pd.merge(browse_conditions_df, browse_interventions_df, on='nct_id')

print("\n\n\n")
print(df.head(5).to_string())

#####	----	****	----	----	****	----	----	****	----    #####


# access data by downloading pipe-delimited flat files
# url = "https://aact.ctti-clinicaltrials.org/pipe_files"
# response = requests.get(url)
# soup = BeautifulSoup(response.text, "lxml")



# body = soup.find_all('td', attrs={'class': 'file-archive'}) #Find all
# for el in body:
# 	link = el.find('a')
# 	print(link)
# 	print("\n")

