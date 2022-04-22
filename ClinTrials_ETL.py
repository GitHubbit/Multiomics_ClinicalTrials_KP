#!/usr/bin/env python3

import os
import pandas as pd
# shell command: brew install postgresql 
import psycopg2
import pandas.io.sql as sqlio


# connect to DB and get the column names of the table
con = None
con = psycopg2.connect(database="aact")
con.rollback()
cursor = con.cursor()

con.autocommit = True # SQL statement is treated as a transaction and is automatically committed right after it is executed
sql = '''SELECT * FROM ctgov.conditions;'''
cursor.execute(sql)
column_names = [desc[0] for desc in cursor.description]
tuples = cursor.fetchall()
con.close()
df = pd.DataFrame(tuples, columns=column_names)

print(df.head(5))