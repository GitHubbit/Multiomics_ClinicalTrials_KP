#!/usr/bin/env conda run -n ct_extract_env python
import ClinTrials_ETL

# def get_release(version_date):
# 	return ClinTrials_ETL.version_date

# get_release()

def get_release(self):
    # hard-coded release for demo purpose
    # "self" is a dumper instance, see:
    # https://github.com/biothings/biothings.api/blob/master/biothings/hub/dataload/dumper.py
    return "1.0"