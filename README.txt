Note: had started to do this in R at first (folder: ClinTrial_ETL_R), but preference of 
others was Python (script is called ClinicalTrials_ETL.py)

Downloaded static database from AACT (Clinical Trials.gov, https://aact.ctti-clinicaltrials.org/download)
Unpack it using Postgresql software
store as zip in current folder
access it in pandas data tables in python (folder: ClinTrial_ETL_Python)

created dir called metamap (mkdir metamap) to contain all MetaMap related programs

Download MetaMap software mapping tool from to run locally:
https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/run-locally/MainDownload.html
	- downloaded public_mm_linux_main_2020.tar.bz2 locally
	- then copied it from local to hypatia server via scp (scp /Downloads/)
	- run bunzip -c public_mm_linux_main_2020.tar.bz2 | tar xvf - ...
	- Follow successive installation instructions from NIH NLM:
		https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/documentation/Installation.html
You must have Java 1.8 or JDK
MetaMap requires 2 servers (SKR Part of Speech Tagger and WSD Word Sense Disambiguation) for our purposes:

We start and stop those servers by calling it from the Python script

Git clone wrapper to run MetaMap in Python script (ClinicalTrials_ETL.py)
git clone https://github.com/AnthonyMRios/pymetamap

MetaMap 2016 requires removal of non-ASCII char. MetaMap 2020 does not. Original script was developed on Mac/OSX, for which there's no MetaMap 2020, only MetaMap 2016 (at time of development, June 2022)
I have included as part of Python script functionality that allows for ClinTrials columns to be processed with non-ASCII character removal for use with MetaMap 2016 (which is not downloaded on Linux machine, btw...if for some reason, you want to use MetaMap 2016 on Linux, it must be downloaded and installed)
Otherwise, on Linux, MetaMap 2020 will be triggered with no text processing required (function: see select_cols_to_preprocess(), the "if-else" statement)

STEPS:
1) extracted data from ClinicalTrials.gov database (specifically, the following tables)
	- conditions [columns = ]
	- interventions [columns = ]
	- browse_conditions [columns = ]
	- browse interventions [columns = ]
2) Pre-process those columns to get non-ASCII char removed (optional, do not use for MetaMap 2020+). Get unique list of Conditions and Interventions and store in list to feed to MetaMap to get matching CUIs, to turn into CURIEs later
3) Start MetaMap SKR Part of Speech Tagger and WSD Word Sense Disambiguation
4) Feed unique Conditions and Interventions to MetaMap python wrapper (pymetamap), return top matching CUI for all concepts and link back to original queried diseases and interventions 
5) Filter out rows where CUI concept is an exact match (barring upper or lower case differences) for original queried disease or intervention (future releases will allow for matched concepts that may conceptually be correct (human subject matter expert), but may not be exact matches)
	- ex: 	diabetes type 2 is the orignal search query
			MetaMap returns: Type II Diabetes Mellitus
			Since "Type II Diabetes Mellitus" =/= "diabetes type 2" (ignoring case)
			We do not include this concept and corresponding relationship in the current version of the KG 
			
			However, if the original search query is "parkinson's disease",
			And MetaMap returns "Parkinson's Disease",
			This record is included since these two concepts are exact matches ignoring case
			
			Future releases will strive to accomodate conceptually correct, inexact matches
5) Create dataframes corresponding to required edges and nodes csv outputs
	edges file fields: subject, predicate, object, subject_name, object_name, category, nctid, nctid_curie
	nodes file fields: id, name, category
6) 



Description of files and outputs:
ClinTrials_KG_nodes_v01_2.csv = 
ClinTrials_KG_edges_v01_2.csv = 
ClinTrials_ETL.py = 
Clinical_Trials_playground.ipynb = 

metamapped_conditions.txt
metamapped_interventions.txt = 
readable_metamapped_conditions.txt
readable_metamapped_interventions.txt

ct_extract_env.yml 
Packages_to_exlude_from_yml.txt