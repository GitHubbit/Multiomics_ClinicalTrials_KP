Note: had started to do this in R at first (folder: ClinTrial_ETL_R), but preference of 
others was Python (script is called ClinicalTrials_ETL.py)

Downloaded static database from AACT (Clinical Trials.gov, https://aact.ctti-clinicaltrials.org/download)
Unpack it using Postgresql software
store as zip in current folder
access it in pandas data tables in python (folder: ClinTrial_ETL_Python)

created dir called metamap to shove all MetaMap related progs into
mkdir "metamap"

Download MetaMap software mapping tool from to run locally:
https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/run-locally/MainDownload.html
	- downloaded public_mm_linux_main_2020.tar.bz2 locally
	- then copied it from local to hypatia server via scp
	- run bunzip -c public_mm_linux_main_2020.tar.bz2 | tar xvf - ...
Follow successive installation instructions from NIH NLM:
https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/documentation/Installation.html
You must have Java 1.8
MetaMap requires 2 servers (SKR Part of Speech Tagger and WSD Word Sense Disambiguation) for our purposes:
We start and stop those servers by calling it from the Python script

Git clone wrapper to run MetaMap in Python script (ClinicalTrials_ETL.py)
git clone https://github.com/AnthonyMRios/pymetamap

MetaMap 2016 requires removal of non-ASCII char. MetaMap 2020 does not. Original script was developed on Mac/OSX, for which there's no MetaMap 2020, only MetaMap 2016 (at time of development, June 2022)
I have included as part of Python script functionality that allows for ClinTrials columns to be processed with non-ASCII character removal for use with MetaMap 2016 (which is not downloaded on Linux machine, btw...if for some reason, you want to use MetaMap 2016 on Linux, it must be downloaded and installed)
Otherwise, on Linux, MetaMap 2020 will be triggered with no text processing required (function: see select_cols_to_preprocess(), the "if-else" statement)


