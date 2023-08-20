## Instructions

This online tool is designed to analyse data and interpret the results of FECRT studies using a standardised and automated procedure.  Please follow the steps below to analyse your data.  You can return to these instructions at any time by clicking the "Instructions" tab above.

### Step 1:  Data entry

When you click the "Data" tab above, you will be presented with three options to enter your data:

- File upload:  select this option to upload your data in Excel and/or CSV format (see details on data formats below)
- Direct entry:  select this option to copy/paste your data directly into a web form
- Demonstration:  select this option to automatically generate some data for testing/demonstration purposes

#### File upload

When uploading data, it is possible to provide multiple CSV and/or Excel files, each specifying a different dataset to analyse.  Alternatively, multiple tabs within the same Excel sheet can be used for different datasets.  Each dataset should be formatted as follows:

- For paired data:  the two columns "PreTreatment" and "PostTreatment" must be present, and each should contain EPG of the corresponding animal pre- and post-treatment.  Any missing observations can be left blank (but these must be missing at random, and not due to e.g. systematic non-treatment of animals with 0 EPG at pre-treatment).
- For unpaired data:  the two columns "Group" and "EPG" must be present, where "Group" contains the text "Treatment" or "Control" and EPG gives the observed EPG.  Alternatively, "Control" and "Treatment" columns can be given, with blank cells used if the number of animals in each group does not match.

Any files/sheets not matching one of these options will be ignored.  The column names must be given in row 1 and the data must start in row 2, with no blank rows.  Do not include any other numbers (such as averages) in any row below the data columns.  The name of the group will be taken from the name of the CSV file or the sheet name for Excel files, so file and sheet names must be unique.  Note that all groups must have the same parameter options (i.e. multiplication factor, target efficacy and grey zone) - if this is not the case then you will need to divide your data into separate analyses accordingly.

#### Direct entry

For direct data entry, you should first select the study type (paired or unpaired, aka case/control), and then select the number of animals.  After clicking on the "Initialise data entry" button you can then enter your data into the cells below.  All data entered must be valid numbers - if an observation is missing then leave the corresponding cell blank.

#### Demonstration

If you want to quickly see what the tool can do without having to enter data, then you can use the demonstration option.  This randomly generates some data (for a paired study with 20 animals), and pre-fills the other options you need to run the analysis.


### Step 2:  Select Parameters

Once you have entered data, you can select the relevant parameters for your analysis using the "Parameters" tab.  If you are following the recommendations given in the 2023 WAAVP guideline for diagnosing anthelmintic resistance using the faecal egg count reduction test in ruminants, horses and swine (Kaplan et al, 2023) then you only need to select the following options:

- Species and anthelmintic:  the host/parasite species and anthelmintic that you have used for your study. Note that this only covers the situations covered by Kaplan et al. (2023) - a more flexible option for specifying custom efficacy targets will be added soon.
- Multiplication factor:  the counting sensitivity of the laboratory method used (e.g. 50 for McMaster). This must be a number greater than zero, and will be used to divide the specified EPG before analysis (which must result in a whole number i.e. count for each data point). Multiplication factor may differ between pre- and post-treatment (or treatment and control), but it is currently not possible for multiplication factor to vary within each observation set (as may be needed for some automated FEC methods).
- Country/region:  an optional text entry giving the country/region where the animals are kept, which will be included in the report if specified.
- Study identifier:  an optional text entry giving the study ID for reference purposes, which will be included in the report if specified.


### Step 3:  Results

Once you have entered data and provided valid parameters, you can run the analysis and obtain your results using the "Results" tab. First, click on the "Click to Calculate" button (if this button is not shown, then you should instead follow the instructions to correct the problem with the data/parameters that you entered). After a short wait, you should then see the headline result of your analysis (for each dataset provided) appear. You can then select to download a full report (including more information such as summary statistics, and comparisons of results generated by different analysis methods), either in Word or PDF format.

If you wish, you may go back and change parameters and then re-run the analysis.  You can also edit your data if you have entered this directly, or replace the files uploaded with new datasets.  Due to technical reasons, it is not currently possible to change the data entry method to direct data entry (or demonstration) after a file has already been uploaded (but you can refresh your browser to re-load the tool if you need to do this).


### About

This online tool was developed by Matt Denwood at the University of Copenhagen, with input from Ray Kaplan, Martin K. Nielsen, Bruno Levecke, Stig Milan Thamsborg, and Iain McKendrick.  We acknowledge funding from NordForsk (the DigiVet project on digitalisation of livestock data:  https://www.dcs.gla.ac.uk/~jenright/digivet_website/) and the Scottish Government Rural and Environment Science and Analytical Services Division (RESAS) Strategic Research Programme, for support of this work.  If you have comments/questions/suggestions, then please feel free to get in touch [by emailing Matt Denwood](md@sund.ku.dk).

For more information on the WAAVP guidelines see [Kaplan et al. (2023)](https://pubmed.ncbi.nlm.nih.gov/37121092/).  To cite this tool in peer-reviewed publications (and for more details on the underlying statistical principles) see [Denwood et al. (2023)](https://pubmed.ncbi.nlm.nih.gov/36621042/).  A version of this tool is hosted at [the fecrt.com website](https://www.fecrt.com).

### Version

This software is part of the [open-source bayescount package](https://github.com/ku-awdc/bayescount), version 1.1.0-1 (2023-08-20).
