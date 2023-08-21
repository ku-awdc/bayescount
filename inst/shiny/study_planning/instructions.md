## Instructions

This online tool allows you to determine the minimum viable sample size for a FECRT study, given a set of user-supplied parameters such as expected efficacy and the size of the grey zone. This allows some additional flexibility around the standard guidelines provided by Kaplan et al. (2023).  You can return to these instructions at any time by clicking the "Instructions" tab above.


### Step 1:  Select Parameters

You can select the relevant parameters for your analysis using the "Parameters" tab.  If you are designing a FECRT study according to the guidelines given in the 2023 WAAVP guideline for diagnosing anthelmintic resistance using the faecal egg count reduction test in ruminants, horses and swine (Kaplan et al, 2023) then you only need to select the following options:

- Host/parasite species  the host/parasite species that you have used for your study. Note that this only covers the situations covered by Kaplan et al. (2023) - a more flexible option for specifying custom efficacy targets will be added soon.
- Target efficacy:  the expected arithmetic mean efficacy for the anthelmintic used (in %).
- Non-inferiority margin:  the desired width of the grey zone (in % points) - changing this updates the lower efficacy target shown, for reference.
- Expected pre-treatment EPG:  the expected (arithmetic) mean EPG in untreated animals (i.e. either pre-treatment or control animals).
- Multiplication factor:  the counting sensitivity of the laboratory method used (e.g. 50 for McMaster). This must be a number greater than zero, and will be used to divide the specified EPG before analysis.
- Maximum sample size:  the largest sample size to consider for the calculation (this should be in excess of the largest feasible sample size in practice).


### Step 2:  Run Calculation and View Results

Once you have provided valid parameters, you can run the sample size calculation using the "Calculation" tab.  After a short delay, some feedback text will be shown to indicate either that a suitable sample size was determined, or that the maximum sample size specified was too small.  If a suitable sample size was not found, then you should go back and change parameters and then re-run the calculation.  Otherwise, the minimum sample size determined to achieve a statistical power of at least 80% for each test (resistance and susceptibility) will be shown along with a graphical illustration of the estimated statistical power as a function of sample size.


### About

This online tool was developed by Matt Denwood at the University of Copenhagen, with input from Ray Kaplan, Martin K. Nielsen, Bruno Levecke, Stig Milan Thamsborg, and Iain McKendrick.  We acknowledge funding from NordForsk (the DigiVet project on digitalisation of livestock data:  https://www.dcs.gla.ac.uk/~jenright/digivet_website/) and the Scottish Government Rural and Environment Science and Analytical Services Division (RESAS) Strategic Research Programme, for support of this work.  If you have comments/questions/suggestions, then please feel free to get in touch [by emailing Matt Denwood](md@sund.ku.dk).

For more information on the WAAVP guidelines see [Kaplan et al. (2023)](https://pubmed.ncbi.nlm.nih.gov/37121092/).  To cite this tool in peer-reviewed publications (and for more details on the underlying statistical principles) see [Denwood et al. (2023)](https://pubmed.ncbi.nlm.nih.gov/36621042/).  A version of this tool is hosted at [the fecrt.com website](https://www.fecrt.com).

### Version

This software is part of the [open-source bayescount package](https://github.com/ku-awdc/bayescount), version 1.1.0 (2023-08-21).
