#### Overview

This App can be used to predict a milestone event time (date of certain
number of events occurs) for clinical trial monitoring. Please see
reference below for the methods used behind this App.

#### Data

-   The Milestone App supports .csv, .tsv, .txt, .xlsx, .xls, .sas7bdat
    files.
-   Any uploaded data should contain 3 columns;
    1.  Column 1: Days of patient enrollment since first patient
        enrolled (i.e. date of patient enrolled minus the date of first
        patient enrolled). Note, the first patient enrolled for the
        study would have 0 day.
    2.  Column 2: Days of patient last known event status time since
        first patient enrolled (i.e. date of patient last know event
        status minus the date of first patient enrolled).
    3.  Column 3: Event indicator (1=event, 0=censored).

#### Disclaimer

Mayo Clinic does not save any of the data or output.
