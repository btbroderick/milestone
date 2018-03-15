Overview
--------

This App is used to predict a milestone event date for Clinical Trial
recruitment. The methods used in this App can be found at
[link](www.rstudio.com)

References
----------

[1](www.rstudio.com) Last, First M., and First M. Last. "Article Title."
Journal Title Series Volume. Issue (Year Published): Page(s). Print.

[2](www.rstudio.com) Last, First M., and First M. Last. "Article Title."
Journal Title Series Volume. Issue (Year Published): Page(s). Print.

Data
----

-   The Milestone App supports .csv, .tsv, .txt, .xlsx, .xls, .sas7bdat
    files.
-   Any uploaded data should contain 3 columns;
    1.  Days from the opening of the Clinical Trial a patient was
        enrolled
    2.  Days from opening of Clinical Trial patient had event (0 if no
        event)
    3.  An event indicator (0 or 1).

<!-- -->

    # A tibble: 10 x 3
       onstudy event_time event_indicator
         <int>      <int>           <int>
     1      82         86               1
     2      32          0               0
     3      63          0               0
     4     166          0               0
     5      56          0               0
     6     143          0               0
     7     163          0               0
     8     303        308               1
     9     321          0               0
    10     263        268               1
