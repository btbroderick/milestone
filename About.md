#### Overview

This App is used to predict a milestone event date for Clinical Trial
recruitment. The methods used in this App can be found at
[link](www.rstudio.com)



#### References

[1](www.rstudio.com) Last, First M., and First M. Last. "Article Title."
Journal Title Series Volume. Issue (Year Published): Page(s). Print.

[2](www.rstudio.com) Last, First M., and First M. Last. "Article Title."
Journal Title Series Volume. Issue (Year Published): Page(s). Print.



#### Data

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
     1     262        267               1
     2     301        306               1
     3     257        262               1
     4      38         41               1
     5      10          0               0
     6     344          0               0
     7     332        337               1
     8     151        156               1
     9     242        247               1
    10      80          0               0
