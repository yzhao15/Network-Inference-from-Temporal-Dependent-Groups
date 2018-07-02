# Network-Inference-from-Temporal-Dependent-Groups
Data:

chimp.xlsx -- Grouped data compiled from the Kibale Chimpanzee Project(https://kibalechimpanzees.wordpress.com/)

Columns    -- groups.
Rows       -- chimps. Last row indicate whether or not this group is identified at the same time of the previous group: 1-no/0-yes

Functions:

main_chimp -- main file
EM         -- Algorithm for the temporal-dependent hub model
E_step     -- E step of the temporal-dependent hub model
M_step     -- M step of the temporal-dependent hub model
HM_EM      -- Algorithm for the classical hub model
PlotGreyScale_try/PlotGreyScale_try_colorbar -- plot the estimated adjacency matrix with or without the colorbar
PlotGreyScale_nonSquare_try -- plot the data
cal_jaccard -- compute the Jaccard index
