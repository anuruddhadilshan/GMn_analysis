# How to use the GMn analysis scripts by Anuruddha Rathnayake - 01/04/2022

1. Run the 'define_tightcuts.C' script - this will show reasonable thresholds that can be applied to each variable in the root tree.
2. Run 'pasrse_gmn_rootfiles.C' to parse out obvious bad events and make the size of our analysis root files smaller, by taking the cut thresholds from step 1 as the inputs.
3. Run 'elastix_4.C' to make delta plots.

--
More steps to be included as the new analysis steps get added.