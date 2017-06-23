%%%% 06/23/2017
%%%% locate index in the saved wcoh variable
%%%% wcoh is an 1-d array, which saved the channel pair functional
%%%% connection metrics. This function is used to find the index on that
%%%% array based on the index of channels.
%%%% wcoh was saved as an 1-d array for saving disk space.

function ind_wcoh = lz_ind_loc_wcoh(iCh, jCh, nCh)

ind_wcoh = jCh - iCh + .5 * ( (2*nCh-iCh) * (iCh-1) );