Codes to extract paddle-to-paddle corrections:
1) parameter_extraction.cpp -> do the actual extraction
2) parameter_test.cpp       -> check that the parameters actually align the paddles
3) mean_time.cpp            -> Copy of 190228_band_time_spectra/mean_time.cpp, this time applying the paddle-to-paddle corrections extracted from 1)

We align all the bars in each layer to a reference bar in that layer. Then, as a second step, we align each layer to a reference layer.
To be more specific, we align all the bars from layer i to bar corresponding to (sector: 2, layer: i, component: 1).
Then, we align these reference bars (sector: 2, layer: i, component: 1) for each layer to bar (sector: 2, layer: 5, component: 1).
