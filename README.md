# monier2018reconstruction codes

[Reconstruction of partially sampled multi-band images – Application to STEM-EELS imaging](https://ieeexplore.ieee.org/document/8447232).

These codes aim at reproducing the paper results. Note that the random generator seeds may be different from those used in the paper. Performances may slightly differ also.

**All codes were written under Matlab 2017a**.

## Generating the synthetic dataset

To generate the synthetic data, please run `synth_data/synthetic/generate_synth_data.m`. This will generate the data file at location `synth_data/synthetic_data.mat`.

## Performance w.r.t. the regularization parameters
 
Please run `1. Regularization parameter/search_parameter_SSS.m` and `1. Regularization parameter/search_parameter_SNN.m` for 3S and S2N methods.

Note that the S2N codes are particularly time-consuming !

## Performance w.r.t. noise level and sampling ratio

Please run `2. pix sigma performances/pix_sigma_SSS.m` and `2. pix sigma performances/pix_sigma_SNN.m` for 3S and S2N methods.

## Reconstruction vs. denoising w.r.t. an unmixing task

All codes for these experiments are located in the `3. Reconstruction performances (unmixing)` folder. Please run codes in the `0. Create data/` folder to create the synthetic data for full sampling and partial sampling protocols. This generates the data to be analyzed in the `1. Data2use/` folder.

Last, run the `2. Processing code/process_performances.m` to analyze the data. Everything is finally saved under the `3. Final figures/` folder.

## The real-data example

All codes for these experiments are located in the `4. Real data` folder. Please run codes in the `0. Create data/` folder to create the denoised and reconstructed data. This generates the data to be analyzed in the `1. Data2use/` folder.

Last, run the `2. Processing code/process_performances.m` to analyze the data. Everything is finally saved under the `3. Final figures/` folder.

Please note that the unmixing procedure assumes that the sample contains 5 components, but only the 3 most significative maps and spectra are kept. During the execution of the `2. Processing code/process_performances.m` code, the Matlab console will ask for the sorted indexes of the 3 significative maps to keep (e.g. `[1 5 3]`). This permits the code to re-arrange the maps and spectra in a correct way.

## Author and license

These codes were written by [Etienne Monier](https://etienne-monier.github.io/) and are distributed under the MIT license.

