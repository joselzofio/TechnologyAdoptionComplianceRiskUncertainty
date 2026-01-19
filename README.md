# Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty

Version 4.0.0

This set of two MATLAB scripts calculates the investment thresholds, along with optimal emissions and violation levels from the simulations reported in Arguedas, C., Peinado, F., and Zofío, J.L. (2026) "Incentives for Green Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty", Resource and Energy Economics, Environmental and Resource Economics, 89(6), art. nº 6.

The first script (1.Technology_Adoption_Compliance_Risk_with_Tech_Certainty_10_10_2025.m) implements the simulations of the baseline model reported in section 4.2 of the article without technological uncertainty about the new technology.

The second script (2.Technology_Adoption_Compliance_Risk_with_Tech_Uncertainty_10_10_2025.m) implements the simulations reported in section 4.2 of the article enhancing the baseline model with technological uncertainty about the efficiency of the new technology, which may have either higher or lower abatement costs than the old technology.   

Both scripts have been run in version R2024b of MATLAB and use the function 'vpasolve' from the Symbolic Math Toolbox to determine optimal emissions. For the case of technological certainty and partial technological uncertainty see equation (2) of the article, minimizing firms' expected disutility in terms of compliance and non-compliance (equation (1)).  

## Usage

To run the scripts, add the m-file to the MATLAB path.

The simulations can be performed with the different parametrizations of the inspection probabilities and fines for non-compliance presented in Table 1 of the article. Users can adjust all the functions and parameters used in the simulations, including the abatement cost functions, the general hyperbolic disutility function (gamma, mu, and nu), tax value (tau), uncertainty about the efficiency of the new technology (alpha), and the inspection probability (pi) and fine functions. 

The scripts guide the user step by step through the process of solving the optimization problems presented in the article.

## NOTE

In the accompanying compressed (zip) file, users can find the scripts to run the simulations of a previous draft of the study: Arguedas, C., Peinado, F., and Zofío, J.L. (2023) "Incentives for Green Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty", Dptm. of Economics, Universidad Autónoma de Madrid. The working paper is available at SSRN: https://ssrn.com/abstract=5200584 or http://dx.doi.org/10.2139/ssrn.5200584. A detailed description of those simulations is presented in the included Working_Paper_README.md file. 


## Authors

Carmen Arguedas <br>
Fernando Peinado <br>
José L. Zofío


