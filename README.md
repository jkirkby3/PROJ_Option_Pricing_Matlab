# Option Pricing PROJ Method (Exotic/Vanilla Options)
Option pricing (exotic/vanilla derivatives) based on an efficient and general Fourier transform pricing framework - the PROJ method (short for Frame Projection). The modules are organized by Pricing Method, then by Model, and then by Contract Type. Each contract has a run script, which starts with "Script_", e.g. "Script_BarrierOptions.m". 
Monte Carlo and other pricing libraries are also provided to support R&D.

<b>Pricing methods supported:</b>
<ul>
  <li> PROJ (General Purpose Fourier Method) </li>
  <li> Monte Carlo </li>
  <li> Analytical </li>
  <li> Fourier (Carr-Madan, PROJ) </li>
  <li> PDE/Finite Difference </li>
  <li> Lattice/Tree </li>
</ul>  
  
<b>Contract types supported:</b>
<ul>
  <li> European Options </li>
  <li> Barrier Options (Single/Double barrier, with early excercise, and rebates) </li>
  <li> Asian Options (Discrete/Continuous)</li>
  <li> Discrete Variance Swaps, Variance/Volatility Options </li>
  <li> Bermudan/American early-exercise Options </li>
  <li> Parisian Options (Cumulative and resetting Parisian barrier options) </li>
  <li> Cliquets/Equity Indexed Annuities (Additive/Multiplicative)</li>
  <li> Forward Starting Options </li>
  <li> Step (Soft Barrier) Options </li>
  <li> Fader Options (To be added) </li>
  <li> Swing Options (To be added) </li>
  <li> Lookback/Hindsight Options (To be added) </li>
  <li> Credit default swaps (To be added) </li>
 </ul>
  
<b>Models supported:</b>
<ul>
  <li> Diffusions (Black-Scholes-Merton) </li>
  <li> Jump Diffusions (Merton Jump, Kou double exponential, Mixed-Normal)  </li>
  <li> General Levy processes (CGMY/KoBoL, Normal-Inverse-Gaussian (NIG), Variance Gamma, Meixner) </li>
  <li> SABR </li>
  <li> Stochastic Local Volatility </li>
  <li> Regime switching jump diffusions </li>
  <li> Time-changed processes </li>
  <li> Stochastic Volatility (Heston/Bates, Hull-White, 4/2, 3/2, alpha-hypergeometric, Jacobi, Schobel-Zhu, Stein-Stein, Scott, tau/2)   </li>
</ul>

<b>Acknowledgement:</b>
These pricing libraries have been built in collaboration with:
<ul>
  <li><a href="https://www.researchgate.net/profile/Justin_Kirkby"> Justin Lars Kirkby </a> </li>
  <li><a href="https://www.researchgate.net/profile/Duy_Nguyen125">Duy Nguyen </a> </li>
  <li><a href="https://www.researchgate.net/profile/Zhenyu_Cui"> Zhenyu Cui </a> </li>
  <li><a href="https://www2.isye.gatech.edu/people/faculty/Shijie_Deng/"> Shijie Deng </a> </li>
 </ul>

<b> Supporting Research Articles: </b>
<ul>
  <li> 
    <a href="https://www.researchgate.net/publication/271529024_Efficient_Option_Pricing_by_Frame_Duality_with_the_Fast_Fourier_Transform"> Efficient Option Pricing by Frame Duality with the Fast Fourier Transform. SIAM J. Financial Math (2015)</a>
    </li>
  <li> <a href="https://www.researchgate.net/publication/290607662_An_Efficient_Transform_Method_for_Asian_Option_Pricing">An Efficient Transform Method for Asian Option Pricing. SIAM J. Financial Math (2016) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/315888055_A_general_framework_for_discretely_sampled_realized_variance_derivatives_in_stochastic_volatility_models_with_jumps">A general framework for discretely sampled realized variance derivatives in stochastic volatility models with jumps. European J. Operational Research (2017) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/317056519_A_unified_approach_to_Bermudan_and_Barrier_options_under_stochastic_volatility_models_with_jumps"> A unified approach to Bermudan and Barrier options under stochastic volatility models with jumps. J. Econ. Dynamics and Control (2017)</a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/271529064_Static_Hedging_and_Pricing_of_Exotic_Options_With_Payoff_Frames"> Static Hedging and Pricing of Exotic Options With Payoff Frames. Mathematical Finance (2018) </a>
  </li>
  <li> <a href="https://www.researchgate.net/publication/320720280_American_and_Exotic_Option_Pricing_with_Jump_Diffusions_and_Other_Levy_Processes"> American and Exotic Option Pricing with Jump Diffusions and Other Levy Processes. J. Computational Finance (2018) </a>
  </li>
  <li> <a href="https://www.researchgate.net/publication/271529039_Robust_Barrier_Option_Pricing_by_Frame_Projection_Under_Exponential_Levy_Dynamics"> Robust Barrier Option Pricing by Frame Projection Under Exponential Levy Dynamics. Applied Mathematical Finance (2018) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/320264970_Robust_option_pricing_with_characteristic_functions_and_the_B-spline_order_of_density_projection"> Robust option pricing with characteristic functions and the B-spline order of density projection, J. Compuational Finance (2017) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/314260670_Equity-linked_annuity_pricing_with_cliquet-style_guarantees_in_regime-switching_and_stochastic_volatility_models_with_jumps"> Equity-linked annuity pricing with cliquet-style guarantees in regime-switching and stochastic volatility models with jumps. Insurance: Mathematics and Economics (2017)</a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/327411363_A_General_Framework_for_Time-Changed_Markov_Processes_and_Applications"> A General Framework for Time-Changed Markov Processes and Applications. European J. Operational Research (2018) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/324731726_A_General_Valuation_Framework_for_SABR_and_Stochastic_Local_Volatility_Models"> A General Valuation Framework for SABR and Stochastic Local Volatility Models. SIAM J. Financial Mathematics (2018) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/330401656_Continuous-Time_Markov_Chain_and_Regime_Switching_Approximations_with_Applications_to_Options_Pricing"> Continuous-Time Markov Chain and Regime Switching Approximations with Applications to Options Pricing. IMA Volumes on Mathematics (2019)</a> 
  </li>
  <li>
    <a href="https://smartech.gatech.edu/bitstream/handle/1853/59138/KIRKBY-DISSERTATION-2016.pdf"> Frame and Fourier Methods for Exotic Option Pricing and Hedging. Georgia Institute of Technology (2016). </a>
  </li>
</ul>


