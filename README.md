# Option Pricing PROJ Method (Exotic/Vanilla Options)
Option pricing (exotic/vanilla derivatives) based on an efficient and general Fourier transform pricing framework - the PROJ method (short for Frame Projection). The modules are organized by Pricing Method, then by Model, and then by Contract Type. Each contract has a run script, which starts with "Script_", e.g. "Script_BarrierOptions.m". 
Monte Carlo and other pricing libraries are also provided to support R&D.

<b>Pricing methods supported:</b>
<ul>
  <li> PROJ (General Purpose Fourier Method) </li>
  <li> CTMC Approximation </li>
  <li> Monte Carlo </li>
  <li> Analytical </li>
  <li> Fourier (PROJ, Carr-Madan, CONV, Lewis, COS, Mellin Transform, Hilbert Transform) </li>
  <li> PDE/Finite Difference </li>
  <li> Lattice/Tree </li>
</ul>  
  
<b>Models supported:</b>
<ul>
  <li> Diffusions (Black-Scholes-Merton) </li>
  <li> Multi-Dimensional Diffusions (Black-Scholes Multi-Asset) </li>
  <li> Jump Diffusions (Merton Jump, Kou double exponential, Mixed-Normal)  </li>
  <li> General Levy processes (CGMY/KoBoL, Normal-Inverse-Gaussian (NIG), Variance Gamma, Meixner, FMLS, TS, Bilateral Gamma) </li>
  <li> SABR </li>
  <li> Stochastic Local Volatility (SLV) </li>
  <li> Regime switching jump diffusions </li>
  <li> Time-changed processes </li>
  <li> Stochastic Volatility (Heston, Hull-White, 4/2, 3/2, alpha-hypergeometric, Jacobi, Schobel-Zhu, Stein-Stein, Scott, tau/2)   </li>
  <li> Stochastic Volatility With Jumps (e.g. Bates, HKDE) </li>
</ul>

<b>Contract types supported (single underlying):</b>
<ul>
  <li> European Options </li>
  <li> Barrier Options (Single/Double barrier, and rebates) </li>
  <li> Return Barrier Options </li>
  <li> Asian Options (Discrete/Continuous)</li>
  <li> Discrete Variance Swaps, Variance/Volatility Options </li>
  <li> Bermudan/American early-exercise Options </li>
  <li> Parisian Options (Cumulative and resetting Parisian barrier options) </li>
  <li> Cliquets/Equity Indexed Annuities (Additive/Multiplicative)</li>
  <li> Equity Linked Death Benefits / Guaranteed Minimum Death Benefits (GMDB) </li>
  <li> Forward Starting Options </li>
  <li> Step (Soft Barrier) Options </li>
  <li> Lookback/Hindsight Options </li>
  <li> Credit default swaps / default probabilities </li>
  <li> Swing Options (Fixed Rights, Linear Recovery & Constant Recovery) </li>
  <li> Fader/Range-Accrual Options  </li>
  <li> Multi-Dimensional Payoffs, European/Bermudan/Barrier (Spread, Exchange, Best/Worst-of, Basket) </li>
  <li> Risk Measures suchs as Expected Shortfall and VaR computations </li>
 </ul>

<b>Contract types supported (multi underlying):</b>
<ul>
  <li> European / Barrier / Bermudan Options </li>
  <li> Spread, Exchange, Best-of, Worst-of, Basket (Geometric/Arthmetic) </li>
</ul>

<b>Acknowledgement:</b>
These pricing libraries have been built in collaboration with:
<ul>
  <li><a href="https://www.researchgate.net/profile/Justin_Kirkby"> Justin Lars Kirkby </a> </li>
  <li><a href="https://www.researchgate.net/profile/Duy_Nguyen125">Duy Nguyen </a> </li>
  <li><a href="https://www.researchgate.net/profile/Zhenyu_Cui"> Zhenyu Cui </a> </li>
  <li><a href="https://www.researchgate.net/profile/Zhimin_Zhang3"> Zhimin Zhang </a> </li>
  <li><a href="https://www.researchgate.net/profile/Shi-jie_Deng"> Shijie Deng </a> </li>
  <li><a href="https://www.researchgate.net/profile/Jean_Philippe_Aguilar"> Jean-Philippe Aguilar </a> </li>
</ul>

<b> Supporting Research Articles: </b>

<b> I) Levy Models, Jump Diffusions, Black Scholes </b>
<ul>
  <li> 
    <a href="https://www.researchgate.net/publication/271529024_Efficient_Option_Pricing_by_Frame_Duality_with_the_Fast_Fourier_Transform"> Efficient Option Pricing by Frame Duality with the Fast Fourier Transform. SIAM J. Financial Math (2015)</a>
  </li>
  <li> <a href="https://www.researchgate.net/publication/290607662_An_Efficient_Transform_Method_for_Asian_Option_Pricing">An Efficient Transform Method for Asian Option Pricing. SIAM J. Financial Math (2016) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/271529064_Static_Hedging_and_Pricing_of_Exotic_Options_With_Payoff_Frames"> Static Hedging and Pricing of Exotic Options With Payoff Frames. Mathematical Finance (2018) </a>
  </li>
  <li> <a href="https://www.researchgate.net/publication/320720280_American_and_Exotic_Option_Pricing_with_Jump_Diffusions_and_Other_Levy_Processes"> American and Exotic Option Pricing with Jump Diffusions and Other Levy Processes. J. Computational Finance (2018) </a>
  </li>
  <li> <a href="https://www.researchgate.net/publication/271529039_Robust_Barrier_Option_Pricing_by_Frame_Projection_Under_Exponential_Levy_Dynamics"> Robust Barrier Option Pricing by Frame Projection Under Exponential Levy Dynamics. Applied Mathematical Finance (2018) </a> 
  </li>
  <li> <a href="https://www.researchgate.net/publication/320264970_Robust_option_pricing_with_characteristic_functions_and_the_B-spline_order_of_density_projection"> Robust option pricing with characteristic functions and the B-spline order of density projection, J. Compuational Finance (2017) </a> 
  </li>
  <li>
    <a href="https://www.researchgate.net/publication/334822473_Valuing_equity-linked_death_benefits_in_general_exponential_Levy_models"> Valuing Equity-Linked Death Benefits in General Exponential Levy Models. J. Comput. and Appl. Math. (2019). </a>
  </li>
  <li> <a href="https://www.researchgate.net/publication/336281657_Swing_Option_Pricing_by_Dynamic_Programming_with_B-Spline_Density_Projection"> Swing Option Pricing By Dynamic Programming with B-Spline Density Projection, IJTAF, Forthcoming (2020)</a> 
  </li>
    <li>
    <a href="https://smartech.gatech.edu/bitstream/handle/1853/59138/KIRKBY-DISSERTATION-2016.pdf"> Frame and Fourier Methods for Exotic Option Pricing and Hedging. Georgia Institute of Technology (2016). </a>
  </li>
</ul>

<b> II) Stochastic Volatility, Markov Chains, and Regime Switching </b>
<ul>
    <li> <a href="https://www.researchgate.net/publication/315888055_A_general_framework_for_discretely_sampled_realized_variance_derivatives_in_stochastic_volatility_models_with_jumps">A general framework for discretely sampled realized variance derivatives in stochastic volatility models with jumps. European J. Operational Research (2017) </a> 
  </li>
    <li> <a href="https://www.researchgate.net/publication/317056519_A_unified_approach_to_Bermudan_and_Barrier_options_under_stochastic_volatility_models_with_jumps"> A unified approach to Bermudan and Barrier options under stochastic volatility models with jumps. J. Econ. Dynamics and Control (2017)</a> 
  </li>
    <li> <a href="https://www.researchgate.net/publication/314260670_Equity-linked_annuity_pricing_with_cliquet-style_guarantees_in_regime-switching_and_stochastic_volatility_models_with_jumps"> Equity-linked annuity pricing with cliquet-style guarantees in regime-switching and stochastic volatility models with jumps. Insurance: Mathematics and Economics (2017)</a> 
  </li>
    <li> <a href="https://www.researchgate.net/publication/330401656_Continuous-Time_Markov_Chain_and_Regime_Switching_Approximations_with_Applications_to_Options_Pricing"> Continuous-Time Markov Chain and Regime Switching Approximations with Applications to Options Pricing. IMA Volumes on Mathematics (2019)</a> 
  </li>
    <li> <a href="https://www.researchgate.net/publication/334716223_Full-fledged_SABR_Through_Markov_Chains?_sg=wav6ifhPa8HmCsvRzHVqYYPU2VHZKMTyP-1ZX_eeuqYZz5cpfKqZ0OCTODC9Ci1aY8j99amKGjbwZnaf1q1k2cTmLdIfxamtOAx_pXs8.W0biWEvGq-ILKu2DgzAI35-BBXMZp3bN1jBLDKKfSg_FgFd9ci8xXqXQIKbA5UoPE6sUA9GrpH8ByrP8-Xx1aA"> Full-Fledged SABR Through Markov Chains, Wilmott (2019) </a> 
  </li>
</ul>

<b> III) Stochastic Local Volatility (SABR, Quadratic SLV, etc) </b>
<ul>
  <li> <a href="https://www.researchgate.net/publication/324731726_A_General_Valuation_Framework_for_SABR_and_Stochastic_Local_Volatility_Models"> A General Valuation Framework for SABR and Stochastic Local Volatility Models. SIAM J. Financial Mathematics (2018) </a> 
  </li>
</ul>

<b> IV) Time-Changed Processes</b>
<ul>
    <li> <a href="https://www.researchgate.net/publication/327411363_A_General_Framework_for_Time-Changed_Markov_Processes_and_Applications"> A General Framework for Time-Changed Markov Processes and Applications. European J. Operational Research (2018) </a> 
  </li>
</ul>

<b> V) Multi-Dimensional Diffusions </b>
<ul>
  <li>  <a href="https://www.researchgate.net/publication/342174203_A_general_continuous_time_Markov_chain_approximation_for_multi-asset_option_pricing_with_systems_of_correlated_diffusions"> A general continuous time Markov chain approximation for multi-asset option pricing with systems of correlated diffusions (2020) </a> 
  </li>
</ul>


