# PROJ_Option_Pricing_Matlab
Option pricing libraries (Matlab) based on the PROJ method (short for Frame Projection), an efficient and general Fourier transform based option pricing framework for Vanilla and Exotic options. The modules are organized by Model, and then by contract type. Each contract has a run script, which starts with "Script_", e.g. "Script_BarrierOptions.m".

<b>Contract types suppoerted:</b>
<ul>
  <li> European Options </li>
  <li> Barrier Options (Single/Double barrier, with early excercise, and rebates) </li>
  <li> Asian Options (Discrete/Continuous)</li>
  <li> Discrete Variance Swaps, Variance/Volatility Options </li>
  <li> Bermudan/American early-exercise Options </li>
  <li> Parisian Options (Cumulative and resetting Parisian barrier options) </li>
  <li> Cliquets/Equity Indexed Annuities (Additive/Multiplicative)</li>
  <li> Forward Starting Options </li>
  <li> Step/Fader Options (To be added) </li>
  <li> Swing Options (To be added) </li>
  <li> Lookback/Hindsight Options (To be added) </li>
  <li> Credit default swaps (To be added) </li>
 </ul>
  
<b>Models supported:</b>
<ul>
  <li> Diffusions (Black-Scholes-Merton) </li>
  <li> Jump Diffusions (Merton Jump, Kou double exponential, Mixed-Normal)  </li>
  <li> General Levy processes (CGMY/KoBoL, Meixner, Variance Gamma) </li>
  <li> Regime switching jump diffusions </li>
  <li> Time-changed processes </li>
  <li> Stochastic Volatility (Heston/Bates, Hull-White, 4/2, 3/2, alpha-hypergeometric, Jacobi, Schobel-Zhu, Stein-Stein, Scott, tau/2)   </li>
</ul>

