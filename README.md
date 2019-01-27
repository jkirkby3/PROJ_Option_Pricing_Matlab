# Option Pricing PROJ Method (Exotic/Vanilla Options)
Option pricing (exotic/vanilla derivatives) based on an efficient and general Fourier transform pricing framework - the PROJ method (short for Frame Projection). The modules are organized by Model, and then by contract type. Each contract has a run script, which starts with "Script_", e.g. "Script_BarrierOptions.m".

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

<b>Acknowledgement:</b>
These pricing methods libraries have been built in collaboration with
<ul>
  <li> <a href="https://www.researchgate.net/profile/Justin_Kirkby"> Justin Lars Kirkby </a> </li>
  <li><a href="https://www.researchgate.net/profile/Duy_Nguyen125">Duy Nguyen </a> </li>
  <li><a href="https://www.researchgate.net/profile/Zhenyu_Cui"> Zhenyu Cui </a> </li>
  <li><a href="https://www2.isye.gatech.edu/people/faculty/Shijie_Deng/"> Shijie Deng </a> </li>
 </ul>

<b> Supporting Research Articles: </b>
<ul>
  <li> 
    <a href="https://www.researchgate.net/publication/271529024_Efficient_Option_Pricing_by_Frame_Duality_with_the_Fast_Fourier_Transform"> Efficient Option Pricing by Frame Duality with the Fast Fourier Transform, SIAM J. Financial Math (2015)</a>
    </li>
  <li> <a href=""> </a> </li>
</ul>


