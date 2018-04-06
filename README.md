# enzkin
ENZyme KINetics is the study of the chemical reactions that are catalysed by enzymes.<br/>
In enzyme kinetics the reaction rate is measured and the effects of
varying the conditions of the reaction investigated. Studying an enzyme
kinetics in this way can reveal the catalytic mechanism of this enzyme, its
role in metabolism, how its activity is controlled, and how a drug or a poison
might inhibit the enzyme.
Michaelis–Menten kinetics approximately describes the kinetics of many
enzymes. It is named after Leonor Michaelis and Maud Menten. This kinetic
model is relevant to situations where very simple kinetics can be assumed,
(i.e. there is no intermediate or product inhibition, and there is no
allostericity or cooperativity).
The Michaelis–Menten equation relates the initial reaction rate v0  to the
substrate concentration S. The corresponding graph is a rectangular
hyperbolic function; the maximum rate is described as Vmax (asymptote); the
concentration of substrate where the v0 is the half of Vmax is the
Michaelis-Menten costant (Km).
To determine the maximum rate of an enzyme mediated reaction, a series of
experiments is carried out with varying substrate concentration and the
initial rate of product formation is measured. 'Initial' here is taken to mean
that the reaction rate is measured after a relatively short time period,
during which complex builds up but the substrate concentration remains
approximately constant and the quasi-steady-state assumption will hold.
Accurate values for Km and Vmax can only be determined by non-linear
regression of Michaelis-Menten data.
The Michaelis-Menten equation can be linearized using several techniques.
ENZKIN uses 6 regression models (2 non-linear and 4 linear) to obtain the
kinetic parameters.

Syntax: 	enzkinout=enzkin(S,v)
     
    Inputs:
          S - data array of substrate concentrations
          v - data array of measured initial velocity
    Outputs:
          - Vmax and Km estimation by:
               ° Michaelis-Menten non linear regression
               ° loglog non linear regression
               ° Lineweaver-Burk linear regression
               ° Hanes-Woolf linear regression
               ° Eadie-Hofstee linear regression
               ° Scatchard linear regression
          - for the linear regressions, all regression data are summarized
          - Plots

The function requires another function of mine MYREGR. If it is not present on
the computer, enzkin will try to download it from FEX

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2010). Enzkin: a tool to estimate Michaelis-Menten kinetic
parameters
http://www.mathworks.com/matlabcentral/fileexchange/26653
