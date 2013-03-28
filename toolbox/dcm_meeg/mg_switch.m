function s = mg_switch(V)
% switching output s: determined by voltage (V) depedent magnesium blockade
% parameters as per Durstewitz, Seamans & Sejnowski 2000

s = 1.50265./(1 + 0.33*exp(-0.06.*V));
          