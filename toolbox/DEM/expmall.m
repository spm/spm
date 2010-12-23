function dx = expmall(J,f,t,EP)

EP(EP==1) = J;
EP(EP==2) = f;


dx = expm(EP*t);

