import numpy as np
import matplotlib.pyplot as plt

# Percorso del file di testo
file_path = 'output'

# Caricamento dati
# Si presume che i numeri complessi siano nel formato 'a+bj' o 'a-bj'
dati = np.loadtxt(file_path, dtype=complex)

# Estrazione delle colonne
t    = dati[:, 0].real         # Prima colonna (reale)
env  = dati[:, 1].real         # Seconda colonna (reale)
psi0 = dati[:, 2]              # Terza colonna (complessa)
psi1 = dati[:, 3]              # Quarta colonna (complessa)
psi2 = dati[:, 4]              # Quinta colonna (complessa)
psi3 = dati[:, 5]              # Sesta colonna (complessa)

fig , ax = plt.subplots(1 , 1 , figsize = (8,8) ,  dpi = 300)
#append data to the figure
ax.plot(t , abs(psi0)[:]**2 , label = r"State $|0\rangle$")
ax.plot(t , abs(psi1)[:]**2 , label = r"State $|1\rangle$")


ax.plot(t , abs(psi2)[:]**2 , label = r"State $|2\rangle$")
ax.plot(t , abs(psi3)[:]**2 , label = r"State $|3\rangle$")

ax.plot(np.nan , ls = "--", c = "grey", alpha = 0.7,  label = "Envelope function" )
ax.legend(loc = 0)

ax1 = ax.twinx()
ax1.plot(t , env , ls = "--", c = "grey", alpha = 0.7)

#personalization
ax.set_xlabel(r"Time $(\mu s)$")
ax.set_ylabel("$|\Psi|^2$")
ax1.set_ylabel("Envelope intensity")

#save fiugre
plt.savefig("output.png")

