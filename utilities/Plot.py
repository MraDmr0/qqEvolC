import numpy as np
import matplotlib.pyplot as plt
import sys

# Percorso del file di testo
prefix = str(sys.argv[1])
file_path = prefix + ".txt"

# Caricamento dati
dati = np.loadtxt(file_path, dtype=complex)

# Estrazione colonne
t = dati[:, 0].real      # Prima colonna (reale)
env = dati[:, 1].real    # Seconda colonna (reale)
psis = dati[:, 2:]       # Tutte le colonne successive (stati) come array 2D

fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=300)

# Ciclo per plottare tutti gli stati
for i, psi in enumerate(psis.T):  # psis.T per iterare sulle colonne
    ax.plot(t, abs(psi)**2, label=fr"State $|{i}\rangle$")

# Envelope function
ax.plot(np.nan, ls="--", c="grey", alpha=0.7, label="Envelope function")
ax.legend(loc=0)

# Secondo asse y per l'envelope
ax1 = ax.twinx()
ax1.plot(t, env, ls="--", c="grey", alpha=0.7)

# Personalizzazione assi
ax.set_xlabel(r"Time $(\mu s)$")
ax.set_ylabel(r"$|\Psi|^2$")
ax1.set_ylabel("Envelope intensity")

# Salvataggio figura
outfile = prefix + ".png"
plt.savefig(outfile)

