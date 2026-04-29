import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 🧩 PARAMÈTRES UTILISATEUR
# ============================================================

# Limites en x (communes aux deux plots)
# Mets None pour laisser le code choisir automatiquement
X_LIMITS = (0, 16)

# Limites en y (communes aux deux plots)
# Mets None pour auto
Y_LIMITS = (-100, 100)
# ============================================================


def load_csv(filename):
    """Lit un CSV et renvoie (t, data, noms_de_colonnes)."""
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    with open(filename, "r") as f:
        columns = f.readline().strip().split(",")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    t = data[:, 0]
    return t, data, columns


def set_limits(ax, datas, cols, xlim=None, ylim=None):
    """Fixe les limites des axes à partir des données ou des valeurs utilisateur."""
    # --- X ---
    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        all_x = np.concatenate([d[:, 0] for d in datas if d is not None])
        ax.set_xlim(np.nanmin(all_x), np.nanmax(all_x))

    # --- Y ---
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        all_y = np.concatenate([d[:, 1:].ravel() for d in datas if d is not None])
        ax.set_ylim(np.nanmin(all_y), np.nanmax(all_y))


def style_axes(ax, title=None, ylabel=None, xlabel=None):
    """Applique le style visuel standard."""
    if title:
        ax.set_title(title)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.set_axisbelow(True)
    ax.grid(True, linewidth=0.8, zorder=0)


# -------------------- Chargement des fichiers --------------------
# LIGNE 1
t1_L1, data1_L1, cols1_L1 = load_csv("solutions.csv")

# LIGNE 2
t1_L2, data1_L2, cols1_L2 = load_csv("./solutions_témoin/WED1.csv")

# -------------------- Figure 2x1 --------------------
fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, sharey=True)

# ------------------ PLOT 1 ------------------
for i, col in enumerate(cols1_L1[1:9], start=1):
    axes[0].scatter(t1_L1, data1_L1[:, i], s=2, label=col, zorder=3)
style_axes(axes[0], ylabel="Nathan - Solutions")
axes[0].legend(fontsize=8)

# ------------------ PLOT 2 ------------------
for i, col in enumerate(cols1_L2[1:], start=1):
    axes[1].scatter(t1_L2, data1_L2[:, i], s=2, label=col, zorder=3)
style_axes(axes[1], xlabel="x", ylabel="Jonathan - Solutions")
axes[1].legend(fontsize=8)

# ------------------ Limites communes ------------------
set_limits(axes[0], [data1_L1, data1_L2], [cols1_L1, cols1_L2], xlim=X_LIMITS, ylim=Y_LIMITS)
set_limits(axes[1], [data1_L1, data1_L2], [cols1_L1, cols1_L2], xlim=X_LIMITS, ylim=Y_LIMITS)

plt.tight_layout(h_pad=1.2)
plt.show()


## ============================================================
## 🧩 PARAMÈTRES UTILISATEUR
## ============================================================
#
## Limites en x (communes à tous les plots)
#X_LIMITS = [
#    (0, 5),   # Ligne 1
#    (0, 15),      # Ligne 2 -> auto d'après les données de la ligne
#    (0, 150),    # Ligne 3
#]
#
## Limites en y : une par ligne (gauche et droite synchronisées)
## Format : [(ymin, ymax), (ymin, ymax), (ymin, ymax)]
#Y_LIMITS = [
#    (-10, 10),     # Ligne 1
#    (-100, 100),     # Ligne 2
#    (-100, 100),     # Ligne 3
#]
## Mets None à la place pour laisser le code calculer automatiquement :
## Y_LIMITS = [None, None, None]
## ============================================================
#
#
#def load_csv(filename):
#    """Lit un CSV et renvoie (t, data, noms_de_colonnes)."""
#    data = np.loadtxt(filename, delimiter=",", skiprows=1)
#    with open(filename, "r") as f:
#        columns = f.readline().strip().split(",")
#
#    if data.ndim == 1:
#        data = data.reshape(1, -1)
#
#    t = data[:, 0]
#    return t, data, columns
#
#
#def set_row_ylims(ax_left, ax_right, datas, cols, manual_limit=None):
#    """Fixe des ylims identiques pour la ligne (2 axes) à partir des données ou des valeurs utilisateur."""
#    if manual_limit is not None:
#        ymin, ymax = manual_limit
#        ax_left.set_ylim(ymin, ymax)
#        ax_right.set_ylim(ymin, ymax)
#        return
#
#    mins, maxs = [], []
#    for data, c in zip(datas, cols):
#        if data is None:
#            continue
#        y = data[:, 1:]
#        if y.size:
#            mins.append(np.nanmin(y))
#            maxs.append(np.nanmax(y))
#    if mins and maxs:
#        y_min, y_max = float(np.min(mins)), float(np.max(maxs))
#        pad = 0.05 * (y_max - y_min) if y_max != y_min else 1.0
#        ax_left.set_ylim(y_min - pad, y_max + pad)
#        ax_right.set_ylim(y_min - pad, y_max + pad)
#
#
## -------------------- Chargement des fichiers --------------------
## LIGNE 1
#t1_L1, data1_L1, cols1_L1 = load_csv("solution.csv")
#t2_L1, data2_L1, cols2_L1 = load_csv("./solutions_témoin/JED1.csv")
#
## LIGNE 2
#t1_L2, data1_L2, cols1_L2 = load_csv("solution.csv")
#t2_L2, data2_L2, cols2_L2 = load_csv("./solutions_témoin/JED1.csv")
#
## LIGNE 3
#t1_L3, data1_L3, cols1_L3 = load_csv("solution.csv")
#t2_L3, data2_L3, cols2_L3 = load_csv("./solutions_témoin/JED1.csv")
#
## -------------------- Figure 3x2 --------------------
#fig, axes = plt.subplots(3, 2, figsize=(12, 6), sharex='col', sharey='row')
#
#
#def style_axes(ax, title=None, ylabel=None, xlabel=None):
#    if title:
#        ax.set_title(title)
#    if ylabel:
#        ax.set_ylabel(ylabel)
#    if xlabel:
#        ax.set_xlabel(xlabel)
#    ax.set_axisbelow(True)
#    ax.grid(True, linewidth=0.8, zorder=0)
#
#
## ------------------ LIGNE 1 ------------------
#for i, col in enumerate(cols1_L1[1:], start=1):
#    axes[0, 0].scatter(t1_L1, data1_L1[:, i], s=2, label=col, zorder=3)
#style_axes(axes[0, 0], title="Nathan - Solution", ylabel="Solutions")
#axes[0, 0].legend(fontsize=8)
#
#for i, col in enumerate(cols2_L1[1:], start=1):
#    axes[0, 1].scatter(t2_L1, data2_L1[:, i], s=2, label=col, zorder=3)
#style_axes(axes[0, 1], title="Jonathan - Solution")
#axes[0, 1].legend(fontsize=8)
#
#set_row_ylims(axes[0, 0], axes[0, 1], [data1_L1, data2_L1], [cols1_L1, cols2_L1], manual_limit=Y_LIMITS[0])
#
## ------------------ LIGNE 2 ------------------
#for i, col in enumerate(cols1_L2[1:], start=1):
#    axes[1, 0].scatter(t1_L2, data1_L2[:, i], s=2, label=col, zorder=3)
#style_axes(axes[1, 0], ylabel="Ligne 2")
#axes[1, 0].legend(fontsize=8)
#
#for i, col in enumerate(cols2_L2[1:], start=1):
#    axes[1, 1].scatter(t2_L2, data2_L2[:, i], s=2, label=col, zorder=3)
#style_axes(axes[1, 1])
#axes[1, 1].legend(fontsize=8)
#
#set_row_ylims(axes[1, 0], axes[1, 1], [data1_L2, data2_L2], [cols1_L2, cols2_L2], manual_limit=Y_LIMITS[1])
## ------------------ LIGNE 3 ------------------
#for i, col in enumerate(cols1_L3[1:], start=1):
#    axes[2, 0].scatter(t1_L3, data1_L3[:, i], s=2, label=col, zorder=3)
#style_axes(axes[2, 0], ylabel="Ligne 3", xlabel="x")
#axes[2, 0].legend(fontsize=8)
#
#for i, col in enumerate(cols2_L3[1:], start=1):
#    axes[2, 1].scatter(t2_L3, data2_L3[:, i], s=2, label=col, zorder=3)
#style_axes(axes[2, 1], xlabel="x")
#axes[2, 1].legend(fontsize=8)
#
#set_row_ylims(axes[2, 0], axes[2, 1], [data1_L3, data2_L3], [cols1_L3, cols2_L3], manual_limit=Y_LIMITS[2])
#
## ------------------ Réglages communs ------------------
#for r in range(3):  # 3 lignes
#    # marges
#    axes[r, 0].margins(x=0.02, y=0.05)
#    axes[r, 1].margins(x=0.02, y=0.05)
#
#    # xlims
#    if isinstance(X_LIMITS, list):
#        lim = X_LIMITS[r]
#        if lim is not None:
#            axes[r, 0].set_xlim(*lim)
#            axes[r, 1].set_xlim(*lim)
#    elif X_LIMITS is not None:
#        # rétro-compat : tuple global (xmin, xmax)
#        axes[r, 0].set_xlim(*X_LIMITS)
#        axes[r, 1].set_xlim(*X_LIMITS)
#
#plt.tight_layout(h_pad=1.2)
#plt.show()