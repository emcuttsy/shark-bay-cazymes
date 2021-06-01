import pandas as pd
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import seaborn as sns


def levels_heatmap(levels, df, palette='rocket', extend='max', figsize=(8,6), dpi=80, flipxy=False):
    df[df!=0] = df[df!=0] + 0.5
    if extend == 'both':
       colnum = len(levels)+1
    else:
       colnum = len(levels)

    colors = sns.color_palette(palette, colnum)
    cmap, norm = matplotlib.colors.from_levels_and_colors(levels, colors, extend=extend)

    if flipxy:
       df = df.T

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    im = ax.imshow(df, cmap=cmap, norm=norm, aspect='auto', interpolation='none')
    ax.set(xticks=range(df.shape[1]), yticks=range(df.shape[0]),
       xticklabels=df.columns, yticklabels=df.index)
    ax.tick_params(axis="x", rotation=90)
    fig.colorbar(im, ax=ax, spacing="proportional")
    return fig, ax