"""
references
    1. Devpractice Team. Python. Визуализация данных. Matplotlib. Seaborn. Mayavi. - devpractice.ru. 2020. - 412 с.: ил.
    2. A Guide to Formatting with f-strings in Python
class matplotlib.lines.Line2D(xdata, ydata, linewidth=None, linestyle=None, color=None, marker=None, markersize=None,
markeredgewidth=None, markeredgecolor=None, markerfacecolor=None, markerfacecoloralt='none', fillstyle=None,
antialiased=None, dash_capstyle=None, solid_capstyle=None, dash_joinstyle=None, solid_joinstyle=None, pickradius=5,
drawstyle=None, markevery=None, **kwargs)[source]
"""
import sys

import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np

from src.hex.hex import Hex


@dataclass
class Data2D:
    x: list
    y: list


@dataclass
class Line2D:
    width: float = 1
    style: str = '-'  # '-' or 'solid';'--' or 'dashed';'-.' or 'dashdot';':' or 'dotted'; 'none', 'None', ' ', or ''
    color: str = 'black'


@dataclass
class Marker:
    char: str = 'o'  # = '' -no marker
    size: float = 3
    edgewidth: float = None
    edgecolor: str = None
    facecolor = None
    facecoloralt = 'none'
    fillstyle = None


@dataclass
class Curve:
    xy: Data2D
    line: Line2D
    marker: Marker


@dataclass
class Text:
    x_label: str | None
    y_label: str | None
    title: str | None


@dataclass
class FigData:
    text: Text
    curve: list[Curve]


def plot_t_with_cross_sections(hex: Hex):
    # canvas
    cm = 1 / 2.54  # converter: centimeters in inches
    fig = plt.figure(
        figsize=(33 * cm, 18.5 * cm))  # set canvas with size (w, h) in cm  {33*cm, 18.5*cm for Asus Zenbook}
    fig.suptitle("QT-DIAGRAM")  # canvas title (common for all subplots)
    # следующие 2 стр.: получить положение верх.лев.точки окна и сдвинуть на "+100+10" (как-то это работает)
    thismanager = plt.get_current_fig_manager()
    thismanager.window.wm_geometry("+100+10")

    gs = fig.add_gridspec(ncols=2, nrows=2, hspace=0)  # set canvas grid: opt.arg. hspace=0

    # link grid with axs; axs ALWAYS have double indexing: 2x2: axs[0,0], axs[0,1], axs[1,0], axs[1,1]
    axs = gs.subplots(sharex=True)

    # set x-ticks manually (when automatically it not nice): index_min = 0, index_max = n - 1, step = 1
    # when axs = gs.subplots(sharex=True) this setting embrace ALL x-axes in canvas.
    plt.xticks(np.arange(0, len(hex.pack.hot.flow[0].t_k), 1))

    # axes labels and grids (x-axe is shared upon all subplots)
    axs[0, 0].set_ylabel("Temperature, K")
    # axs[0,0].set_xlabel("HEX cross-section index")
    axs[0, 0].grid(True)

    axs[1, 0].set_ylabel("dT[i,j] = Thot[i] - Tcold[j], K")
    axs[1, 0].set_xlabel("HEX cross-section index")
    axs[1, 0].grid(True)

    # axs[0, 1].set_ylabel("Total Q transferred, W")
    # axs[0,1].set_xlabel("HEX cross-section index")
    # axs[0, 1].grid(True)

    axs[0, 1].set_ylabel("Q per flow, W")
    axs[0, 1].set_xlabel("HEX cross-section index")
    axs[0, 1].grid(True)

    c_list_short = ["r", "b"]
    c_list_full = ["r", "b", "g", "c", "m", "y", "k", "w"]
    ls_list = ['-', '--', '-.', ':']
    mr_list = ['o', 's', '*', '+', 'x', '>', '<', 'd', '-']
    label_str = ['<~ hot ', '~> cold', '<~ total', '~> total']

    # -------------------------------------------- thot/cold(flow[i],n) -------------------------------------------------
    for i, pack in enumerate([hex.pack.hot, hex.pack.cold]):
        for j in range(len(pack.flow)):
            label = f'{label_str[i]} [{j}]'
            axs[0, 0].plot(pack.flow[j].t_k, alpha=0.7, label=label, c=c_list_short[i], ls=ls_list[0],
                           marker=mr_list[j], ms=4.0, mec=c_list_short[i], mfc='w', )
    # legend position loc= 0 - best; 1 - upper right; 2 - upper left. see [1] p. 37
    axs[0, 0].legend(loc=2, fontsize=8, title_fontsize=9,
                     title='         flow')
    # arrow with text: cold flow
    bbox_properties = dict(boxstyle='rarrow, pad=0.14', ec='b', fc='b', ls='-', lw=1)
    axs[0, 0].text(10.0, 148.0, 'cold flows', fontsize=8, fontweight='bold', color='w', bbox=bbox_properties)
    # arrow with text: hot flow
    bbox_properties = dict(boxstyle='larrow, pad=0.14', ec='r', fc='r', ls='-', lw=1)
    axs[0, 0].text(8.0, 315.0, 'hot flows', fontsize=8, fontweight='bold', color='w', bbox=bbox_properties)

    # -------------------------------------------- dt(n) ---------------------------------------------------------------
    k = 0
    for i in range(len(hex.pack.hot.flow)):
        for j in range(len(hex.pack.cold.flow)):
            label = f'[{i}] - [{j}]    ' \
                    f'{hex.dt_min_inter_flows[i][j]:^7.2f} ' \
                    f'{hex.dt_max_inter_flows[i][j]:^7.2f} ' \
                    f'{hex.dt_avr_inter_flows[i][j]:^7.2f}'
            axs[1, 0].plot(hex.dt_inter_flows[i][j], alpha=0.7, label=label, c=c_list_full[6], ls=ls_list[k],
                           marker=mr_list[k], ms=4.0, mec=c_list_full[6], mfc='w', )
            k += 1
    # legend position loc= 0 - best; 1 - upper right; 2 - upper left. see [1] p. 37
    axs[1, 0].legend(loc=2, fontsize=8, title_fontsize=9, title='      (i) - (j)              dt\n'
                                                                '                  ---------------------------\n'
                                                                '                   min      max     avr')
    # relative position: h, v; 0, 0 - almost center of the drawing; h<0 = move left; w>0
    # axs[1, 0].text(7, 130, 'dTmin = ???', fontsize=10, color='k')  # fontweight='b',

    zeroes = [1 for i in range(len(hex.pack.hot.dq_w))]
    # -------------------------------------------- Qsum(n) -------------------------------------------------------------
    # label = 'hot flows: emitted'
    # axs[0, 1].plot(hex.pack.hot.dq_w, alpha=0.7, label=label, c='r', ls=ls_list[0],
    #                marker=mr_list[0], ms=4.0, mec='r', mfc='w', )
    # label = 'cold flows: absorbed'
    # axs[0, 1].plot(hex.pack.cold.dq_w, alpha=0.7, label=label, c='b', ls=ls_list[0],
    #                marker=mr_list[0], ms=4.0, mec='b', mfc='w', )
    # axs[0, 1].legend(loc=2, title='        Total amount of heat ')
    # axs[0, 1].plot(zeroes, alpha=0.3, c='k', ls=ls_list[0])

    # -------------------------------------------- Qhot/cold(flow[i],n) -------------------------------------------------
    for i, pack in enumerate([hex.pack.hot, hex.pack.cold]):
        for j in range(len(pack.flow)):
            label = f'{label_str[i]} [{j}]'
            axs[0, 1].plot(pack.flow[j].dq_w, alpha=0.7, label=label, c=c_list_short[i], ls=ls_list[0],
                           marker=mr_list[j], ms=4.0, mec=c_list_short[i], mfc='w', )
    if len(hex.pack.hot.flow)>1:
        label = f'{label_str[2]}'
        axs[0, 1].plot(hex.pack.hot.dq_w, alpha=0.7, label=label, c='r', ls=ls_list[0],
                   marker=mr_list[2], ms=6.0, mec='r', mfc='r', )
    if len(hex.pack.cold.flow)>1:
        label = f'{label_str[3]}'
        axs[0, 1].plot(hex.pack.cold.dq_w, alpha=0.7, label=label, c='b', ls=ls_list[0],
                   marker=mr_list[2], ms=6.0, mec='b', mfc='b', )

    axs[0, 1].legend(loc=2, fontsize=8, title_fontsize=9,  title='            flow(i)')
    axs[0, 1].plot(zeroes, alpha=0.3, c='k', ls=ls_list[0])
    plt.show()


class PlotFig(ABC):
    @abstractmethod
    def __init__(self, fig_data: list[FigData]):
        self.n_figs = len(fig_data)
        self.fig = list()
        for data in fig_data:
            self.fig.append(FigData(text=data.text, curve=data.curve))
        # self.fig[i].text.x_label, self.fig[i].text.y_label, self.fig[i].text.title
        # self.fig[i].curve[j].x, self.fig[i].curve[j].y
        # self.fig[i].curve[j].line
        # self.fig[i].curve[j].marker

    @abstractmethod
    def plot_fig(self) -> None:
        pass


class PlotFig1Row1Col(PlotFig):

    def __init__(self, fig_data: list[FigData]):
        super().__init__(fig_data)

    def plot_fig(self) -> None:
        # Create just a figure and only one subplot
        fig, ax = plt.subplots()
        i = 0
        n_curves = len(self.fig[i].curve)
        for j in range(n_curves):
            ax.plot(self.fig[i].curve[j].xy.x, self.fig[i].curve[j].xy.y,
                    linewidth=self.fig[i].curve[j].line.width,
                    linestyle=self.fig[i].curve[j].line.style,
                    color=self.fig[i].curve[j].line.color,
                    marker=self.fig[i].curve[j].marker.char,
                    markersize=self.fig[i].curve[j].marker.size)

        ax.grid(True)
        ax.set_title(self.fig[i].text.title)
        ax.set_xlabel(self.fig[i].text.x_label)
        ax.set_ylabel(self.fig[i].text.y_label)

        plt.show()


class PlotFig2Rows1Col(PlotFig):
    """
    self.fig[0] - upper fig (row-1,col-1)
    self.fig[1] - lowe fig (row-2,col-1)
    common title saved in self.fig[0].text.title  => self.fig[1].text.title - has no title
    """

    def __init__(self, fig_data: list[FigData]):
        super().__init__(fig_data)

    def plot_fig(self) -> None:
        # Create figure and grid inside fig: grid of 1 column and 2 rows; no vertical gap: hspace=0
        fig = plt.figure(figsize=(7, 7))  # (w, h) in inch.
        gs = fig.add_gridspec(2, 1, hspace=0)
        # create a list of subplots: axs within fig and with shared x-axe (it is possible to share y: sharey=True
        axs = gs.subplots(sharex=True)  # axs is a list: axs[0]-upper fig (row-1,col-1); axs[1]-lower fig (row-2,col-1)

        # both figures share same title:
        fig.suptitle(self.fig[0].text.title)

        # axs[0]-upper fig (row-1,col-1)
        for i in range(len(self.fig)):
            axs[i].set_ylabel(self.fig[i].text.y_label)
            axs[i].set_xlabel(self.fig[i].text.x_label)
            axs[i].grid(True)
            for j in range(len(self.fig[i].curve)):
                axs[i].plot(self.fig[i].curve[j].xy.x, self.fig[i].curve[j].xy.y,
                            linewidth=self.fig[i].curve[j].line.width,
                            linestyle=self.fig[i].curve[j].line.style,
                            color=self.fig[i].curve[j].line.color,
                            marker=self.fig[i].curve[j].marker.char,
                            markersize=self.fig[i].curve[j].marker.size)
        plt.show()
