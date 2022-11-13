from __future__ import annotations

import math
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Any
from urllib.request import urlopen

import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
from Bio.AlignIO import MultipleSeqAlignment
from matplotlib import colors
from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection
from matplotlib.colors import is_color_like
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from pymsaviz.config import COLOR_SCHEMES, AxesType


class MsaViz:
    """Multiple Sequence Alignment Visualization"""

    def __init__(
        self,
        msa: str | Path | MultipleSeqAlignment,
        format: str = "fasta",
        color_scheme: str | None = None,
        start: int = 1,
        end: int | None = None,
        wrap_length: int | None = None,
        wrap_space_size: float = 3.0,
        show_label: bool = True,
        show_seq_char: bool = True,
        show_grid: bool = False,
        show_count: bool = False,
        show_consensus: bool = False,
        consensus_color: str = "#1f77b4",
        consensus_size: float = 2.0,
    ):
        """
        Parameters
        ----------
        msa : str | Path | MultipleSeqAlignment
            MSA file, MSA file URL, MSA object
        format : str, optional
            Alignment file format (e.g. `fasta`, `phylip`, `clustal`, `emboss`, etc...)
        color_scheme : str | None, optional
            Color scheme. If None, `Zappo`(AA) or `Nucleotide`(NT) is set.
            [`Clustal`|`Zappo`|`Taylor`|`Flower`|`Blossom`|`Sunset`|`Ocean`|
            `Hydrophobicity`|`HelixPropensity`|`StrandPropensity`|`TurnPropensity`|
            `BuriedIndex`|`Nucleotide`|`Purine/Pyrimidine`|`None`]
        start : int, optional
            Start position of visualization (one-based coordinates)
        end : int | None, optional
            End position of visualization (one-based coordinates)
        wrap_length : int | None, optional
            Wrap sequence length. If None, no wrapping sequence.
        wrap_space_size: float, optional
            Space size between wrap MSA plot area
        show_label : bool, optional
            If True, show label
        show_seq_char : bool, optional
            If True, show sequence character
        show_grid : bool, optional
            If True, show grid
        show_count : bool, optional
            If True, show seq char count without gap on right side
        show_consensus : bool, optional
            If True, show consensus sequence
        consensus_color : str, optional
            Consensus identity bar color
        consensus_size : float, optional
            Consensus identity bar height size
        """
        # Load MSA
        if isinstance(msa, MultipleSeqAlignment):
            self._msa = msa
        elif isinstance(msa, str) and msa in ("http://", "https://"):
            content = urlopen(msa).read().decode("utf-8")
            self._msa = AlignIO.read(StringIO(content), format)
        else:
            self._msa: MultipleSeqAlignment = AlignIO.read(msa, format)
        self._msa_summary = SummaryInfo(self._msa)

        # Check & Set start, end position
        end = self.alignment_length if end is None else end
        if not 1 <= start <= end <= self.alignment_length:
            err_msg = f"start={start}, end={end} is invalid MSA range "
            err_msg += f"(1 <= start <= end <= {self.alignment_length})"
            raise ValueError(err_msg)
        self._start, self._end = start - 1, end
        self._length = self._end - self._start

        # Set user-specified plot configs
        if wrap_length in (0, None) or wrap_length > self._length:
            self._wrap_length = self._length
        else:
            self._wrap_length = wrap_length
        self._wrap_space_size = wrap_space_size
        self._show_seq_char = show_seq_char
        self._show_label = show_label
        self._show_grid = show_grid
        self._show_count = show_count
        self._show_consensus = show_consensus
        self._consensus_color = consensus_color
        self._consensus_size = consensus_size
        self._highlight_positions = None
        self._pos2marker_kws: dict[int, dict[str, Any]] = {}
        self._pos2text_kws: dict[int, dict[str, Any]] = {}

        # Plot configs (constant value)
        self._ticks_interval = 10
        self._x_unit_size_ratio = 0.7

        # Set color scheme
        if color_scheme is None:
            color_scheme = "Zappo" if self._is_aa_msa() else "Nucleotide"
        if color_scheme not in self.available_color_schemes():
            err_msg = f"color_scheme='{color_scheme}' is invalid.\n"
            err_msg += f"Available color scheme = {self.available_color_schemes()}"
            raise ValueError(err_msg)
        self._color_scheme = COLOR_SCHEMES[color_scheme]

    ############################################################
    # Property
    ############################################################

    @property
    def msa(self) -> MultipleSeqAlignment:
        """Multiple Sequence Alignment object (BioPython)"""
        return self._msa

    @property
    def msa_count(self) -> int:
        """MSA count"""
        return len(self._msa)

    @property
    def alignment_length(self) -> int:
        """Alignment length"""
        return self._msa.get_alignment_length()

    @property
    def id_list(self) -> list[str]:
        """MSA ID list"""
        return [rec.id for rec in self._msa]

    @property
    def seq_list(self) -> list[str]:
        """MSA sequence list"""
        return [str(rec.seq) for rec in self._msa]

    @property
    def wrap_num(self) -> int:
        """Wrap number"""
        if self._wrap_length is None:
            return 0
        else:
            return math.ceil(self._length / self._wrap_length) - 1

    @property
    def consensus_seq(self) -> str:
        """Consensus sequence"""
        return str(self._msa_summary.dumb_consensus(threshold=0))

    @property
    def color_scheme(self) -> dict[str, str]:
        """Color scheme"""
        return self._color_scheme

    ############################################################
    # Public Method
    ############################################################

    @staticmethod
    def available_color_schemes() -> list[str]:
        """Get available color schemes

        Returns
        -------
        color_scheme_names : list[str]
            Available color schemes
        """
        return list(COLOR_SCHEMES.keys())

    def set_custom_color_scheme(self, color_scheme: dict[str, str]) -> None:
        """Set user-defined custom color scheme (Overwrite color scheme setting)

        Parameters
        ----------
        color_scheme : dict[str, str]
            Custome color scheme dict (e.g. `{"A": "red", "R": "#F01505", ...}`)
        """
        if isinstance(color_scheme, dict):
            if not all(map(is_color_like, color_scheme.values())):
                err_msg = f"color_scheme={color_scheme} contains invalid color code."
                raise ValueError(err_msg)
            self._color_scheme = color_scheme
        else:
            raise ValueError(f"color_scheme='{color_scheme}' is not dict type.")

    def set_highlight_pos(self, positions: list[tuple[int, int] | int]) -> None:
        """Set user-defined highlight MSA positions

        Parameters
        ----------
        positions : list[tuple[int, int] | int]
            Highlight positions. int and tuple range mixture positions can be specified.
            (e.g. If `[1, 5, (10, 13), 18]` is set, `1, 5, 10, 11, 12, 13, 18`
            positions are highlighted)
        """
        self._highlight_positions = self._parse_positions(positions)

    def set_highlight_pos_by_ident_thr(
        self, min_thr: float = 0, max_thr: float = 100
    ) -> None:
        """Set highlight MSA positions by consensus identity threshold

        Parameters
        ----------
        min_thr : float, optional
            Min identity threshold for highlight position selection
        max_thr : float, optional
            Max identity threshold for highlight position selection
        """
        ident_list = self._get_consensus_identity_list()
        highlight_positions: list[int] = []
        for pos, ident in enumerate(ident_list):
            if min_thr <= ident <= max_thr:
                highlight_positions.append(pos)
        self._highlight_positions = highlight_positions

    def add_markers(
        self,
        positions: list[tuple[int, int] | int],
        marker: str = "v",
        color: str = "black",
        size: float = 6,
    ) -> None:
        """Add markers on specified positions

        Parameters
        ----------
        positions : list[tuple[int, int] | int]
            Marker positions. int and tuple range mixture positions can be specified.
            (e.g. If `[1, 5, (10, 13), 18]` is set, markers are plotted on
            `1, 5, 10, 11, 12, 13, 18` positions)
        marker : str, optional
            Marker type of matplotlib.
            See https://matplotlib.org/stable/api/markers_api.html for details.
        color : str, optional
            Marker color
        size : float, optional
            Marker size
        """
        for pos in self._parse_positions(positions):
            self._pos2marker_kws[pos] = dict(
                marker=marker,
                color=color,
                markersize=size,
                clip_on=False,
            )

    def add_text_annotation(
        self,
        range: tuple[int, int],
        text: str,
        text_color: str = "black",
        text_size: float = 10,
        range_color: str = "black",
    ) -> None:
        """Add text annotation in specified range

        Parameters
        ----------
        range : tuple[int, int]
            Annotation start-end range tuple
        text : str
            Annotation text
        text_color : str, optional
            Text color
        text_size : float, optional
            Text size
        range_color : str, optional
            Annotation range line color
        """
        # Add annotation text
        start, end = range[0] - 1, range[1]
        x = (start + end) / 2
        pos = int(x)
        self._pos2text_kws[pos] = dict(
            x=x,
            y=self.msa_count + 0.75,
            s=text,
            color=text_color,
            size=text_size,
            ha="center",
            va="bottom",
        )
        # Add annotation range line markers
        marker_size = 10 * (self._x_unit_size_ratio / 0.7)
        self.add_markers([range], marker="_", color=range_color, size=marker_size)

    def plotfig(self, dpi: int = 100) -> Figure:
        """Plot figure

        Parameters
        ----------
        dpi : int, optional
            Figure DPI

        Returns
        -------
        fig : Figure
            Figure
        """
        # Setup plot figure configs
        y_unit_size = 0.20
        x_unit_size = y_unit_size * self._x_unit_size_ratio

        ax_type2y_size = {
            AxesType.MSA: self.msa_count * y_unit_size,
            AxesType.SPACE: y_unit_size * 1.5,
            AxesType.CONSENSUS: y_unit_size * self._consensus_size,
            AxesType.WRAP_SPACE: y_unit_size * self._wrap_space_size,
        }

        plot_ax_types = []
        for wrap_idx in range(self.wrap_num + 1):
            plot_ax_types.append(AxesType.MSA)
            if self._show_consensus:
                plot_ax_types.append(AxesType.SPACE)
                plot_ax_types.append(AxesType.CONSENSUS)
            if wrap_idx != self.wrap_num:
                plot_ax_types.append(AxesType.WRAP_SPACE)

        y_size_list = [ax_type2y_size[t] for t in plot_ax_types]
        figsize = (self._wrap_length * x_unit_size, sum(y_size_list))
        fig: Figure = plt.figure(figsize=figsize, dpi=dpi, tight_layout=True)
        gs = GridSpec(nrows=len(plot_ax_types), ncols=1, height_ratios=y_size_list)
        gs.update(left=0, right=1, bottom=0, top=1, hspace=0, wspace=0)

        # Plot figure
        wrap_cnt = 0
        for idx, plot_ax_type in enumerate(plot_ax_types):
            ax: Axes = fig.add_subplot(gs[idx])
            if not isinstance(ax, Axes):
                raise TypeError("Error: Not matplotlib Axes class instance.")

            start = self._start + self._wrap_length * wrap_cnt
            end = self._start + self._wrap_length * (wrap_cnt + 1)
            end = self._end if end > self._end else end

            if plot_ax_type == AxesType.MSA:
                self._plot_msa(ax, start, end)
            elif plot_ax_type == AxesType.CONSENSUS:
                self._plot_consensus(ax, start, end)
            elif plot_ax_type == AxesType.SPACE:
                ax.axis("off")
            elif plot_ax_type == AxesType.WRAP_SPACE:
                ax.axis("off")
                wrap_cnt += 1
            else:
                raise NotImplementedError(f"plot_ax_type='{plot_ax_type}' is invalid.")

        return fig

    def savefig(
        self,
        savefile: str | Path,
        dpi: int = 100,
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : str | Path
            Save file
        dpi : int, optional
            DPI
        pad_inches : float, optional
            Padding inches
        """
        fig = self.plotfig(dpi=dpi)
        fig.savefig(
            fname=savefile,
            dpi=dpi,
            pad_inches=pad_inches,
        )

    ############################################################
    # Private Method
    ############################################################

    def _plot_msa(
        self, ax: Axes, start: int | None = None, end: int | None = None
    ) -> None:
        """Plot MSA

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to be plotted
        start : int | None, optional
            Start position. If None, `0` is set.
        end : int | None, optional
            End position. If None, `alignment_length` is set.
        """
        # Set xlim, ylim
        start = 0 if start is None else start
        end = self.alignment_length if end is None else end
        ax.set_xlim(start, start + self._wrap_length)
        ax.set_ylim(0, self.msa_count)

        # Set spines & tick params (Only show bottom ticklables)
        for pos in ("left", "right", "top", "bottom"):
            ax.spines[pos].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

        # Plot alignment position every 10 chars on xticks
        tick_ranges = range(start + 1, end + 1)
        xticklabels = list(filter(lambda n: n % self._ticks_interval == 0, tick_ranges))
        xticks = [n - 0.5 for n in xticklabels]
        ax.set_xticks(xticks, xticklabels, size=8)  # type: ignore

        plot_patches = []
        for cnt in range(self.msa_count):
            msa_seq = self.seq_list[cnt]
            y_lower = self.msa_count - (cnt + 1)
            y_center = y_lower + 0.5
            # Plot label text
            if self._show_label:
                msa_id = self.id_list[cnt]
                ax.text(start - 1, y_center, msa_id, ha="right", va="center", size=10)
            # Plot count text
            if self._show_count:
                scale = end - self._start - msa_seq[self._start : end].count("-")
                ax.text(end + 1, y_center, scale, ha="left", va="center", size=10)
            for x_left in range(start, end):
                # Add colored rectangle patch
                seq_char = msa_seq[x_left]
                seq_char_color = self.color_scheme.get(seq_char, "#FFFFFF")
                rect_prop = dict(xy=(x_left, y_lower), width=1, height=1)
                rect = Rectangle(**rect_prop, fc=seq_char_color)
                highlight_positions = self._highlight_positions
                if highlight_positions is None or x_left in highlight_positions:
                    plot_patches.append(rect)
                if self._show_grid:
                    grid = Rectangle(**rect_prop, ec="lightgrey", lw=0.5, fill=False)
                    plot_patches.append(grid)
                # Plot seq char text
                x_center = x_left + 0.5
                if self._show_seq_char:
                    # TODO: Better text y-position adjustment logic
                    y_adj = y_center if seq_char in ("-", "*") else y_center - 0.1
                    ax.text(
                        x_center, y_adj, seq_char, ha="center", va="center", size=10
                    )
                # Plot marker
                if cnt == 0 and x_left in self._pos2marker_kws:
                    marker_kws = self._pos2marker_kws[x_left]
                    ax.plot(x_center, y_center + 1, **marker_kws)
                # Plot text annotation
                if cnt == 0 and x_left in self._pos2text_kws:
                    text_kws = self._pos2text_kws[x_left]
                    ax.text(**text_kws)

        # Plot colored rectangle patch collection (Use collection for speedup)
        ax.add_collection(PatchCollection(plot_patches, match_original=True))

    def _plot_consensus(
        self, ax: Axes, start: int | None = None, end: int | None = None
    ) -> None:
        """Plot consensus seq char & identity bar

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to be plotted
        start : int | None, optional
            Start position. If None, `0` is set.
        end : int | None, optional
            End position. If None, `alignment_length` is set.
        """
        # Set xlim, ylim
        start = 0 if start is None else start
        end = self.alignment_length if end is None else end
        ax.set_xlim(start, start + self._wrap_length)
        ax.set_ylim(0, 100)  # 0 - 100 [%]

        # Plot label text
        if self._show_label and self._consensus_size != 0:
            ax.text(start - 1, 40, "Consensus", ha="right", va="center", size=10)

        # Set spines & tick params
        for pos in ("left", "right", "top", "bottom"):
            ax.spines[pos].set_visible(False)
        ax.tick_params(bottom=False, left=False, labelleft=False, pad=0)

        # Plot consensus seq chars on xticks
        xticks = list(map(lambda n: n + 0.5, range(start, end)))
        xticklabels = list(self.consensus_seq[start:end])
        ax.set_xticks(xticks, xticklabels, size=10)  # type: ignore

        # Plot consensus identity bar
        ident_list = self._get_consensus_identity_list(start, end)
        color_list = self._get_interpolate_colors(self._consensus_color, ident_list)
        ax.bar(xticks, ident_list, width=1, color=color_list, ec="white", lw=0.5)

    def _get_consensus_identity_list(
        self, start: int | None = None, end: int | None = None
    ) -> list[float]:
        """Get consensus identity list

        Parameters
        ----------
        start : int | None, optional
            Start position. If None, `0` is set.
        end : int | None, optional
            End position. If None, `alignment_length` is set.

        Returns
        -------
        consensus_identity_list : list[float]
            Consensus identity list (0 - 100 [%])
        """
        start = 0 if start is None else start
        end = self.alignment_length if end is None else end
        consensus_identity_list = []
        for idx, seq_char in enumerate(self.consensus_seq[start:end], start):
            column_chars = str(self.msa[:, idx])
            counter = Counter(filter(lambda c: c not in ("-", "*"), column_chars))
            count = counter.most_common()[0][1] if len(counter) != 0 else 0
            consensus_identity = (count / self.msa_count) * 100
            consensus_identity_list.append(consensus_identity)
        return consensus_identity_list

    def _get_interpolate_colors(
        self,
        color: str,
        values: list[float],
        vmin: float = 0,
        vmax: float = 100,
    ) -> list[str]:
        """Interpolate colors by size of values

        Parameters
        ----------
        color : str
            Base color for interpolation
        values : list[float]
            Values for interpolation
        vmin : float, optional
            Min value
        vmax : float, optional
            Max value

        Returns
        -------
        interpolated_colors : list[str]
            Interpolated colors based on values
        """
        cmap = colors.LinearSegmentedColormap.from_list("m", ("white", color))
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        return [colors.to_hex(cmap(norm(v))) for v in values]

    def _is_aa_msa(self) -> bool:
        """Check MSA is `aa` or `nt`

        If the ratio of `ATGCUN` char is less than 90%, return True.

        Returns
        -------
        check_result : bool
            Check result
        """
        nt_count, all_count = 0, 0
        for seq in self.seq_list:
            for seq_char in seq:
                if seq_char == "-":
                    continue
                all_count += 1
                if seq_char in "ATGCUN":
                    nt_count += 1
        return nt_count / all_count < 0.9

    def _parse_positions(self, positions: list[tuple[int, int] | int]) -> list[int]:
        """Parse int and tuple range mixture positions

        e.g. `[1, 5, (10, 13), 18]` means `1, 5, 10, 11, 12, 13, 18` positions

        Parameters
        ----------
        positions : list[tuple[int, int] | int]
            int and tuple range mixture positions (one-based coordinates)

        Returns
        -------
        result_positions : list[int]
            Parse result int positions (zero-based coordinates)
        """
        result_positions: list[int] = []
        for pos in positions:
            if isinstance(pos, (tuple, list)):
                result_positions.extend(list(range(pos[0] - 1, pos[1])))
            elif type(pos) == int:
                result_positions.append(pos - 1)
            else:
                raise ValueError(f"Invalid positions = {positions}")
        return sorted(set(result_positions))
