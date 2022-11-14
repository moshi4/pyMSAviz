from __future__ import annotations

import argparse
from pathlib import Path

from pymsaviz import MsaViz, __version__


def main():
    """Main function called from CLI"""
    args = get_args()
    run(**args.__dict__)


def run(
    infile: str | Path,
    outfile: str | Path,
    format: str = "fasta",
    color_scheme: str = "Zappo",
    start: int = 1,
    end: int | None = None,
    wrap_length: int | None = None,
    wrap_space_size: float = 3.0,
    show_grid: bool = False,
    show_count: bool = False,
    show_consensus: bool = False,
    consensus_color: str = "#1f77b4",
    consensus_size: float = 2.0,
    sort: bool = False,
    dpi: int = 300,
):
    """Run MSA visualization"""
    mv = MsaViz(
        msa=infile,
        format=format,
        start=start,
        end=end,
        wrap_length=wrap_length,
        wrap_space_size=wrap_space_size,
        color_scheme=color_scheme,
        show_grid=show_grid,
        show_count=show_count,
        show_consensus=show_consensus,
        consensus_color=consensus_color,
        consensus_size=consensus_size,
        sort=sort,
    )
    mv.savefig(outfile, dpi=dpi)


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """
    description = "MSA(Multiple Sequence Alignment) visualization CLI tool"
    parser = argparse.ArgumentParser(
        description=description,
        add_help=False,
        epilog=f"Available Color Schemes:\n{MsaViz.available_color_schemes()}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        type=Path,
        help="Input MSA file",
        metavar="I",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=Path,
        help="Output MSA visualization file (*.png|*.jpg|*.svg|*.pdf)",
        required=True,
        metavar="O",
    )
    default_msa_format = "fasta"
    parser.add_argument(
        "--format",
        type=str,
        help=f"MSA file format (Default: '{default_msa_format}')",
        default=default_msa_format,
        metavar="",
    )
    default_color_scheme = "Zappo"
    parser.add_argument(
        "--color_scheme",
        type=str,
        help=f"Color scheme (Default: '{default_color_scheme}')",
        default=default_color_scheme,
        choices=MsaViz.available_color_schemes(),
        metavar="",
    )
    default_start = 1
    parser.add_argument(
        "--start",
        type=int,
        help=f"Start position of MSA visualization (Default: {default_start})",
        default=default_start,
        metavar="",
    )
    default_end = None
    parser.add_argument(
        "--end",
        type=int,
        help="End position of MSA visualization (Default: 'MSA Length')",
        default=default_end,
        metavar="",
    )
    default_wrap_length = None
    parser.add_argument(
        "--wrap_length",
        type=int,
        help=f"Wrap length (Default: {default_wrap_length})",
        default=default_wrap_length,
        metavar="",
    )
    default_wrap_space_size = 3.0
    parser.add_argument(
        "--wrap_space_size",
        type=float,
        help="Space size between wrap MSA plot area "
        f"(Default: {default_wrap_space_size})",
        default=default_wrap_space_size,
        metavar="",
    )
    parser.add_argument(
        "--show_grid",
        help="Show grid (Default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "--show_count",
        help="Show seq char count without gap on right side (Default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "--show_consensus",
        help="Show consensus sequence (Default: OFF)",
        action="store_true",
    )
    default_consensus_color = "#1f77b4"
    parser.add_argument(
        "--consensus_color",
        type=str,
        help=f"Consensus identity bar color (Default: '{default_consensus_color}')",
        default=default_consensus_color,
        metavar="",
    )
    default_consensus_size = 2.0
    parser.add_argument(
        "--consensus_size",
        type=float,
        help=f"Consensus identity bar height size (Default: {default_consensus_size})",
        default=default_consensus_size,
        metavar="",
    )
    parser.add_argument(
        "--sort",
        help="Sort MSA order by NJ tree constructed from MSA distance matrix "
        "(Default: OFF)",
        action="store_true",
    )
    default_dpi = 300
    parser.add_argument(
        "--dpi",
        type=int,
        help=f"Figure DPI (Default: {default_dpi})",
        default=default_dpi,
        metavar="",
    )
    parser.add_argument(
        "-v",
        "--version",
        version=f"v{__version__}",
        help="Print version information",
        action="version",
    )
    parser.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
