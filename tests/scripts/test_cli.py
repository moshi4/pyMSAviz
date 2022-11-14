import subprocess
from pathlib import Path


def test_cli_default_run(msa_fasta_file: Path, tmp_path: Path):
    """Test CLI with default option"""
    fig_outfile = tmp_path / "test.png"

    cmd = f"pymsaviz -i {msa_fasta_file} -o {fig_outfile}"
    subprocess.run(cmd, shell=True)

    assert fig_outfile.exists()


def test_cli_full_option_run(msa_fasta_file: Path, tmp_path: Path):
    """Test CLI with full option"""
    fig_outfile = tmp_path / "test.png"

    cmd = f"pymsaviz -i {msa_fasta_file} -o {fig_outfile} --format fasta "
    cmd += "--color_scheme Taylor --start 50 --end 250 --wrap_length 100 "
    cmd += "--wrap_space_size 3.0 --show_grid --show_count --show_consensus "
    cmd += "--consensus_color green --consensus_size 2.0 --sort --dpi 100"
    subprocess.run(cmd, shell=True)

    assert fig_outfile.exists()
