from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymsaviz import MsaViz, get_msa_testdata


def test_simple_all_run(msa_fasta_file: Path, tmp_path: Path):
    """Test simple all run (Only check if no error occurs)"""
    mv = MsaViz(msa_fasta_file)

    fig_outfile = tmp_path / "test.png"
    mv.savefig(fig_outfile)

    assert fig_outfile.exists()


def test_all_run_with_options(msa_fasta_file: Path, tmp_path: Path):
    """Test all run with options (Only check if no error occurs)"""
    mv = MsaViz(
        msa_fasta_file,
        color_scheme="Identity",
        wrap_length=50,
        show_label=False,
        show_seq_char=False,
        sort=True,
    )
    mv.set_highlight_pos([1, 5, (10, 13), 18])
    mv.set_highlight_pos_by_ident_thr(min_thr=80, max_thr=100)
    mv.add_markers([50, 51, 52, (60, 70), 80], marker="x", color="blue", size=6)
    mv.add_text_annotation(
        (100, 120), text="test", text_color="blue", text_size=10, range_color="blue"
    )

    fig_outfile = tmp_path / "test.png"
    mv.savefig(fig_outfile)

    assert fig_outfile.exists()


def test_basic_property():
    """Test basic property"""
    msa = MultipleSeqAlignment([])
    id_list = ["first", "second", "third", "fourth"]
    seq_list = [
        "CDNIPGFED",
        "ADNIPGFED",
        "BDNIPGFED",
        "DDNIPGFED",
    ]
    for id, seq in zip(id_list, seq_list):
        msa.append(SeqRecord(Seq(seq), id=id))

    mv = MsaViz(msa)
    assert mv.msa_count == 4
    assert mv.alignment_length == 9
    assert mv.id_list == id_list
    assert mv.seq_list == seq_list
    assert mv.wrap_num == 0
    assert mv.consensus_seq == "XDNIPGFED"


def test_set_custom_color_scheme(dummy_msa: MultipleSeqAlignment):
    """Test set_custom_color_scheme"""
    mv = MsaViz(dummy_msa)
    # Case1: Set correct custom color scheme
    custom_color_scheme = {"A": "red", "T": "blue", "G": "green", "C": "orange"}
    mv.set_custom_color_scheme(custom_color_scheme)
    assert mv.color_scheme == custom_color_scheme

    # Case2: Set invalid custom color scheme
    invalid_color_scheme = {"A": "invalid", "T": "blue", "G": "green", "C": "orange"}
    with pytest.raises(ValueError):
        mv.set_custom_color_scheme(invalid_color_scheme)


def test_set_custom_color_func(msa_fasta_file: Path, tmp_path: Path):
    """Test set_custom_color_func"""
    mv = MsaViz(msa_fasta_file)

    def custom_color_func(
        row_pos: int, col_pos: int, seq_char: str, msa: MultipleSeqAlignment
    ) -> str:
        if col_pos < 60 and seq_char != "-":
            return "salmon"
        if col_pos >= 60 and 1 <= row_pos <= 4:
            return "lime"
        return "white"

    mv.set_custom_color_func(custom_color_func)

    fig_outfile = tmp_path / "test.png"
    mv.savefig(fig_outfile)

    assert fig_outfile.exists()


def test_consensus_identity():
    """Test consensus identity calculation"""
    msa = MultipleSeqAlignment([])
    # Test MSA summary
    # 1: 'ABCDE'(All different char) => 'X' [20 %]
    # 2: 'GGGGG'(All 'G') => 'G' [100 %]
    # 3: '-----'(All gaps) => 'X' [0 %]
    # 4: '--V--'(one char & gaps) => 'V' [20 %]
    # 5: '-AAAC'('A' is most common) => 'A' [60 %]
    # 6: 'RRTTI'('R' & 'T' is most common) => 'X' [40 %]
    # 7: 'XXAX-'('X' is most common) => 'X' [60 %]
    seq_list = [
        "AG---RX",
        "BG--ARX",
        "CG-VATA",
        "DG--ATX",
        "EG--CI-",
    ]
    for seq in seq_list:
        msa.append(SeqRecord(Seq(seq)))

    # Test consensus seq & identity
    mv = MsaViz(msa)
    assert mv.consensus_seq == "XGXVAXX"
    consensus_ident_list = mv._get_consensus_identity_list()
    assert consensus_ident_list == [20, 100, 0, 20, 60, 40, 60]


def test_is_aa_msa():
    """Test `aa` or `nt` MSA check"""
    # Case1: AA MSA
    aa_msa = MultipleSeqAlignment(
        [
            SeqRecord(Seq("MFLTALLCRGRI")),
            SeqRecord(Seq("MFLT---TRGVI")),
        ]
    )
    assert MsaViz(aa_msa)._is_aa_msa() is True

    # Case2: NT MSA
    nt_msa = MultipleSeqAlignment(
        [
            SeqRecord(Seq("ATGC--TGCA")),
            SeqRecord(Seq("AAGCTCTGCA")),
        ]
    )
    assert MsaViz(nt_msa)._is_aa_msa() is False


def test_parse_positions(dummy_msa: MultipleSeqAlignment):
    """Test parse_positions"""
    mv = MsaViz(dummy_msa)
    # Case1: int value
    assert mv._parse_positions([1]) == [0]
    # Case2: int values
    assert mv._parse_positions([1, 5, 10, 20]) == [0, 4, 9, 19]
    # Case3: tuple range
    assert mv._parse_positions([(5, 9)]) == [4, 5, 6, 7, 8]
    # Case4: int values & tuple range
    assert mv._parse_positions([1, 5, (10, 13), 18]) == [0, 4, 9, 10, 11, 12, 17]


def test_get_msa_testdata():
    """Test get_msa_testdata"""
    assert get_msa_testdata().exists()
    assert get_msa_testdata("HIGD2A.fa").exists()
    with pytest.raises(ValueError):
        get_msa_testdata("invalid_name")
