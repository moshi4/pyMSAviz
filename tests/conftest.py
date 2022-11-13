from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def testdata_dir() -> Path:
    """Test data directory"""
    return Path(__file__).parent / "testdata"


@pytest.fixture
def msa_fasta_file(testdata_dir: Path) -> Path:
    """MSA fasta file"""
    return testdata_dir / "example.faa"


@pytest.fixture
def dummy_msa() -> MultipleSeqAlignment:
    """Dummy MSA object"""
    return MultipleSeqAlignment([SeqRecord(Seq("ATGC")), SeqRecord(Seq("ATGC"))])
