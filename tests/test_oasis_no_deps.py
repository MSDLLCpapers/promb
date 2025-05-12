from promb import init_db
import pytest

from promb.db import PrombDB


@pytest.fixture
def human_db() -> PrombDB:
    return init_db('human-oas')


def test_is_human(human_db: PrombDB):
    assert human_db.contains('QVQLVQSGV') == True
    assert human_db.contains('ELVISLIVE') == False


def test_compute_peptide_content(human_db: PrombDB):
    score = human_db.compute_peptide_content('QVQLQQSGAELARPGASVKMSCKAS')
    assert score == 8/17


def test_find_peptides_matching_template(human_db: PrombDB):
    peptides = human_db.find_peptides_matching_template('QVQL*QSGA')
    assert peptides == ['QVQLAQSGA', 'QVQLCQSGA', 'QVQLDQSGA', 'QVQLEQSGA', 'QVQLFQSGA', 'QVQLGQSGA', 'QVQLHQSGA', 'QVQLIQSGA', 'QVQLKQSGA', 'QVQLLQSGA', 'QVQLMQSGA', 'QVQLPQSGA', 'QVQLQQSGA', 'QVQLRQSGA', 'QVQLSQSGA', 'QVQLTQSGA', 'QVQLVQSGA', 'QVQLWQSGA']
    peptides = human_db.find_peptides_matching_template('*VQL*QSGA')
    assert peptides == ['DVQLAQSGA', 'DVQLEQSGA', 'DVQLGQSGA', 'DVQLLQSGA', 'DVQLMQSGA', 'DVQLVQSGA', 'EVQLAQSGA', 'EVQLDQSGA', 'EVQLEQSGA', 'EVQLFQSGA', 'EVQLGQSGA', 'EVQLIQSGA', 'EVQLKQSGA', 'EVQLLQSGA', 'EVQLMQSGA', 'EVQLPQSGA', 'EVQLQQSGA', 'EVQLRQSGA', 'EVQLSQSGA', 'EVQLTQSGA', 'EVQLVQSGA', 'EVQLWQSGA', 'HVQLAQSGA', 'HVQLEQSGA', 'HVQLFQSGA', 'HVQLGQSGA', 'HVQLIQSGA', 'HVQLLQSGA', 'HVQLMQSGA', 'HVQLPQSGA', 'HVQLQQSGA', 'HVQLVQSGA', 'QVQLAQSGA', 'QVQLCQSGA', 'QVQLDQSGA', 'QVQLEQSGA', 'QVQLFQSGA', 'QVQLGQSGA', 'QVQLHQSGA', 'QVQLIQSGA', 'QVQLKQSGA', 'QVQLLQSGA', 'QVQLMQSGA', 'QVQLPQSGA', 'QVQLQQSGA', 'QVQLRQSGA', 'QVQLSQSGA', 'QVQLTQSGA', 'QVQLVQSGA', 'QVQLWQSGA']
    peptides = human_db.find_peptides_matching_template('*VQL*QSG*')
    assert len(peptides) == 225
    peptides = human_db.find_peptides_matching_template('*V*L*QSG*')
    assert len(peptides) == 594


def test_find_point_mutant_peptides(human_db: PrombDB):
    peptides = human_db.find_point_mutant_peptides('QVQLQQSGA')
    assert len(peptides) == 71


