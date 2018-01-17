from fmbiopy.fmpaths import get_bowtie2_indices
from fmscripts.index_fasta import main

def test_index_fasta(dat):
    fasta = dat['tiny']['assemblies'][0]
    main(fasta)
    expected_outputs = [fasta.with_suffix('.fai', depth=0),
                        fasta.with_suffix('.dict')]
    expected_outputs += get_bowtie2_indices(fasta.with_suffix(''))
    for output in expected_outputs:
        assert output.exists()
