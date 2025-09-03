from AminoExtract import main


def test_translate_aa():
    args = [
        "--input",
        "tests/unit/data/ESIB_EQA_2024_SARS1_01_consensus.fa",
        "--features",
        "tests/unit/data/ESIB_EQA_2024_SARS1_01_consensus.gff",
        "--output",
        "tests/unit/data/ESIB_EQA_2024_SARS1_01_consensus.translated.fa",
        "-ft",
        "all",
        "-n",
        "ESIB_EQA_2024_SARS1_01",
        "--keep-gaps",
    ]

    a = main(provided_args=args)
