"""
This module provides functions and classes for parsing and validating command-line arguments,
processing primer data, calculating amplicon start and end positions, and computing mean coverage
for each amplicon.

The script requires the following input files:
- A BED file with primers as given by AmpliGone.
- A TSV file with coverages as given by TrueConsense.
- A sample ID.
- An output CSV file for average coverage per amplicon.

The module includes the following components:
- ReadDirection: An enum class to standardize read directions.
- standardize_direction: Function to standardize read direction strings.
- open_tsv_file: Function to open a TSV file and return its contents as a pandas DataFrame.
- split_primer_names: Function to split primer names in the DataFrame into separate columns.
- calculate_amplicon_start_end: Function to calculate the start and end positions of amplicons.
- create_amplicon_names_list: Function to create a list of unique amplicon names.
- write_output: Function to write the calculated amplicon sizes to a CSV file.
- main: Main function to execute the amplicon coverage calculation script.

Example usage:
--------------
To use this module, import it and call the main function with the appropriate arguments.

    from amplicon_covs import main

    if __name__ == "__main__":
        main()

Command-line usage:
-------------------
To run the script from the command line:
$ python amplicon_covs.py --primers primers.bed --coverages coverages.tsv --key sample1 --output output.csv
"""

from enum import Enum

import pandas as pd
from amplicon_arg_parser import parse_args


class ReadDirection(Enum):
    """
    An enum class to standardize read directions.

    Attributes
    ----------
    FORWARD : list of str
        List of strings representing forward read directions.

    REVERSE : list of str
        List of strings representing reverse read directions.

    Notes
    -----
    Add any new read directions to the respective value lists.
    """

    FORWARD = ["FW", "F", "1", "LEFT", "POSITIVE", "FORWARD", "PLUS"]
    REVERSE = ["RV", "R", "2", "RIGHT", "NEGATIVE", "REVERSE", "MINUS"]


def standardize_direction(direction: str) -> str:
    """
    Standardizes the given read direction string to a standardized direction name.

    Also needs the ReadDirection class to be defined.

    Parameters
    ----------
    direction : str
        The read direction string to be standardized.

    Returns
    -------
    str
        The standardized direction name, either 'FORWARD' or 'REVERSE'.

    Raises
    ------
    ValueError
        If the given direction does not match any recognized read direction.

    Examples
    --------
    >>> standardize_direction('F')
    'FORWARD'
    >>> standardize_direction('RIGHT')
    'REVERSE'
    """
    for read_direction in ReadDirection:  # ReadDirection.FORWARD, ReadDirection.REVERSE
        if direction.upper() in read_direction.value:
            return read_direction.name
    raise ValueError(
        f"Direction {direction} does not contain a recognized read direction."
        f"You can use {ReadDirection.FORWARD.value} for forward reads and {ReadDirection.REVERSE.value} for reverse reads."
    )


def open_tsv_file(filename: str, index_col: int | None = None) -> pd.DataFrame:
    """
    Opens a TSV file and returns its contents as a pandas DataFrame.

    This function checks if the file is properly tab-separated and contains no NaN values.

    When the input is parsed by amplicon_arg_parser, the file is guaranteed to exist and be readable. So there are no checks for that.

    Parameters
    ----------
    filename : str
        The path to the TSV file to be opened.

    index_col : int or None, optional
        Column to set as index (default is None).

    Returns
    -------
    pd.DataFrame
        The contents of the TSV file as a pandas DataFrame.

    Raises
    ------
    ValueError
        If the file contains NaN values.

    Examples
    --------
    >>> type(open_tsv_file('example.tsv'))
    <class 'pandas.core.frame.DataFrame'>
    """
    df = pd.read_csv(filename, sep="\t", header=None, index_col=index_col)
    if df.isnull().to_numpy().any():
        raise ValueError(f"File {filename} contains NaN values.")
    return df


def split_primer_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Splits the primer names in the DataFrame into separate columns.

    This function processes each row of the DataFrame by splitting the values in the 'split_name_list' column
    into 'name', 'count', 'alt', and 'direction' columns.

    Parameters
    ----------
    df : pd.DataFrame
        A pandas DataFrame containing primer data. The primer names should be in the third column (index 3).

    Returns
    -------
    pd.DataFrame
        The modified DataFrame with 'name', 'count', 'alt', and 'direction' columns populated.

    Raises
    ------
    ValueError
        If the 'split_name_list' column does not contain the expected number of elements (3 or 4).

    Examples
    --------
    >>> df = pd.DataFrame({3: ["NAME_NUMBER_DIRECTION", "NAME_NUMBER_ALT_DIRECTION"]})
    >>> split_primer_names(df)
           name   count   alt  direction
    0      NAME  NUMBER  None  DIRECTION
    1      NAME  NUMBER   ALT  DIRECTION
    """

    def _process_primer_row(row: pd.Series) -> pd.Series:
        """
        Processes a row of primer data by splitting the values from the 'split_name_list' column into their own column.

        Parameters
        ----------
        row : pd.Series
            A pandas Series representing a row of primer data. The 'split_name_list' column should contain
            a list of strings split from the primer name.

        Returns
        -------
        pd.Series
            The modified row with 'name', 'count', 'alt', and 'direction' columns populated.

        Raises
        ------
        ValueError
            If the 'split_name_list' column does not contain the expected number of elements (3 or 4).

        Examples
        --------
        >>> row = pd.Series({"split_name_list": ["NAME", "NUMBER", "DIRECTION"]})
        >>> process_primer_row(row)
        name         NAME
        count      NUMBER
        alt          None
        direction  DIRECTION
        dtype: object
        """
        split_names = row["split_name_list"]
        if len(split_names) == 4:
            row["name"], row["count"], row["alt"], row["direction"] = row[
                "split_name_list"
            ]
        elif len(split_names) == 3:
            row["name"], row["count"], row["alt"], row["direction"] = (
                split_names[0],
                split_names[1],
                None,
                split_names[2],
            )
        else:
            raise ValueError(
                f"Primer name {row[3]} does not contain the expected number of underscores. "
                "It should either be NAME_NUMBER_DIRECTION or NAME_NUMBER_ALT_DIRECTION."
            )
        return row

    df["split_name_list"] = df[3].str.split("_", expand=False)
    df = df.apply(_process_primer_row, axis=1)
    return df.drop(columns=["split_name_list"])


def calculate_amplicon_start_end(primers: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the start and end positions of amplicons based on primer data.

    This function processes the primer data to ensure that the start and end of an amplicon
    do not contain overlapping primers. It creates new rows with the minimum start and maximum end
    positions for each amplicon that has ALT primers, and then calculates the start and end positions
    for each amplicon.

    Parameters
    ----------
    primers : pd.DataFrame
        A pandas DataFrame containing primer data. The DataFrame should have columns for 'count',
        'alt', 'direction', and primer start and end positions.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame with columns 'amplicon_number', 'start', and 'end', representing the
        calculated start and end positions for each amplicon.

    Examples
    --------
    >>> primers = pd.DataFrame({
    ...     'count': [1, 1, 2, 2],
    ...     'alt': [None, 'ALT', None, 'ALT'],
    ...     'direction': ['FORWARD', 'FORWARD', 'REVERSE', 'REVERSE'],
    ...     1: [10, 20, 30, 40],
    ...     2: [50, 60, 70, 80],
    ...     3: ['name1', 'name2', 'name3', 'name4']
    ... })
    >>> calculate_amplicon_start_end(primers)
       amplicon_number  start  end
    0                1     10   60
    1                2     30   80
    """

    def _create_minmax_row(group: pd.DataFrame, direction: str) -> pd.Series:
        """
        Creates a new row with the minimum of the primer start and the maximum of the primer end.
        This way the start and end of an amplicon will never contain overlapping primers.

        Parameters
        ----------
        group : pd.DataFrame
            A pandas DataFrame containing the group of primer data.

        direction : str
            The direction of the read (e.g., 'FORWARD' or 'REVERSE').

        Returns
        -------
        pd.Series
            A pandas Series representing the new row with '_MINMAX' appended to the third column,
            the minimum value of the column with the primer start (column[1]), the maximum value of primer end (column[2]),
            and the direction.

        Examples
        --------
        >>> group = pd.DataFrame({
        ...     1: [10, 20, 30],
        ...     2: [40, 50, 60],
        ...     3: ['SC2_LEFT', 'SC2_LEFT_ALT1', 'SC2_LEFT_ALT2'],
        ...     4: ['FORWARD', 'FORWARD', 'FORWARD']
        ... })
        >>> _create_minmax_row(group, 'FORWARD')
        1          10
        2          60
        3    SC2_LEFT_MINMAX
        4      FORWARD
        dtype: object
        """
        new_row = group.iloc[0].copy()
        new_row[3] = new_row[3] + "_MINMAX"
        new_row[1] = group.loc[group[1].idxmin(), 1]
        new_row[2] = group.loc[group[2].idxmax(), 2]
        new_row[4] = direction
        return new_row

    df = pd.DataFrame(primers.loc[:, "count"].unique(), columns=["amplicon_number"])

    for amplicon_number in df["amplicon_number"]:
        amplicon_number_group = primers[primers["count"] == amplicon_number]
        if "ALT" in list(amplicon_number_group["alt"]):
            fw_group = amplicon_number_group[
                amplicon_number_group["direction"] == ReadDirection.FORWARD.name
            ]
            rv_group = amplicon_number_group[
                amplicon_number_group["direction"] == ReadDirection.REVERSE.name
            ]

            new_fw_row = _create_minmax_row(fw_group, ReadDirection.FORWARD.name)
            new_rv_row = _create_minmax_row(rv_group, ReadDirection.REVERSE.name)

            primers.loc[len(primers)] = new_fw_row
            primers.loc[len(primers)] = new_rv_row

            idxes_to_drop = fw_group.index.append(rv_group.index)
            primers = primers.drop(idxes_to_drop)

        # x = amplicon number
        # rows[1] = start primer
        # rows[2] = end primer
        # .values[0] = start position primer
        # .values[1] = end position primer
        if int(amplicon_number) == 1:  # make exceptions for first and last amplicon
            df.loc[df["amplicon_number"] == amplicon_number, "start"] = primers[
                primers["count"] == amplicon_number
            ][1].values[
                0
            ]  # no modification needed
        else:
            df.loc[df["amplicon_number"] == amplicon_number, "start"] = (
                primers[primers["count"] == str(int(amplicon_number) - 1)][2].values[1]
                + 1
            )  # add 1 to avoid overlapping primers
        if int(amplicon_number) == len(df["amplicon_number"]):
            df.loc[df["amplicon_number"] == amplicon_number, "end"] = primers[
                primers["count"] == amplicon_number
            ][2].values[1]
        else:
            df.loc[df["amplicon_number"] == amplicon_number, "end"] = primers[
                primers["count"] == str(int(amplicon_number) + 1)
            ][1].values[
                0
            ]  # no modification needed for end position
    return df


def create_amplicon_names_list(primers: pd.DataFrame) -> list[str]:
    """
    Create a list of unique amplicon names based on the primers DataFrame.

    This function generates a list of amplicon names by appending unique counts to the base amplicon name.

    Parameters
    ----------
    primers : pd.DataFrame
        A pandas DataFrame containing primer information. It must have columns "name" and "count".

    Returns
    -------
    list[str]
        A list of unique amplicon names.

    Examples
    --------
    >>> primers = pd.DataFrame({"name": ["amplicon1", "amplicon1"], "count": [1, 2]})
    >>> create_amplicon_names_list(primers)
    ['amplicon1_1', 'amplicon1_2']
    """
    amplicon_name = primers.loc[0, "name"]
    return [f"{amplicon_name}_{x}" for x in primers["count"].unique()]


def write_output(amplicon_sizes: pd.DataFrame, output_file: str) -> None:
    """
    Writes the calculated amplicon sizes to a CSV file.

    Parameters
    ----------
    amplicon_sizes : pd.DataFrame
        A pandas DataFrame containing amplicon size data.

    output_file : str
        The path to the output file to be written.
    """
    amplicon_sizes.to_csv(output_file, sep=",", index=True)


def main(args: list[str] | None = None) -> None:
    """
    Main function to execute the amplicon coverage calculation script.

    This function parses command-line arguments, reads input files, processes primer data,
    calculates amplicon start and end positions, computes mean coverage for each amplicon,
    and writes the results to an output file.

    Parameters
    ----------
    args : list[str], optional
        A list of command-line arguments to parse. If None, the arguments are taken from sys.argv.
        This parameter is used for testing purposes.

    Workflow
    --------
    1. Parse command-line arguments using parse_args.
    2. Open and read the primers file using open_tsv_file.
    3. Split primer names into separate columns using split_primer_names.
    4. Standardize primer directions using standardize_direction.
    5. Calculate amplicon start and end positions using calculate_amplicon_start_end.
    6. Open and read the coverages file using open_tsv_file.
    7. Calculate mean coverage for each amplicon.
    8. Create a list of unique amplicon names using create_amplicon_names_list.
    9. Write the final DataFrame to the output file using write_output.

    Examples
    --------
    To run the script from the command line:
    $ python amplicon_covs.py --primers primers.bed --coverages coverages.tsv --key sample1 --output output.csv
    """
    flags = parse_args(args)
    primers = open_tsv_file(flags.primers)
    primers = split_primer_names(primers)
    primers["direction"] = primers.loc[:, "direction"].apply(standardize_direction)
    amplicon_sizes = calculate_amplicon_start_end(primers)

    coverages = open_tsv_file(flags.coverages, index_col=0)

    def calculate_mean(input_array: pd.Series, coverages: pd.DataFrame) -> pd.Series:
        """
        Calculates the mean coverage for a given amplicon.

        This function computes the mean of all nucleotide coverages for the specified amplicon.
        It takes into account both complete and partial reads to measure primer efficiency.

        Parameters
        ----------
        input_array : pd.Series
            A pandas Series representing a row from the amplicon sizes DataFrame. It should contain
            the start and end positions of the amplicon.

        coverages : pd.DataFrame
            A pandas DataFrame containing coverage data. The index should represent nucleotide positions.

        Returns
        -------
        pd.Series
            A pandas Series with the mean coverage for the specified amplicon.

        Examples
        --------
        >>> input_array = pd.Series({"start": 10, "end": 20})
        >>> coverages = pd.DataFrame({"coverage": [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]})
        >>> calculate_mean(input_array, coverages)
        coverage    55.0
        dtype: float64
        """
        return coverages.iloc[
            int(input_array.iloc[1]) - 1 : int(input_array.iloc[2])
        ].mean()

    amplicon_sizes["coverage"] = amplicon_sizes.apply(
        lambda x: calculate_mean(x, coverages), axis=1
    )

    amplicon_sizes["amplicon_names"] = create_amplicon_names_list(primers)

    final_df = pd.DataFrame(
        [amplicon_sizes["coverage"].values],
        columns=[amplicon_sizes["amplicon_names"]],
        index=[flags.key],
    )

    write_output(final_df, flags.output)


if __name__ == "__main__":
    main()
