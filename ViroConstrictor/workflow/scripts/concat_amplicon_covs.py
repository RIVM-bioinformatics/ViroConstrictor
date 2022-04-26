import argparse

import pandas as pd

args = argparse.ArgumentParser()

args.add_argument(
    "--input",
    type=str,
    nargs="+",
    default="1",
    metavar="File1 File2 File3",
    required=True,
)

args.add_argument(
    "--output",
    type=str,
    metavar="File",
    help="Output file name",
    required=True,
)

flags = args.parse_args()


def Make_frames(files):
    frames = []
    for i in files:
        df = pd.read_csv(i, sep=",")
        df.set_index("name", inplace=True)
        frames.append(df)
    return frames


def conc(frames):
    return pd.concat(frames, sort=True)


if __name__ == "__main__":

    Frame = conc(Make_frames(flags.input))

    Frame = Frame.fillna("NA")

    Frame.to_csv(flags.output, sep=",")
