import argparse
import sys

import pandas as pd

args = argparse.ArgumentParser()

args.add_argument(
    "--primers",
    metavar="File",
    type=str,
    help="input file with primers as given by AmpliGone",
    required=True,
)

args.add_argument(
    "--coverages",
    metavar="File",
    type=str,
    help="Input file with coverages as given by TrueConsense",
    required=True,
)

args.add_argument("--key", metavar="String", type=str, help="Sample ID", required=True)

args.add_argument(
    "--output",
    metavar="File",
    type=str,
    help="Output file with average coverage per amplicon",
    required=True,
)

flags = args.parse_args()


def split_frames(df):
    left = ["LEFT", "PLUS", "POSITIVE"]
    right = ["RIGHT", "MINUS", "NEGATIVE"]

    leftdf = pd.DataFrame(columns=df.columns)
    rightdf = pd.DataFrame(columns=df.columns)

    for x in df.itertuples():
        if any(l in x[1] for l in left) is True:
            leftdf = leftdf.append(
                pd.DataFrame(
                    {"name": x.name, "start": x.start, "stop": x.stop}, index=[0]
                )
            )
        if any(r in x[1] for r in right) is True:
            rightdf = rightdf.append(
                pd.DataFrame(
                    {"name": x.name, "start": x.start, "stop": x.stop}, index=[0]
                )
            )

    leftdf.reset_index(inplace=True)
    rightdf.reset_index(inplace=True)

    leftdf = leftdf.drop(columns=["index"])
    rightdf = rightdf.drop(columns=["index"])

    return leftdf, rightdf


def remove_keyword(prname):
    keywords = ["LEFT", "RIGHT", "PLUS", "MINUS", "POSITIVE", "NEGATIVE"]
    sname = prname.split("_")
    for y, z in enumerate(sname):
        if z in keywords:
            sname.pop(y)
    return "_".join(sname)


def remove_alt_keyword(df):
    for x, y in enumerate(df["name"]):
        if "alt" in y:
            z = y.split("_")
            for a, b in enumerate(z):
                if "alt" in b:
                    z.pop(a)
            z = "_".join(z)
            df.at[x, "name"] = z
    return df


def index_to_remove_stops(one, indexone, two, indextwo):
    stop_one = one.get("stop")
    stop_two = two.get("stop")
    if stop_one > stop_two:
        return indextwo
    if stop_two > stop_one:
        return indexone
    return None


def index_to_remove_starts(one, indexone, two, indextwo):
    start_one = one.get("start")
    start_two = two.get("start")
    if start_one < start_two:
        return indextwo
    if start_two < start_one:
        return indexone
    return None


def remove_alt_primer_l(df):
    xx = df.to_dict(orient="records")
    to_rm = []
    lastindex = list(enumerate(xx))[-1][0]
    for a, x in enumerate(xx):
        if a != lastindex:
            if xx[a].get("name") == xx[a + 1].get("name"):
                rm_indx = index_to_remove_stops(xx[a], a, xx[a + 1], a + 1)
                if rm_indx is not None:
                    to_rm.append(rm_indx)
    filtereddf = df.drop(to_rm)
    return filtereddf


def remove_alt_primer_r(df):
    xx = df.to_dict(orient="records")
    to_rm = []
    lastindex = list(enumerate(xx))[-1][0]
    for a, x in enumerate(xx):
        if a != lastindex:
            if xx[a].get("name") == xx[a + 1].get("name"):
                rm_indx = index_to_remove_starts(xx[a], a, xx[a + 1], a + 1)
                if rm_indx is not None:
                    to_rm.append(rm_indx)
    filtereddf = df.drop(to_rm)
    return filtereddf


def Find_NonOverlap(df):
    dd = df.to_dict(orient="records")
    startingpoint = {}
    endingpoint = {}
    lastindex = list(enumerate(dd))[-1][0]
    firstindex = list(enumerate(dd))[0][0]
    for x, v in enumerate(dd):
        t_end = v.get("rightstart")
        if x != firstindex:
            s = dd[x - 1].get("rightstart")
        else:
            s = v.get("leftstop")
        if x != lastindex:
            end_override = dd[x + 1].get("leftstop")
        else:
            end_override = None
        if end_override is not None:
            if end_override in range(s, t_end):
                primerstart = s
                primerend = end_override
            else:
                primerstart = s
                primerend = t_end
        else:
            primerstart = s
            primerend = t_end
        startingpoint[primerstart] = v.get("name")
        endingpoint[primerend] = v.get("name")

    startdf = (
        pd.DataFrame.from_dict(startingpoint, orient="index")
        .reset_index()
        .rename(columns={0: "name", "index": "unique_start"})
    )
    enddf = (
        pd.DataFrame.from_dict(endingpoint, orient="index")
        .reset_index()
        .rename(columns={0: "name", "index": "unique_end"})
    )
    df = pd.merge(df, startdf, on="name", how="inner")
    df = pd.merge(df, enddf, on="name", how="inner")

    return df


def avg(lst):
    return round(sum(lst) / len(lst), 2)


def Average_cov(primers, covs):
    covd = covs.to_dict(orient="records")
    primd = primers.to_dict(orient="records")
    averages = {}

    for x, v in enumerate(primd):
        localcov = []

        prstart = v.get("unique_start")
        prend = v.get("unique_end")
        pr_range = list(range(prstart, prend))
        for i in pr_range:
            localcov.append(covd[i].get("cov"))
        averages[avg(localcov)] = v.get("name")

    avgdf = (
        pd.DataFrame.from_dict(averages, orient="index")
        .reset_index()
        .rename(columns={0: "name", "index": "avg_cov"})
    )
    primers = pd.merge(primers, avgdf, on="name", how="inner")

    return primers


if __name__ == "__main__":
    covs = pd.read_csv(
        flags.coverages, sep="\t", names=["position", "cov"], index_col="position"
    )

    try:
        prims = pd.read_csv(flags.primers)
    except Exception:
        print("Error reading primers file")
        with open(flags.output, "w") as f:
            f.write(
                f"""name,
{flags.key},
"""
            )
        sys.exit()

    if len(prims) <= 0:
        print("Primers file is empty, writing output and exiting...")
        with open(flags.output, "w") as f:
            f.write(
                f"""name,
{flags.key},
                """
            )
            sys.exit()

    lf, rf = split_frames(prims)

    lf["name"] = lf["name"].apply(remove_keyword)
    rf["name"] = rf["name"].apply(remove_keyword)

    lf = remove_alt_primer_l(remove_alt_keyword(lf))
    rf = remove_alt_primer_r(remove_alt_keyword(rf))

    non_overlapping_points = Find_NonOverlap(
        pd.merge(lf, rf, on="name", how="inner")
        .rename(
            columns={
                "start_x": "leftstart",
                "stop_x": "leftstop",
                "start_y": "rightstart",
                "stop_y": "rightstop",
            }
        )
        .drop_duplicates(subset="name")
    )

    with_average = Average_cov(non_overlapping_points, covs)

    with_average = with_average.drop(
        columns=[
            "leftstart",
            "leftstop",
            "rightstart",
            "rightstop",
            "unique_start",
            "unique_end",
        ]
    ).rename(columns={"avg_cov": flags.key})

    with_average = with_average.transpose()

    with_average.to_csv(flags.output, sep=",", index=True, header=False)
