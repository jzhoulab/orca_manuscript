import os
import sys
import numpy as np
import pandas as pd
import torch
import pyranges

if os.getenv("ORCA_PATH"):
    sys.path.append(os.getenv("ORCA_PATH"))
else:
    # change this to the right path if you use a different path
    # or specify the ORCA_PATH environmental variable
    sys.path.append("../../orca/")
sys.path.append("../")
os.makedirs("local_interaction_outputs", exist_ok=True)

# Change to true to include HCTnoc_1m model
process_hctnoc = False

# load resources
import orca_predict

orca_predict.load_resources(models=["1M"])
from orca_predict import *

if process_hctnoc:
    from orca_models_extra import Hctnoc_1M

    hct_1m = Hctnoc_1M()
    hct_1m.eval()
    hct_1m.cuda()


def pred_1m(seq, model):
    pred = model(seq.transpose(1, 2))
    return pred


chrm, start = sys.argv[1], int(sys.argv[2])

sequence = hg38.get_encoding_from_coords(chrm, start, start + 1000000)[None, :, :]
sequence_containsn = np.reshape(sequence == 0.25, (100000, 40)).any(axis=1)

blacklist = pyranges.read_bed("../resources/mcools.blacklist.4000.nleq10.bed")
b = pyranges.PyRanges(
    pd.DataFrame(
        dict(
            Chromosome=chrm,
            Start=start + np.linspace(0, 1000000, 251)[:-1],
            End=start + np.linspace(0, 1000000, 251)[:-1] + 4000,
            Index=np.arange(250),
        )
    )
)

mask = np.ones(250) == 1
rows = blacklist.join(b)
if len(rows) > 0:
    rows_index = np.array(rows.Index)
    mask[rows_index] = False

if np.mean(mask) < 0.5:
    print("More than 50% of positions are blacklisted, aborting...")
    sys.exit()


alldiffsh1esc = []
alldiffshff = []
if process_hctnoc:
    alldiffshctnoc = []

for i in range(3):
    basedesign = np.arange(10000, 90000).reshape(20, 4000)
    for i in range(basedesign.shape[0]):
        basedesign[i, :] = basedesign[i, np.random.permutation(4000)]
    with torch.no_grad():
        pred = (
            pred_1m(torch.FloatTensor(sequence).cuda(), h1esc_1m).squeeze().detach().cpu().numpy()
        )
        predhff = (
            pred_1m(torch.FloatTensor(sequence).cuda(), hff_1m).squeeze().detach().cpu().numpy()
        )
        if process_hctnoc:
            predhct = (
                pred_1m(torch.FloatTensor(sequence).cuda(), hct_1m).squeeze().detach().cpu().numpy()
            )

        diffsh1esc = np.ones(100000)
        diffshff = np.ones(100000)
        if process_hctnoc:
            diffshctnoc = np.ones(100000)
        for j in range(basedesign.shape[1]):
            sequence_mut = sequence.copy()
            for i in basedesign[:, j]:
                sequence_mut[0, i * 10 : (i + 1) * 10, :] = sequence[
                    0, np.random.randint(0, sequence.shape[1], 10), :
                ]
            pred_mut = (
                pred_1m(torch.FloatTensor(sequence_mut).cuda(), h1esc_1m)
                .squeeze()
                .detach()
                .cpu()
                .numpy()
            )
            predhff_mut = (
                pred_1m(torch.FloatTensor(sequence_mut).cuda(), hff_1m)
                .squeeze()
                .detach()
                .cpu()
                .numpy()
            )
            if process_hctnoc:
                predhct_mut = (
                    pred_1m(torch.FloatTensor(sequence_mut).cuda(), hct_1m)
                    .squeeze()
                    .detach()
                    .cpu()
                    .numpy()
                )
            diffsh1esc[basedesign[:, j]] = np.mean(
                np.abs(
                    pred_mut[np.floor(basedesign[:, j] / 400).astype(int), :][:, mask]
                    - pred[np.floor(basedesign[:, j] / 400).astype(int), :][:, mask]
                ),
                axis=1,
            )
            diffshff[basedesign[:, j]] = np.mean(
                np.abs(
                    predhff_mut[np.floor(basedesign[:, j] / 400).astype(int), :][:, mask]
                    - predhff[np.floor(basedesign[:, j] / 400).astype(int), :][:, mask]
                ),
                axis=1,
            )
            if process_hctnoc:
                diffshctnoc[basedesign[:, j]] = np.mean(
                    np.abs(
                        predhct_mut[np.floor(basedesign[:, j] / 400).astype(int), :][:, mask]
                        - predhct[np.floor(basedesign[:, j] / 400).astype(int), :][:, mask]
                    ),
                    axis=1,
                )
        diffsh1esc[sequence_containsn] = 0
        diffshff[sequence_containsn] = 0
        if process_hctnoc:
            diffshctnoc[sequence_containsn] = 0
        alldiffsh1esc.append(diffsh1esc)
        alldiffshff.append(diffshff)
        if process_hctnoc:
            alldiffshctnoc.append(diffshctnoc)


pd.DataFrame(
    {
        "chr": np.repeat(chrm, 100000),
        "start": start + np.arange(100000) * 10,
        "end": start + np.arange(1, 100001) * 10,
        "score": np.fmin(np.fmin(alldiffsh1esc[0], alldiffsh1esc[1]), alldiffsh1esc[2]),
    }
).iloc[10000:90000, :].to_csv(
    "./local_interaction_outputs/local_interaction.v2."
    + chrm
    + "."
    + str(start)
    + ".h1esc_1m.bedgraph",
    sep="\t",
    header=None,
    index=False,
)

pd.DataFrame(
    {
        "chr": np.repeat(chrm, 100000),
        "start": start + np.arange(100000) * 10,
        "end": start + np.arange(1, 100001) * 10,
        "score": np.fmin(np.fmin(alldiffshff[0], alldiffshff[1]), alldiffshff[2]),
    }
).iloc[10000:90000, :].to_csv(
    "./local_interaction_outputs/local_interaction.v2."
    + chrm
    + "."
    + str(start)
    + ".hff_1m.bedgraph",
    sep="\t",
    header=None,
    index=False,
)

if process_hctnoc:
    pd.DataFrame(
        {
            "chr": np.repeat(chrm, 100000),
            "start": start + np.arange(100000) * 10,
            "end": start + np.arange(1, 100001) * 10,
            "score": np.fmin(np.fmin(alldiffshctnoc[0], alldiffshctnoc[1]), alldiffshctnoc[2]),
        }
    ).iloc[10000:90000, :].to_csv(
        "./local_interaction_outputs/local_interaction."
        + chrm
        + "."
        + str(start)
        + ".hctnoc_1m.bedgraph",
        sep="\t",
        header=None,
        index=False,
    )

