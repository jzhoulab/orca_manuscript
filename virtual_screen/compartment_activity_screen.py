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
os.makedirs("compartment_activity_outputs", exist_ok=True)


# load resources
import orca_predict

orca_predict.load_resources(models=[])
from orca_predict import *
from orca_models import HCTnoc

hctnoc = HCTnoc()
hctnoc.eval()
hctnoc.cuda()


sources = pd.read_csv("./swapdiff.screen.list.v3.csv")
targets = pd.read_csv("./swapdiff.screen.targetlist.v3.csv")

target_i = int(sys.argv[1])
chrm = targets["chr"][target_i]
chrlen = [len for c, len in hg38.get_chr_lens() if c == chrm][0]
chrlen_round = chrlen - chrlen % 4000

seq_start = np.fmax(targets["start"][target_i] - 16000000, 0)
seq_start = np.fmin(chrlen_round - 32000000, seq_start)


print(chrm, seq_start, seq_start + 32000000)
target_seq = hg38.get_encoding_from_coords(chrm, seq_start, seq_start + 32000000)

target_seq_start = targets["start"][target_i] - seq_start
target_seq_i = int(np.floor(target_seq_start / 128000))


blacklist = pyranges.read_bed("../resources/mcools.blacklist.4000.nleq10.bed")
b = pyranges.PyRanges(
    pd.DataFrame(
        dict(
            Chromosome=chrm,
            Start=seq_start + np.linspace(0, 32000000, 251)[:-1],
            End=seq_start + np.linspace(0, 32000000, 251)[:-1] + 128000,
            Index=np.arange(250),
        )
    )
)

mask = np.ones(250) == 1
rows = blacklist.join(b)
if len(rows) > 0:
    print(rows)
    rows_index = np.array(rows.Index)
    for i in np.unique(rows_index):
        if np.sum(rows_index == i) > 8:
            print(i, np.sum(rows_index == i))
            # if >25% is blacklisted
            mask[i] = False

if np.mean(mask) < 0.5:
    print("More than 50% of positions are blacklisted, aborting...")
    sys.exit()
else:
    print(np.mean(mask) * 100, "% unblacklisted regions.")


model = hctnoc

with torch.no_grad():
    sequence = torch.FloatTensor(target_seq[None, :, :]).clone()
    sequence_rc = torch.FloatTensor(target_seq[None, ::-1, ::-1].copy())
    seqT = torch.Tensor(sequence.float()).transpose(1, 2).cuda()
    seqT_rc = torch.Tensor(sequence_rc.float()).transpose(1, 2).cuda()
    encoding0 = (model.net0)(seqT)
    encoding0_rc = (model.net0)(seqT_rc)

    def eval_step(level, start, coarse_pred=None, plot=False):
        distenc = torch.log(
            torch.FloatTensor(model.normmats[level][None, None, :, :]).cuda()
        ).expand(sequence.shape[0], 1, 250, 250)
        if coarse_pred is not None:
            pred = model.denets[level].forward(
                encodings[level][:, :, int(start / level) : int(start / level) + 250],
                distenc,
                coarse_pred,
            )
        else:
            pred = model.denets[level].forward(
                encodings[level][:, :, int(start / level) : int(start / level) + 250], distenc
            )
        return pred

    start = 0
    encoding1, encoding2, encoding4, encoding8, encoding16, encoding32 = model.net(encoding0)
    encodings = {
        1: encoding1,
        2: encoding2,
        4: encoding4,
        8: encoding8,
        16: encoding16,
        32: encoding32,
    }
    pred32 = eval_step(32, start)
    encoding1, encoding2, encoding4, encoding8, encoding16, encoding32 = model.net(encoding0_rc)
    encodings = {
        1: encoding1,
        2: encoding2,
        4: encoding4,
        8: encoding8,
        16: encoding16,
        32: encoding32,
    }
    pred32_rc = torch.flip(eval_step(32, start), [2, 3])

    sources["score"] = np.nan
    for k in range(sources.shape[0]):

        sequence2 = sequence.clone()
        sequence2[0, target_seq_start : target_seq_start + 12800, :] = torch.FloatTensor(
            hg38.get_encoding_from_coords(sources["chr"][k], sources["start"][k], sources["end"][k])
        )
        sequence2_rc = sequence_rc.clone()
        sequence2_rc[
            :, 32000000 - target_seq_start - 12800 : 32000000 - target_seq_start, :
        ] = torch.FloatTensor(
            hg38.get_encoding_from_coords(
                sources["chr"][k], sources["start"][k], sources["end"][k]
            )[::-1, ::-1].copy()
        )

        padding = 112000
        ind_padding = int(padding / 4000)

        ind_s = int(np.fmax(np.floor((target_seq_start - padding * 2) / 4000), 0))
        ind_e = int(np.fmin(np.ceil(((target_seq_start + 12800 + padding * 2) / 4000)), 8000))

        ind_padding_offset_s = int(np.fmax((target_seq_start - padding) / 4000 - ind_s, 0))
        ind_padding_offset_e = int(np.fmax(ind_e - (target_seq_start + 12800 + padding) / 4000, 0))

        seqT = (
            torch.Tensor(sequence2[:, ind_s * 4000 : ind_e * 4000, :].float())
            .transpose(1, 2)
            .cuda()
        )
        seqT_rc = (
            torch.Tensor(
                sequence2_rc[:, 32000000 - ind_e * 4000 : 32000000 - ind_s * 4000, :].float()
            )
            .transpose(1, 2)
            .cuda()
        )
        encoding0_mut = encoding0.clone()
        encoding0_mut[:, :, ind_s + ind_padding_offset_s : ind_e - ind_padding_offset_e] = (
            model.net0
        )(seqT)[:, :, ind_padding_offset_s:-ind_padding_offset_e]

        encoding0_rc_mut = encoding0_rc.clone()
        encoding0_rc_mut[
            :, :, 8000 - ind_e + ind_padding_offset_e : 8000 - ind_s - ind_padding_offset_s
        ] = (model.net0)(seqT_rc)[:, :, ind_padding_offset_e:-ind_padding_offset_s]

        def eval_step(level, start, coarse_pred=None, plot=False):
            distenc = torch.log(
                torch.FloatTensor(model.normmats[level][None, None, :, :]).cuda()
            ).expand(sequence.shape[0], 1, 250, 250)
            if coarse_pred is not None:
                pred = model.denets[level].forward(
                    encodings[level][:, :, int(start / level) : int(start / level) + 250],
                    distenc,
                    coarse_pred,
                )
            else:
                pred = model.denets[level].forward(
                    encodings[level][:, :, int(start / level) : int(start / level) + 250], distenc
                )

            return pred

        start = 0
        encoding1, encoding2, encoding4, encoding8, encoding16, encoding32 = model.net(
            encoding0_mut
        )
        encodings = {
            1: encoding1,
            2: encoding2,
            4: encoding4,
            8: encoding8,
            16: encoding16,
            32: encoding32,
        }
        pred32_mut = eval_step(32, start)

        encoding1, encoding2, encoding4, encoding8, encoding16, encoding32 = model.net(
            encoding0_rc_mut
        )
        encodings = {
            1: encoding1,
            2: encoding2,
            4: encoding4,
            8: encoding8,
            16: encoding16,
            32: encoding32,
        }
        pred32_rc_mut = torch.flip(eval_step(32, start), [2, 3])

        sources.loc[k, "score"] = torch.mean(
            torch.abs(
                0.5 * pred32_mut[0, 0, target_seq_i, mask]
                - 0.5 * pred32[0, 0, target_seq_i, mask]
                + 0.5 * pred32_rc_mut[0, 0, target_seq_i, mask]
                - 0.5 * pred32_rc[0, 0, target_seq_i, mask]
            )
        ).item()
        print(
            k,
            torch.mean(
                torch.abs(
                    0.5 * pred32_mut[0, 0, target_seq_i, mask]
                    - 0.5 * pred32[0, 0, target_seq_i, mask]
                    + 0.5 * pred32_rc_mut[0, 0, target_seq_i, mask]
                    - 0.5 * pred32_rc[0, 0, target_seq_i, mask]
                )
            ).item(),
        )


sources.iloc[:, 1:].to_csv(
    "./compartment_activity_outputs/target." + sys.argv[1] + ".bedgraph",
    sep="\t",
    index=False,
    header=False,
)

