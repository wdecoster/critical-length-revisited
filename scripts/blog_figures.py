#!/usr/bin/env python3
"""Figures for the blog post, from the workflow's summary tables (results/pilot/summary.tsv,
results/pilot/distributions.tsv, results/ladder/summary_steps1-2.tsv). Sniffles2 throughout, scored with GIAB's truvari
setup followed by truvari refine ("strict" in the workflow) unless --scoring says otherwise."""

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
STRATA = [("all", "all", "All SVs", BLUE, "o"), ("context", "tandem_repeats", "Tandem repeats", ORANGE, "s"),
          ("context", "segdups", "Segmental duplications", AQUA, "^")]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", required=True)
    parser.add_argument("--distributions", required=True)
    parser.add_argument("--ladder", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--scoring", default="strict", choices=["strict", "lenient"])
    return parser.parse_args()


def kb(x, _):
    return f"{x / 1000:g} kb"


def style(ax, xlabel, ylabel):
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(kb))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda *_: ""))
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED)


LABEL_OFFSET = {"All SVs": 7, "Tandem repeats": -7}  # the two upper lines run close together


def label_end(ax, x, y, text, color):
    ax.annotate(text, (x, y), xytext=(6, LABEL_OFFSET.get(text, 0)), textcoords="offset points", va="center",
                fontsize=9, color=INK)


def n50_figure(dist, scoring, out):
    d = dist[(dist.caller == "sniffles2") & (dist.scoring == scoring) & (dist["shape"] == "standard")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, cov in zip(axes, [15, 30]):
        for kind, stratum, name, color, marker in STRATA:
            x = d[(d.coverage == cov) & (d.kind == kind) & (d.stratum == stratum)].sort_values("achieved_n50")
            ax.plot(x.achieved_n50, x.f1, color=color, marker=marker, markersize=6, linewidth=2, label=name)
            label_end(ax, x.achieved_n50.iloc[-1], x.f1.iloc[-1], name, color)
        style(ax, "Read N50 of the library", "F1" if cov == 15 else "")
        ax.set_title(f"{cov}x", color=INK, loc="left", fontsize=11)
        ax.set_xlim(4000, 400000)
    axes[0].set_ylim(0.5, 1.0)
    axes[0].legend(frameon=False, loc="lower right", fontsize=9)
    fig.suptitle("SV discovery against the read N50 of the library (HG002, ONT R10, Sniffles2)", color=INK,
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fixed_figure(pilot, scoring, out):
    d = pilot[(pilot.read_set == "ulk") & (pilot.caller == "sniffles2") & (pilot.scoring == scoring)
              & (pilot.coverage == 15)]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for kind, stratum, name, color, marker in STRATA:
        x = d[(d.kind == kind) & (d.stratum == stratum)].groupby("length").f1.mean().sort_index()
        ax.plot(x.index, x.values, color=color, marker=marker, markersize=6, linewidth=2, label=name)
        label_end(ax, x.index[-1], x.values[-1], name, color)
    style(ax, "Read length (all reads cut to this length)", "F1")
    ax.set_xlim(800, 600000)
    ax.set_ylim(0.25, 1.0)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    ax.set_title("SV discovery against a fixed read length (15x)", color=INK, loc="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def then_now_figure(pilot, ladder, scoring, out):
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    series = [
        (ladder[(ladder.genome == "HG00733_2020") & (ladder.coverage == 40) & (ladder.setting == "default")],
         "2020 pipeline, HG00733 (40x simulated)", MUTED, "o", "--"),
        (ladder[(ladder.genome == "HG002_Q100") & (ladder.coverage == 20) & (ladder.setting == "default")],
         "2020 pipeline, HG002 Q100 (20x simulated)", ORANGE, "s", "--"),
    ]
    for x, name, color, marker, ls in series:
        x = x.sort_values("length")
        ax.plot(x.length, x.f1, color=color, marker=marker, markersize=6, linewidth=2, linestyle=ls, label=name)
    now = pilot[(pilot.read_set == "ulk") & (pilot.caller == "sniffles2") & (pilot.scoring == scoring)
                & (pilot.kind == "all") & (pilot.coverage == 15)].groupby("length").f1.mean().sort_index()
    ax.plot(now.index, now.values, color=BLUE, marker="^", markersize=6, linewidth=2,
            label="Now: real R10 reads cut to length, HG002 (15x)")
    style(ax, "Read length (all reads cut to this length)", "F1")
    ax.set_xlim(800, 700000)
    ax.set_ylim(0.2, 1.0)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    ax.set_title("Then and now", color=INK, loc="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": "white", "axes.facecolor": "white"})
    pilot = pd.read_table(args.pilot)
    pilot["length"] = pilot["length"].astype(int)
    dist = pd.read_table(args.distributions)
    ladder = pd.read_table(args.ladder)
    os.makedirs(args.outdir, exist_ok=True)
    n50_figure(dist, args.scoring, os.path.join(args.outdir, f"{args.prefix}n50.png"))
    fixed_figure(pilot, args.scoring, os.path.join(args.outdir, f"{args.prefix}fixed_length.png"))
    then_now_figure(pilot, ladder, args.scoring, os.path.join(args.outdir, f"{args.prefix}then_and_now.png"))


if __name__ == "__main__":
    main()
