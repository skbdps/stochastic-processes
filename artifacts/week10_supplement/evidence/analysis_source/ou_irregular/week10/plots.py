"""Figures rendered only from bound Week 10 and original E4 tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..experiments.provenance import file_digest
from ..experiments.storage import atomic_json


def render(package):
    package = Path(package)
    output = package/"figures"
    if output.exists():
        raise FileExistsError("figure output must be new")
    output.mkdir()
    inputs = {}

    def read(name):
        inputs[name] = file_digest(package/name)
        return pd.read_csv(package/name)

    def save(fig, stem):
        fig.savefig(output/(stem+".png"), dpi=180, facecolor="white")
        fig.savefig(output/(stem+".pdf"), metadata={"CreationDate": None, "ModDate": None})
        plt.close(fig)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": .17})
    original = read("span/span_unadjusted.csv")
    adjusted = read("span/span_standardized.csv")
    horizon = read("analysis/matched_iid_comparisons.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.3))
    for offset, frame, label, color in zip((-.18, 0, .18), (original, adjusted, horizon),
            ("Original fixed N", "Span-standardized original", "New fixed T = 499.5"),
            ("#4C566A", "#D9822B", "#176B87")):
        rows = frame[(frame.estimator == "pfml") & (frame.parameter == "mu")].sort_values("ratio")
        for ax, metric, scale in ((axes[0], "rmse_ratio", 1), (axes[1], "coverage_difference", 100)):
            estimate = rows[metric].to_numpy()*scale
            lower = rows[metric+"_mc_lower"].to_numpy()*scale
            upper = rows[metric+"_mc_upper"].to_numpy()*scale
            ax.errorbar(np.arange(3)+offset, estimate, yerr=[estimate-lower, upper-estimate],
                        fmt="o", capsize=3, color=color, label=label)
            ax.set_xticks(range(3), ["1", "5", "20"])
            ax.set_xlabel("Arrival intensity ratio r")
    axes[0].axhline(1, color="black", ls=":", lw=1)
    axes[1].axhline(0, color="black", ls=":", lw=1)
    axes[0].set_ylabel("PFML mean RMSE ratio (MMPP / IID)")
    axes[1].set_ylabel("PFML mean coverage difference (percentage points)")
    axes[0].set_title("Precision contrast")
    axes[1].set_title("95% Wald coverage contrast")
    fig.suptitle("PFML mean contrasts across span controls", y=.97, fontsize=14)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(.5,.08), ncol=3)
    fig.text(.5,.025,"Pointwise 95% Monte Carlo intervals. Standardization conditions on selected support and estimated weights.",
             ha="center", fontsize=9)
    fig.subplots_adjust(bottom=.28, top=.84, wspace=.34)
    save(fig,"S1_span_controls")

    fixedn = read("analysis/original_e4_mechanism.csv")
    fixedt = read("analysis/horizon_mechanism.csv")
    order = [("regular",0),("iid",1),("mmpp",1),("iid",5),("mmpp",5),("iid",20),("mmpp",20)]
    labels = ["Regular", "IID\nr=1", "MMPP\nr=1", "IID\nr=5", "MMPP\nr=5", "IID\nr=20", "MMPP\nr=20"]
    fig, axes = plt.subplots(2,2,figsize=(12.5,8.2),sharey="row")
    for col, frame, title in ((0,fixedn,"Original fixed N = 1000"),(1,fixedt,"New fixed T = 499.5")):
        ordered = pd.DataFrame([frame[(frame.mechanism==m)&(frame.ratio==r)].iloc[0] for m,r in order])
        obs = "observed_sd" if col==0 else "observed_pfml_sd"
        pred = "fixed_rho_mean_sd" if col==0 else "predicted_fixed_rho_sd"
        se = "observed_rms_se" if col==0 else "observed_pfml_rms_se"
        cov = "observed_coverage" if col==0 else "observed_pfml_coverage"
        axes[0,col].plot(range(7),ordered[obs],"o-",color="#176B87",label="Empirical fitted-mean SD")
        axes[0,col].plot(range(7),ordered[pred],"s--",color="#D9822B",label="Fixed-coefficient covariance prediction")
        axes[0,col].plot(range(7),ordered[se],"^:",color="#6E6E6E",label="RMS reported standard error")
        axes[1,col].plot(range(7),ordered[cov]*100,"o-",color="#176B87",label="Observed coverage")
        axes[1,col].plot(range(7),ordered.normal_working_coverage_approximation*100,
                         "s--",color="#D9822B",label="Population-working Gaussian approximation")
        axes[1,col].axhline(95,color="black",ls=":",lw=1)
        axes[0,col].set_title(title)
        for row in range(2):
            axes[row,col].set_xticks(range(7),labels)
        axes[0,col].set_ylabel("Mean-estimation SD / reported SE")
        axes[1,col].set_ylabel("95% Wald coverage (%)")
    axes[0,0].legend(fontsize=8,loc="upper left")
    axes[1,0].legend(fontsize=8,loc="lower left")
    fig.suptitle("Higher-lag covariance predicts the PFML mean precision loss",fontsize=14)
    fig.text(.5,.02,"Covariance prediction is exact for a fixed-coefficient linear diagnostic; fitted-estimator and coverage comparisons are approximations.",
             ha="center",fontsize=9)
    fig.subplots_adjust(top=.91,bottom=.10,hspace=.32,wspace=.23)
    save(fig,"S2_mean_mechanism")

    parameters = read("analysis/parameter_summary.csv")
    paths = read("analysis/schedule_predictions.csv.gz")
    fig,axes=plt.subplots(2,3,figsize=(13,7.8))
    colors={"exact":"#176B87","pfml":"#D9822B","euler":"#7A5195"}
    for j,parameter in enumerate(("theta","mu","sigma")):
        for method in ("exact","pfml","euler"):
            for mechanism,style in (("mmpp","o-"),("iid","s--")):
                selected=parameters[(parameters.parameter==parameter)&(parameters.estimator==method)
                                    &(parameters.mechanism==mechanism)].sort_values("ratio")
                axes[0,j].plot(range(3),selected.coverage_operational*100,style,
                              color=colors[method],label=method.upper()+" "+mechanism.upper(),lw=1.2,ms=4)
        axes[0,j].set_xticks(range(3),["1","5","20"])
        axes[0,j].set_title({"theta": r"$\theta$", "mu": r"$\mu$", "sigma": r"$\sigma$"}[parameter]+" coverage")
        axes[0,j].set_xlabel("Arrival intensity ratio r")
        axes[0,j].axhline(95,color="black",ls=":",lw=1)
        axes[0,j].set_ylim((86,101) if parameter=="mu" else (-3,103))
        axes[0,j].set_ylabel("Operational coverage (%)")
    for j,ratio in enumerate((1,5,20)):
        selected=paths[paths.ratio==ratio]
        bins=np.linspace(selected.actual_n.min(),selected.actual_n.max(),24)
        for mechanism,color in (("mmpp","#176B87"),("iid","#D9822B")):
            axes[1,j].hist(selected[selected.mechanism==mechanism].actual_n,bins=bins,
                           histtype="step",linewidth=1.8,color=color,label=mechanism.upper())
        axes[1,j].axvline(1000,color="black",ls=":",lw=1)
        axes[1,j].set_title(f"Observation counts at r={ratio}")
        axes[1,j].set_xlabel("N including observations at 0 and T")
        axes[1,j].set_ylabel("Independent paths")
    axes[1,0].legend(fontsize=9)
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc="lower center",ncol=6,fontsize=9,bbox_to_anchor=(.5,.035))
    fig.suptitle("Fixed-horizon context: all parameters and random observation counts",fontsize=14)
    fig.text(.5,.015,"Each cell: 2,000 paths at T=499.5. Count reference: original N=1000; Poisson-window expected N=1001.",ha="center",fontsize=9)
    fig.subplots_adjust(top=.91,bottom=.15,hspace=.42,wspace=.29)
    save(fig,"S3_fixed_horizon_context")
    manifest={"inputs":inputs,"files":{p.name:file_digest(p) for p in sorted(output.iterdir())},
              "scope":"descriptive mechanism curves; pointwise MC intervals only in S1"}
    atomic_json(output/"figure_manifest.json",manifest,immutable=True)
    return manifest
