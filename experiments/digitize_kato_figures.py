"""Digitise the data markers of Kato (2022) Figs 4-7 from the PDF's embedded raster images.

The figures are 305 dpi JPEGs inside the PDF. This extracts them with `pdfimages`,
calibrates the y axis (and x axis for Figs 6-7) against the labelled major ticks, finds
the coloured marker discs, and prints / saves their data coordinates. It replaces reading
values off a 110 dpi render by eye.

    python experiments/digitize_kato_figures.py [--pdf journal.pone.0275113.pdf]
"""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

PALETTE = {  # RGB sampled from the markers of Fig 4
    "green": (1, 77, 1), "navy": (2, 1, 170), "purple": (130, 1, 127),
    "lilac": (170, 85, 169), "blue": (2, 1, 251),
}
FIGS = {"fig4": 1, "fig5": 2, "fig6": 3, "fig7": 4}  # index in `pdfimages` output for pages 7-11


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def runs(mask1d):
    """(start, end) of consecutive True runs."""
    idx = np.flatnonzero(mask1d)
    if idx.size == 0:
        return []
    cut = np.flatnonzero(np.diff(idx) > 1)
    starts = np.r_[idx[0], idx[cut + 1]]
    ends = np.r_[idx[cut], idx[-1]]
    return list(zip(starts, ends))


def find_axes(im):
    """Return (x of the y axis, y of the x axis, major tick rows, major tick columns)."""
    dark = im.max(axis=2) < 90
    h, w = dark.shape
    # y axis: column with the longest dark vertical run within the left 30%
    best = (0, 0)
    for x in range(0, int(w * 0.3)):
        r = max((e - s for s, e in runs(dark[:, x])), default=0)
        if r > best[0]:
            best = (r, x)
    ax = best[1]
    # x axis: row with the longest dark horizontal run in the lower 60%
    best = (0, 0)
    for y in range(int(h * 0.4), h):
        r = max((e - s for s, e in runs(dark[y, :])), default=0)
        if r > best[0]:
            best = (r, y)
    ay = best[1]
    # y ticks point inward: dark extent to the right of the axis line (~4 px wide)
    ext = np.zeros(h)
    for y in range(h):
        n = 0
        x = ax + 4
        while x < w and dark[y, x]:
            n += 1
            x += 1
        ext[y] = n
    ext[max(ay - 8, 0):] = 0     # the x-axis line itself is not a tick
    major_y = [(s + e) / 2 for s, e in runs(ext >= 0.7 * ext.max())]
    ext = np.zeros(w)
    for x in range(w):
        n = 0
        y = ay - 4
        while y > 0 and dark[y, x]:
            n += 1
            y -= 1
        ext[x] = n
    ext[:ax + 8] = 0             # nor is the y-axis line
    major_x = [(s + e) / 2 for s, e in runs(ext >= 0.7 * ext.max())] if ext.max() > 3 else []
    return ax, ay, major_y, major_x


def _merge_runs(mask1d, gap=6):
    idx = np.flatnonzero(mask1d)
    if idx.size == 0:
        return []
    cut = np.flatnonzero(np.diff(idx) > gap)
    return list(zip(np.r_[idx[0], idx[cut + 1]], np.r_[idx[cut], idx[-1]]))


def markers(im, rgb, tol=60, r_open=6):
    """Centroids of filled and hollow marker discs of one colour."""
    dist = np.sqrt(((im - np.array(rgb)) ** 2).sum(axis=2))
    m = dist < tol
    out = []
    for hollow in (False, True):
        mm = ndi.binary_fill_holes(ndi.binary_closing(m, disk(2))) if hollow else m
        op = ndi.binary_opening(mm, disk(r_open))
        lab, n = ndi.label(op)
        for k in range(1, n + 1):
            ys, xs = np.nonzero(lab == k)
            if ys.size < 60:
                continue
            cy, cx = ys.mean(), xs.mean()
            rr = disk(6)
            patch = m[int(cy) - 6:int(cy) + 7, int(cx) - 6:int(cx) + 7]
            solid = (patch & rr).sum() / rr.sum() if patch.shape == rr.shape else 0
            out.append(dict(x=cx, y=cy, solid=float(solid), area=int(ys.size)))
    # filled run first; a hollow candidate duplicating a solid one is dropped
    solid_pts = []
    for o in sorted((o for o in out if o["solid"] > 0.7), key=lambda o: -o["area"]):
        if not any(abs(o["x"] - q["x"]) < 8 and abs(o["y"] - q["y"]) < 8 for q in solid_pts):
            solid_pts.append(o)
    keep = list(solid_pts)
    hollow_pts = []
    for o in sorted((o for o in out if o["solid"] <= 0.7), key=lambda o: -o["area"]):
        if any(abs(o["x"] - q["x"]) < 8 and abs(o["y"] - q["y"]) < 8 for q in hollow_pts):
            continue
        hollow_pts.append(o)
    for o in hollow_pts:
        if o["solid"] <= 0.7 and not any(abs(o["x"] - s["x"]) < 8 and abs(o["y"] - s["y"]) < 8 for s in solid_pts):
            keep.append(dict(o, hollow=True))
    for o in keep:
        o.setdefault("hollow", False)
    return sorted(keep, key=lambda o: o["x"])


def digitise(pdf):
    work = Path(tempfile.mkdtemp())
    subprocess.run(["pdfimages", "-png", "-f", "7", "-l", "11", str(pdf), str(work / "img")], check=True)
    files = sorted(work.glob("img-*.png"))
    res = {}
    for name, idx in FIGS.items():
        im = np.array(Image.open(files[idx]).convert("RGB")).astype(float)
        ax, ay, my, mx = find_axes(im)
        # y: calibrate on the centres of the five tick labels 1.0 ... 0.2 (data lines occlude
        # the ticks in Figs 5-7). Labels sit ~2.4% of a tick spacing above their ticks (checked
        # against the unoccluded ticks of Fig 4), so shift by that.
        band = (im[:, max(ax - 105, 0):ax - 18, :].max(axis=2) < 110).sum(axis=1) > 0
        centres = [(s_ + e_) / 2 for s_, e_ in _merge_runs(band) if e_ - s_ > 15]
        fit_y = None
        if len(centres) == 5:
            d = np.diff(centres).mean()
            centres = [c + 0.0244 * d for c in centres]
            fit_y = np.polyfit(centres, [1.0, 0.8, 0.6, 0.4, 0.2], 1)
            ys = centres
        info = dict(axis_x=int(ax), axis_y=int(ay), y_ticks=[round(v, 1) for v in my], x_ticks=[round(v, 1) for v in mx])
        res[name] = dict(info=info, series={})
        if fit_y is None:
            continue
        resid = np.abs(np.polyval(fit_y, ys) - np.array([1.0, 0.8, 0.6, 0.4, 0.2])).max()
        info["y_fit_resid_gini"] = float(resid)
        info["gini_per_px"] = float(abs(fit_y[0]))
        for cname, rgb in PALETTE.items():
            pts = markers(im, rgb)
            res[name]["series"][cname] = [dict(px=round(p["x"], 1), py=round(p["y"], 1),
                                               gini=float(np.polyval(fit_y, p["y"])), hollow=p["hollow"],
                                               solid=round(p["solid"], 2)) for p in pts]
        res[name]["x_marker_px"] = sorted({round(p["px"]) for s in res[name]["series"].values() for p in s})
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default="journal.pone.0275113.pdf")
    ap.add_argument("--out", default="results/kato2022/digitised_figures.json")
    a = ap.parse_args()
    res = digitise(a.pdf)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res, indent=1))
    for fig, d in res.items():
        print(fig, d["info"])
        for c, pts in d["series"].items():
            print("  ", c, [("o" if p["hollow"] else "") + f"{p['gini']:.3f}" for p in pts])
