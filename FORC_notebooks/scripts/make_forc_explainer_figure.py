"""Build the four-panel FORC explainer figure for FORC_intro.md.

The figure walks through what a FORC measurement is and how the distribution
comes out of it, using the Baraboo hematite vein material example distributed
with RockmagPy rather than a schematic:

    (a) one first-order reversal curve within the major hysteresis loop,
        with the reversal field Ha and one measured point M(Ha, Hb) marked;
    (b) the full family of 269 reversal curves that fills the loop;
    (c) the same measurements as a grid in (Hb, Ha) coordinates, coloured by
        magnetization, showing that the data occupy the half-plane Hb >= Ha;
    (d) the FORC distribution, the mixed second derivative of that surface,
        on the same axes, with the rotated Bc and Bu axes drawn on top.

Run from anywhere:

    python FORC_notebooks/scripts/make_forc_explainer_figure.py

Writes ``book/images/forc_explainer.png`` and ``.pdf``.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.patches import FancyArrowPatch

import pmagpy.forc as forc

REPO = Path(__file__).resolve().parents[2]
RAW_FILE = REPO / "example_data" / "FORC" / "baraboo_vein_material_hematite.txt"
OUT_STEM = REPO / "book" / "images" / "forc_explainer"

B_STEP = 0.005          # regridding step, T
SMOOTH_STRENGTH = 1.0   # multiplier on the automated LOESS spans
HIGHLIGHT_HA = -0.30    # reversal field of the curve featured in panel (a), T

MOMENT_SCALE = 1e6      # A m^2 -> micro A m^2
MOMENT_LABEL = r"$M$ ($\mu$A m$^2$)"

CURVE_COLOR = "#1f4aa8"
ACCENT = "#c4172a"


def load_baraboo():
    """Read the example file and compute the gridded surface and distribution.

    Returns:
        tuple: the list of reversal-curve Segments, the reversal fields Ha,
        the applied fields Hb, the gridded magnetization M(Ha, Hb), and the
        FORC distribution rho on the same grid. Fields are in tesla.
    """
    segments, _ = forc.phase1_prepare_segments_dual(
        str(RAW_FILE), export_magic=False, verbose=False
    )
    curves = [s for s in segments if s.kind == "forc"]

    out = forc.process_forc(
        mode="i", path=str(RAW_FILE), do_regrid=True, B_step=B_STEP,
        smooth_strength=SMOOTH_STRENGTH,
        plot_hyst=False, plot_rho=False, verbose=False,
    )
    return curves, out["Ha_vals_used"], out["Hb_vals_used"], out["M_grid_used"], out["rho"]


def major_loop_from_curves(curves):
    """Recover both branches of the major loop from the measurements alone.

    The lowest reversal curve begins at the most negative field reached and so
    traces the ascending branch. The descending branch needs no reconstruction
    either: every FORC begins on it, at the point (Ha, M(Ha, Ha)), so the locus
    of curve start points is the measured descending branch.

    Args:
        curves: Reversal-curve Segments, in file order.

    Returns:
        tuple: field and moment of the ascending branch, then of the
        descending branch, each sorted by increasing field.
    """
    # The ascending branch is the outer envelope of the family: every FORC
    # lies inside the major loop, so at any applied field the lowest moment
    # reached across all curves is the branch itself. Taking the envelope
    # rather than the single lowest curve matters because each curve spans a
    # fixed field window, and the lowest one stops well short of the maximum
    # field that other curves reach.
    H_all = np.unique(np.concatenate([np.asarray(s.H, float) for s in curves]))
    stack = np.full((len(curves), H_all.size), np.nan)
    for k, seg in enumerate(curves):
        H, M = np.asarray(seg.H, float), np.asarray(seg.M, float)
        inside = (H_all >= H.min()) & (H_all <= H.max())
        stack[k, inside] = np.interp(H_all[inside], H, M)
    with np.errstate(invalid="ignore"):
        M_env = np.nanmin(stack, axis=0)
    ok = np.isfinite(M_env)
    H_up, M_up = H_all[ok], M_env[ok]

    starts = sorted(((float(s.Ha), float(s.M[0])) for s in curves
                     if s.Ha is not None and len(s.M)), key=lambda t: t[0])
    H_dn = np.array([h for h, _ in starts])
    M_dn = np.array([m for _, m in starts])

    # No curve begins above the highest reversal field, so the descending
    # branch is unmeasured there and the loop is left open at the top. A
    # hysteresis loop is inversion symmetric, M(-H) = -M(H), so that segment
    # can be taken from the ascending branch. It is used only above the
    # highest reversal field, where there is nothing to contradict it; at the
    # join the two agree to 0.08% of the peak moment for this specimen.
    gap = H_up > H_dn[-1]
    if gap.any():
        H_fill = H_up[gap]
        M_fill = -np.interp(-H_fill, H_up, M_up)
        H_dn = np.concatenate([H_dn, H_fill])
        M_dn = np.concatenate([M_dn, M_fill])
    return H_up, M_up, H_dn, M_dn


def _zero_axes(ax):
    """Mark H = 0 and M = 0, which orient the eye on a hysteresis plot."""
    ax.axhline(0.0, color="0.78", lw=0.7, zorder=-5)
    ax.axvline(0.0, color="0.78", lw=0.7, zorder=-5)


def panel_a(ax, curves):
    """Draw one reversal curve inside the major loop, with Ha and M(Ha, Hb)."""
    _zero_axes(ax)
    H_up, M_up, H_dn, M_dn = major_loop_from_curves(curves)
    ax.plot(1e3 * H_up, MOMENT_SCALE * M_up, color="0.25", lw=1.4)
    ax.plot(1e3 * H_dn, MOMENT_SCALE * M_dn, color="0.25", lw=1.4)
    # symmetric about mu_0 H = 0, so the zero axis sits at the centre
    span = 1e3 * max(abs(H_up.min()), abs(H_dn.max())) * 1.04
    ax.set_xlim(-span, span)

    # baseline for the dropped guide lines, just under the lowest moment drawn
    floor = MOMENT_SCALE * float(np.min(M_up)) * 0.86

    curve = min(curves, key=lambda s: abs(s.Ha - HIGHLIGHT_HA))
    Hc_, Mc_ = np.asarray(curve.H, float), np.asarray(curve.M, float)
    ax.plot(1e3 * Hc_, MOMENT_SCALE * Mc_, color=CURVE_COLOR, lw=2.0, zorder=4)

    # the reversal field: where this curve leaves the descending branch
    Ha = float(curve.Ha)
    ax.plot([1e3 * Ha], [MOMENT_SCALE * Mc_[0]], "o", ms=6, mfc="white",
            mec=CURVE_COLOR, mew=1.6, zorder=6)
    ax.plot([1e3 * Ha, 1e3 * Ha], [floor, MOMENT_SCALE * Mc_[0]],
            color=CURVE_COLOR, ls=":", lw=1.0)
    ax.annotate(r"$H_a$", xy=(1e3 * Ha, floor), ha="center", va="bottom",
                fontsize=11, color=CURVE_COLOR)

    # one measured point partway along the curve
    k = int(0.45 * len(Hc_))
    Hb, Mb = Hc_[k], MOMENT_SCALE * Mc_[k]
    ax.plot([1e3 * Hb], [Mb], "o", ms=6, color=ACCENT, zorder=7)
    ax.plot([1e3 * Hb, 1e3 * Hb], [floor, Mb], color=ACCENT, ls=":", lw=1.0)
    ax.annotate(r"$H_b$", xy=(1e3 * Hb, floor), ha="center", va="bottom",
                fontsize=11, color=ACCENT)

    # Label the measured point, colouring each field symbol to match the
    # marker it refers to: Ha is the reversal field of this whole curve, not
    # something read off the vertical axis at this point.
    def _part(text, color):
        return TextArea(text, textprops=dict(color=color, fontsize=10))

    label = HPacker(children=[_part(r"$M($", "0.15"),
                              _part(r"$H_a$", CURVE_COLOR),
                              _part(r"$,$ ", "0.15"),
                              _part(r"$H_b$", ACCENT),
                              _part(r"$)$", "0.15")],
                    align="baseline", pad=0, sep=0)
    ax.add_artist(AnnotationBbox(
        label, (1e3 * Hb, Mb), xybox=(46, -22), xycoords="data",
        boxcoords="offset points", frameon=False, box_alignment=(0.0, 0.5),
        arrowprops=dict(arrowstyle="-", color="0.45", lw=0.8)))

    # direction of measurement
    j = int(0.72 * len(Hc_))
    ax.add_patch(FancyArrowPatch(
        (1e3 * Hc_[j], MOMENT_SCALE * Mc_[j]), (1e3 * Hc_[j + 6], MOMENT_SCALE * Mc_[j + 6]),
        arrowstyle="-|>", mutation_scale=13, color=CURVE_COLOR, zorder=8))

    ax.set_xlabel(r"$\mu_0 H$ (mT)")
    ax.set_ylabel(MOMENT_LABEL)
    ax.set_title("one first-order reversal curve", fontsize=10, loc="left")


def panel_b(ax, curves):
    """Draw the family of reversal curves that fills the major loop."""
    _zero_axes(ax)
    for s in curves:
        ax.plot(1e3 * np.asarray(s.H, float), MOMENT_SCALE * np.asarray(s.M, float),
                color="0.15", lw=0.3, alpha=0.6)
    ax.set_xlabel(r"$\mu_0 H$ (mT)")
    ax.set_ylabel(MOMENT_LABEL)
    ax.set_title(f"all {len(curves)} first-order reversal curves", fontsize=10,
                 loc="left")


def _grid_mesh(Ha_vals, Hb_vals):
    """Return pcolormesh edge arrays in (Hb, Ha) millitesla."""
    Ha_e, Hb_e = np.meshgrid(forc.centers_to_edges(Ha_vals),
                             forc.centers_to_edges(Hb_vals), indexing="ij")
    return 1e3 * Hb_e, 1e3 * Ha_e


def _attach_colorbar(ax, mappable, label):
    """Place a colorbar flush against an equal-aspect axes.

    Passing ``ax=`` to ``figure.colorbar`` reserves space from the axes' full
    bounding box, which for an equal-aspect panel is wider than the drawn
    square, leaving the bar floating well clear of the plot. An inset tracks
    the shrunk box instead.
    """
    cax = ax.inset_axes([1.03, 0.0, 0.045, 1.0])
    cb = ax.figure.colorbar(mappable, cax=cax)
    cb.set_label(label, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    return cb


def panel_c(ax, curves, limits):
    """Scatter the individual measurements in (Hb, Ha) coordinates.

    Every measured point of every reversal curve is one sample of the surface,
    plotted at its own (Hb, Ha) and coloured by its moment. Showing the points
    rather than an interpolated surface makes the sampling visible: the data
    are a set of horizontal rows, one per reversal curve, and they exist only
    where Hb >= Ha.
    """
    Hb = np.concatenate([np.asarray(s.H, float) for s in curves])
    Ha = np.concatenate([np.full(len(s.H), float(s.Ha)) for s in curves])
    M = np.concatenate([np.asarray(s.M, float) for s in curves])
    pcm = ax.scatter(1e3 * Hb, 1e3 * Ha, c=MOMENT_SCALE * M, s=0.8,
                     cmap="viridis", linewidths=0, rasterized=True)
    lo, hi = limits
    ax.plot([lo, hi], [lo, hi], color="0.25", lw=1.0, ls="--")
    ax.annotate("no data:\n$H_b < H_a$", xy=(-0.55 * hi, 0.42 * hi),
                fontsize=9, color="0.35", ha="center", va="center")

    # The other empty corner is not a physical boundary but an acquisition
    # one. This file's header sets the bias-axis limits (its Hb1 and Hb2 tags)
    # to -250 and +100 mT, and every curve stops where Bu = (Hb + Ha)/2
    # reaches +100 mT, which is the line drawn here.
    bu_max = 100.0
    ax.plot([lo, hi], [2 * bu_max - lo, 2 * bu_max - hi],
            color="0.45", lw=0.9, ls="-.")
    # Bu > bu_max lies ABOVE this line, and Hb >= Ha keeps us right of the
    # diagonal, so the empty wedge is the upper-right of the panel.
    ax.annotate("no data:\n" r"$B_u$ beyond the" "\nacquisition limit",
                xy=(0.66 * hi, 0.44 * hi), fontsize=9, color="0.35",
                ha="center", va="center")
    ax.annotate(r"$H_b = H_a$", xy=(0.62 * hi, 0.62 * hi), fontsize=9,
                color="0.25", rotation=45, ha="center", va="bottom")
    _square(ax, limits)
    ax.set_xlabel(r"$\mu_0 H_b$ (mT)")
    ax.set_ylabel(r"$\mu_0 H_a$ (mT)")
    ax.set_title(r"$M(H_a, H_b)$ measurement grid", fontsize=10, loc="left")
    _attach_colorbar(ax, pcm, MOMENT_LABEL)


def panel_d(ax, Ha_vals, Hb_vals, rho, limits):
    """Draw the FORC distribution on the same axes, with rotated Bc, Bu axes."""
    X, Y = _grid_mesh(Ha_vals, Hb_vals)
    vmax = forc.rho_window_vmax(Ha_vals, Hb_vals, rho, pct=100,
                                Bu_min=-0.15, Bu_max=0.075,
                                Bc_min=0.0, Bc_max=0.6)
    pcm = ax.pcolormesh(X, Y, rho / vmax, cmap=forc.get_forc_cmap(2),
                        vmin=-1, vmax=1, shading="auto")
    lo, hi = limits

    # Bc = (Hb - Ha)/2 and Bu = (Hb + Ha)/2, so in these axes Bc increases
    # along (1, -1) and Bu along (1, 1). The physical boundary Hb = Ha is the
    # Bc = 0 line, which is why it runs parallel to the Bu axis.
    ax.plot([lo, hi], [lo, hi], color="0.35", lw=0.9, ls="--", zorder=5)
    ax.plot([lo, hi], [-lo, -hi], color="0.55", lw=0.7, ls=":", zorder=5)

    for (dx, dy), label in (((1, 1), r"$B_u$"), ((1, -1), r"$B_c$")):
        tip = 0.93 * hi / np.sqrt(2)
        ax.add_patch(FancyArrowPatch(
            (0, 0), (dx * tip, dy * tip), arrowstyle="-|>", mutation_scale=9,
            color="0.3", lw=0.8, alpha=0.8, zorder=6))
        ax.annotate(label, xy=(dx * tip, dy * tip),
                    xytext=(5, 4 if dy > 0 else -12), textcoords="offset points",
                    fontsize=10, color="0.2", zorder=7)
    # No separate "Bc = 0" / "Bu = 0" labels: those lines are the Bu and Bc
    # axes themselves, so labelling both would say the same thing twice.
    ax.annotate(r"$H_b = H_a$", xy=(-0.66 * hi, -0.66 * hi), fontsize=8,
                color="0.4", rotation=45, ha="center", va="bottom")

    _square(ax, limits)
    ax.set_xlabel(r"$\mu_0 H_b$ (mT)")
    ax.set_ylabel(r"$\mu_0 H_a$ (mT)")
    ax.set_title(r"FORC distribution $\rho = -\frac{1}{2}\,"
                 r"\partial^2 M / \partial H_a \partial H_b$",
                 fontsize=10, loc="left")
    _attach_colorbar(ax, pcm, r"$\rho$ / max")


def _square(ax, limits):
    """Apply equal, shared limits so the 45 degree rotation reads correctly."""
    lo, hi = limits
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")


def main():
    curves, Ha_vals, Hb_vals, M_grid, rho = load_baraboo()

    # a window that holds the part of the half-plane where rho is resolved
    limits = (-450.0, 450.0)

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.2), constrained_layout=True)
    # the colorbars are inset axes, which constrained_layout does not account
    # for, so widen the gap between columns by hand
    fig.get_layout_engine().set(w_pad=0.10, wspace=0.10)
    panel_a(axes[0, 0], curves)
    panel_b(axes[0, 1], curves)
    panel_c(axes[1, 0], curves, limits)
    panel_d(axes[1, 1], Ha_vals, Hb_vals, rho, limits)

    for ax, label in zip(axes.ravel(), "abcd"):
        ax.annotate(f"({label})", xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(-38, 12), textcoords="offset points",
                    fontsize=13, fontweight="bold", va="top")

    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        fig.savefig(OUT_STEM.with_suffix(suffix), dpi=200, bbox_inches="tight")
    print(f"wrote {OUT_STEM.with_suffix('.png')} and .pdf")
    print(f"  {len(curves)} reversal curves, grid {M_grid.shape}, "
          f"peak |M| = {MOMENT_SCALE * np.nanmax(np.abs(M_grid)):.2f} uA m^2")


if __name__ == "__main__":
    main()
