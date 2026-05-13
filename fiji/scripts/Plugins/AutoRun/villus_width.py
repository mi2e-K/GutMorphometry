#!/usr/bin/env python3
"""
villus_width.py  –  Centerline-normal width measurement for intestinal villi.

Usage:
    python villus_width.py <image_tif> <pixel_um> <poly_json> <cl_json>
                           <output_txt> <qc_png> [<label>]

Arguments
---------
image_tif   : cropped TIF of the villus region
pixel_um    : pixel width in µm/px
poly_json   : JSON file with polygon coords [[x,y], ...]  (pixels, crop space)
cl_json     : JSON file with centerline coords [[x,y], ...] (pixels, crop space)
output_txt  : output file  →  mean_mid,w25,w50,w75,wmax,n_valid,qc_flag
qc_png      : QC overlay PNG (saved permanently next to CSV)
label       : optional label string (e.g. "282_L_1_outer_VW1")

Algorithm
---------
1.  Rasterise polygon → coarse mask (spatial constraint).
1b. If bg_rgb is available in musc_params.json: refine mask by OD-distance
    tissue segmentation within the polygon (Otsu on d_bg inside polygon,
    then morphological cleanup + largest component).  Falls back to raw
    polygon when bg_rgb is absent.
2.  Resample centerline at 1-px arc-length intervals; Gaussian smooth.
3.  Compute unit tangent/normal at every point via central differences.
4.  Cast rays in ±normal directions; find first mask exit → half-widths.
5.  Exclude tip/base 10%; collect valid width profile.
6.  Summarise: percentiles, mean over 20-80%, max.
7.  Save QC overlay: ROI contour (cyan), centerline (yellow),
    20/50/80% chords (green), invalid points (red).
    Title indicates whether tissue segmentation was used.
"""

import sys as _sys

if _sys.version_info[0] >= 3:

    import sys
    import json
    import numpy as np
    from scipy.ndimage import gaussian_filter1d
    from scipy.interpolate import interp1d
    import os
    from skimage import io, filters
    from skimage.draw import polygon as sk_polygon
    from skimage.measure import find_contours, label as sk_label
    from skimage.morphology import binary_closing, disk, remove_small_objects
    from scipy.ndimage import map_coordinates, binary_fill_holes

    # ------------------------------------------------------------------
    # Background-reference loading
    # ------------------------------------------------------------------

    def _load_bg_rgb():
        """
        Read bg_rgb from musc_params.json (written by tune_muscularis.py).
        Returns a list [R, G, B] or None if not available.
        The file lives at <Fiji-root>/GutMorphometry_config/musc_params.json.
        """
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fiji_root  = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
        params_path = os.path.join(fiji_root, "GutMorphometry_config",
                                   "musc_params.json")
        if not os.path.exists(params_path):
            return None
        try:
            with open(params_path) as _f:
                p = json.load(_f)
            return p.get("bg_rgb", None)
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Tissue segmentation within polygon
    # ------------------------------------------------------------------

    def _refine_mask(img_rgb, poly_mask, bg_rgb):
        """
        Refine coarse polygon mask to true tissue boundary.

        Within the polygon, compute OD Euclidean distance from the background
        colour.  Pixels with d_bg > Otsu-threshold(d_bg inside polygon) are
        classified as tissue.  Morphological cleanup is then applied:

            binary_closing(disk(3))
            → binary_fill_holes          (include lumen in width boundary)
            → remove_small_objects(64)
            → keep largest component
            → binary_fill_holes
            → intersect with polygon     (never go outside the user's ROI)

        Falls back to poly_mask if segmentation yields an empty mask.
        """
        bg_px = np.clip(np.array(bg_rgb[:3], dtype=np.float64), 1.0, 255.0)
        bg_od = -np.log10(bg_px / 255.0)

        rgb_f = np.clip(img_rgb[:, :, :3].astype(np.float64), 1.0, 255.0)
        od    = -np.log10(rgb_f / 255.0)

        # Euclidean OD distance from background at every pixel
        d_bg = np.sqrt(((od - bg_od) ** 2).sum(axis=2))

        # Otsu threshold on distances inside the polygon
        d_vals = d_bg[poly_mask]
        if d_vals.size < 20:
            return poly_mask
        try:
            thresh = filters.threshold_otsu(d_vals)
        except Exception:
            return poly_mask

        tissue = poly_mask & (d_bg > thresh)

        # Morphological cleanup
        tissue = binary_closing(tissue, disk(3))
        tissue = binary_fill_holes(tissue)
        tissue = remove_small_objects(tissue, min_size=64)

        # Keep largest connected component
        labeled = sk_label(tissue)
        if labeled.max() == 0:
            return poly_mask          # fallback: segmentation gave nothing
        sizes    = np.bincount(labeled.ravel())
        sizes[0] = 0
        tissue   = labeled == sizes.argmax()
        tissue   = binary_fill_holes(tissue)

        # Must stay inside the polygon (no expansion outside user's ROI)
        tissue = tissue & poly_mask

        if not tissue.any():
            return poly_mask          # fallback if result is empty
        return tissue.astype(bool)

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    def _rasterize_polygon(poly_x, poly_y, shape):
        """Binary mask from polygon vertex arrays (pixel coords)."""
        mask = np.zeros(shape[:2], dtype=bool)
        rr, cc = sk_polygon(poly_y, poly_x, shape=shape[:2])
        mask[rr, cc] = True
        return mask

    def _resample_centerline(pts, step=1.0):
        """
        Resample a polyline at equal arc-length intervals of `step` pixels.
        Returns (resampled_pts, s_values) where s_values is cumulative arc length.
        """
        pts = np.asarray(pts, dtype=float)
        diffs = np.diff(pts, axis=0)
        seg_len = np.sqrt((diffs ** 2).sum(axis=1))
        cum_len = np.concatenate([[0.0], np.cumsum(seg_len)])
        total = float(cum_len[-1])
        if total < 2.0:
            return pts, cum_len

        # Remove duplicate s values to make interp1d happy
        _, unique_idx = np.unique(cum_len, return_index=True)
        s_u = cum_len[unique_idx]
        p_u = pts[unique_idx]

        n = max(10, int(total / step))
        s_new = np.linspace(0.0, total, n)
        fx = interp1d(s_u, p_u[:, 0], kind="linear")
        fy = interp1d(s_u, p_u[:, 1], kind="linear")
        resampled = np.column_stack([fx(s_new), fy(s_new)])
        return resampled, s_new

    def _smooth_centerline(pts, total_len):
        """Gaussian smooth. sigma ≈ 3 % of total arc length (min 2 px)."""
        sigma = max(2.0, total_len * 0.03)
        xs = gaussian_filter1d(pts[:, 0], sigma)
        ys = gaussian_filter1d(pts[:, 1], sigma)
        return np.column_stack([xs, ys])

    def _compute_normals(pts):
        """
        Unit tangent and left-hand normal at every point (central differences).
        Returns tx, ty, nx, ny as 1-D arrays of length len(pts).
        """
        dx = np.gradient(pts[:, 0])
        dy = np.gradient(pts[:, 1])
        mag = np.hypot(dx, dy)
        mag = np.where(mag < 1e-9, 1e-9, mag)
        tx =  dx / mag
        ty =  dy / mag
        nx = -ty          # perpendicular to tangent
        ny =  tx
        return tx, ty, nx, ny

    def _ray_exit_distance(mask, px, py, dx, dy, max_dist=800, step=0.5):
        """
        Distance to first background pixel along direction (dx, dy) from (px, py).
        Returns None if the start point is outside the mask or no exit is found
        within max_dist.  Uses sub-pixel linear interpolation.
        """
        h, w = mask.shape
        t = np.arange(0.0, max_dist, step)
        xs = px + t * dx
        ys = py + t * dy
        in_bounds = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
        vals = map_coordinates(
            mask.astype(np.float32), [ys, xs],
            order=1, mode="constant", cval=0.0
        )
        # First step outside mask (or out of image)
        exits = np.where((vals < 0.5) | ~in_bounds)[0]
        if len(exits) == 0:
            return None
        if exits[0] == 0:
            return None   # start already outside
        i = exits[0]
        # Linear sub-pixel interpolation
        v0, v1 = float(vals[i - 1]), float(vals[i])
        t_exit = t[i - 1] + step * (0.5 - v0) / (v0 - v1 + 1e-9)
        return float(max(0.0, t_exit))

    # ------------------------------------------------------------------
    # Core measurement
    # ------------------------------------------------------------------

    def _compute_width_profile(mask, cl_pts, nx, ny, exclude_frac=0.10):
        """
        Boundary-to-boundary width at each centerline sample point.

        Returns
        -------
        widths    : (N,) float array, NaN for invalid points
        left_pts  : (N, 2) exit points in the +normal direction (NaN if invalid)
        right_pts : (N, 2) exit points in the -normal direction (NaN if invalid)
        """
        n = len(cl_pts)
        widths    = np.full(n, np.nan, dtype=float)
        left_pts  = np.full((n, 2), np.nan, dtype=float)
        right_pts = np.full((n, 2), np.nan, dtype=float)

        i_start = int(np.ceil(n * exclude_frac))
        i_end   = int(np.floor(n * (1.0 - exclude_frac)))

        for i in range(i_start, i_end):
            px, py = cl_pts[i, 0], cl_pts[i, 1]

            # Skip if centerline point is outside the mask
            ix, iy = int(round(px)), int(round(py))
            if not (0 <= iy < mask.shape[0] and 0 <= ix < mask.shape[1]):
                continue
            if not mask[iy, ix]:
                continue

            d_pos = _ray_exit_distance(mask, px, py,  nx[i],  ny[i])
            d_neg = _ray_exit_distance(mask, px, py, -nx[i], -ny[i])

            if d_pos is None or d_neg is None:
                continue

            widths[i]    = d_pos + d_neg
            left_pts[i]  = [px + d_pos * nx[i], py + d_pos * ny[i]]
            right_pts[i] = [px - d_neg * nx[i], py - d_neg * ny[i]]

        return widths, left_pts, right_pts

    def _summarize(widths, s_values, pixel_um):
        """
        Summary statistics from width profile.
        Returns None when fewer than 3 valid samples exist.
        """
        valid = ~np.isnan(widths)
        n_valid = int(valid.sum())
        if n_valid < 3:
            return None

        total = float(s_values[-1])

        def _width_at_frac(frac):
            s_tgt = frac * total
            idx   = int(np.argmin(np.abs(s_values - s_tgt)))
            if not np.isnan(widths[idx]):
                return float(widths[idx] * pixel_um)
            # Nearest valid neighbour
            valid_idx = np.where(valid)[0]
            near = valid_idx[np.argmin(np.abs(valid_idx - idx))]
            return float(widths[near] * pixel_um)

        mid_mask = valid & (s_values >= 0.20 * total) & (s_values <= 0.80 * total)
        mean_mid = (float(np.nanmean(widths[mid_mask]) * pixel_um)
                    if mid_mask.any() else float("nan"))

        return {
            "width_mean_mid_um": mean_mid,
            "width_20_um":       _width_at_frac(0.20),
            "width_50_um":       _width_at_frac(0.50),
            "width_80_um":       _width_at_frac(0.80),
            "width_max_um":      float(np.nanmax(widths[valid]) * pixel_um),
            "n_valid":           n_valid,
            "qc_flag":           0,
        }

    # ------------------------------------------------------------------
    # QC overlay
    # ------------------------------------------------------------------

    def _save_qc(img_rgb, mask, cl_smooth, widths, left_pts, right_pts,
                 s_values, stats, label, qc_path, seg_mode="polygon boundary"):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 8))
        fig.patch.set_facecolor("#1a1a1a")
        ax.set_facecolor("#1a1a1a")
        ax.imshow(img_rgb[:, :, :3])

        # ROI boundary contour — cyan
        for contour in find_contours(mask.astype(float), 0.5):
            ax.plot(contour[:, 1], contour[:, 0],
                    color="cyan", linewidth=1.5, label="ROI boundary")

        # Centerline — yellow
        ax.plot(cl_smooth[:, 0], cl_smooth[:, 1],
                color="yellow", linewidth=1.5, label="Centerline")

        total = float(s_values[-1])

        # Dashed boundary lines at 20% and 80% — range used for mean_mid
        # (disabled — kept for reference)
        # h_img, w_img = mask.shape
        # _, _, nx_qc, ny_qc = _compute_normals(cl_smooth)
        # for frac in (0.20, 0.80):
        #     idx = int(np.argmin(np.abs(s_values - frac * total)))
        #     cx, cy = cl_smooth[idx, 0], cl_smooth[idx, 1]
        #     span = max(h_img, w_img)
        #     x0 = cx - span * nx_qc[idx];  y0 = cy - span * ny_qc[idx]
        #     x1 = cx + span * nx_qc[idx];  y1 = cy + span * ny_qc[idx]
        #     ax.plot([x0, x1], [y0, y1],
        #             color="white", linewidth=1.2,
        #             linestyle="--", dashes=(6, 4), alpha=0.7)

        # Width chords at 20 / 50 / 80 % — shades of green
        for frac, col, tag in [(0.20, "#66ffaa", "20%"),
                                (0.50, "#00dd55", "50%"),
                                (0.80, "#009933", "80%")]:
            idx = int(np.argmin(np.abs(s_values - frac * total)))
            if not np.isnan(widths[idx]):
                lp = left_pts[idx]
                rp = right_pts[idx]
                ax.plot([rp[0], lp[0]], [rp[1], lp[1]],
                        color=col, linewidth=2.0)
                mid = cl_smooth[idx]
                ax.text(mid[0] + 4, mid[1], tag,
                        color=col, fontsize=7, va="center",
                        bbox=dict(facecolor="#00000080", pad=1, edgecolor="none"))

        # Invalid points in the measured range — red
        n = len(widths)
        i_start = int(np.ceil(n * 0.10))
        i_end   = int(np.floor(n * 0.90))
        invalid = (np.isnan(widths) &
                   (np.arange(n) >= i_start) &
                   (np.arange(n) <  i_end))
        if invalid.any():
            ax.scatter(cl_smooth[invalid, 0], cl_smooth[invalid, 1],
                       c="red", s=18, zorder=5, label="Invalid")

        # Title
        if stats:
            title = ("{}\n"
                     "mean_mid={:.1f} µm  "
                     "w20={:.1f}  w50={:.1f}  w80={:.1f}  wmax={:.1f}\n"
                     "n_valid={}  qc_flag={}  [{}]").format(
                label or "Villus Width",
                stats["width_mean_mid_um"],
                stats["width_20_um"], stats["width_50_um"],
                stats["width_80_um"], stats["width_max_um"],
                stats["n_valid"], stats["qc_flag"], seg_mode)
        else:
            title = (label or "Villus Width") + "\n(insufficient valid measurements)\n[{}]".format(seg_mode)

        ax.set_title(title, fontsize=9, color="white", pad=6)
        ax.axis("off")

        handles = [
            plt.Line2D([0], [0], color="cyan",    lw=1.5, label="ROI boundary"),
            plt.Line2D([0], [0], color="yellow",  lw=1.5, label="Centerline"),
            plt.Line2D([0], [0], color="#00dd55", lw=2.0, label="Width chord"),
            plt.Line2D([0], [0], marker="o", color="red", lw=0,
                       markersize=5, label="Invalid"),
        ]
        ax.legend(handles=handles, loc="lower right", fontsize=7,
                  facecolor="#333333", labelcolor="white", framealpha=0.8)

        plt.tight_layout()
        plt.savefig(qc_path, dpi=150, bbox_inches="tight",
                    facecolor="#1a1a1a")
        plt.close(fig)

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    def main():
        n = len(sys.argv)
        if n < 7:
            sys.stderr.write(
                "Usage: villus_width.py <image> <pixel_um> <poly_json> "
                "<cl_json> <output_txt> <qc_png> [<label>]\n"
            )
            sys.exit(1)

        image_path = sys.argv[1]
        pixel_um   = float(sys.argv[2])
        poly_path  = sys.argv[3]
        cl_path    = sys.argv[4]
        out_path   = sys.argv[5]
        qc_path    = sys.argv[6]
        label      = sys.argv[7] if n > 7 else ""

        # Load image
        img = io.imread(image_path)
        if img.ndim < 3:
            sys.stderr.write("ERROR: Expected an RGB image.\n")
            sys.exit(1)
        h, w = img.shape[:2]

        # Load coordinates
        with open(poly_path) as f:
            poly_coords = json.load(f)   # [[x, y], ...]
        with open(cl_path) as f:
            cl_coords = json.load(f)     # [[x, y], ...]

        poly_x = np.array([p[0] for p in poly_coords], dtype=float)
        poly_y = np.array([p[1] for p in poly_coords], dtype=float)
        cl_pts = np.array([[p[0], p[1]] for p in cl_coords], dtype=float)

        # Coarse polygon mask
        poly_mask = _rasterize_polygon(poly_x, poly_y, (h, w))

        # Refine to tissue boundary if background reference is available
        bg_rgb = _load_bg_rgb()
        if bg_rgb is not None:
            mask = _refine_mask(img, poly_mask, bg_rgb)
            seg_mode = "tissue-seg (bg_rgb={})".format(bg_rgb)
        else:
            mask = poly_mask
            seg_mode = "polygon boundary (no bg_rgb)"
        sys.stdout.write("  Mask mode: {}\n".format(seg_mode))

        # Resample and smooth centerline
        cl_resampled, s_values = _resample_centerline(cl_pts.tolist(), step=1.0)
        total_len = float(s_values[-1]) if len(s_values) > 1 else 0.0
        if total_len < 5.0:
            sys.stderr.write("ERROR: Centerline too short ({:.1f} px).\n"
                             .format(total_len))
            with open(out_path, "w") as f:
                f.write("nan,nan,nan,nan,nan,0,1")
            sys.exit(0)

        cl_smooth = _smooth_centerline(cl_resampled, total_len)

        # Normals
        _tx, _ty, nx, ny = _compute_normals(cl_smooth)

        # Width profile
        widths, left_pts, right_pts = _compute_width_profile(
            mask, cl_smooth, nx, ny, exclude_frac=0.10
        )

        # Summary statistics
        stats = _summarize(widths, s_values, pixel_um)

        # Write output
        if stats is None:
            sys.stderr.write("WARNING: Insufficient valid measurements.\n")
            with open(out_path, "w") as f:
                f.write("nan,nan,nan,nan,nan,0,1")
        else:
            with open(out_path, "w") as f:
                f.write("{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:d},{:d}".format(
                    stats["width_mean_mid_um"],
                    stats["width_20_um"],
                    stats["width_50_um"],
                    stats["width_80_um"],
                    stats["width_max_um"],
                    stats["n_valid"],
                    stats["qc_flag"],
                ))

        # QC image
        _save_qc(img, mask, cl_smooth, widths, left_pts, right_pts,
                 s_values, stats, label, qc_path, seg_mode=seg_mode)

        if stats:
            sys.stdout.write(
                "Width mean_mid={:.1f} µm  w20={:.1f}  w50={:.1f}  "
                "w80={:.1f}  wmax={:.1f}  n_valid={}\n".format(
                    stats["width_mean_mid_um"],
                    stats["width_20_um"], stats["width_50_um"],
                    stats["width_80_um"], stats["width_max_um"],
                    stats["n_valid"],
                )
            )

    if __name__ == "__main__":
        main()
