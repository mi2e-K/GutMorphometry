#!/usr/bin/env python3
"""
musc_thickness.py  --  Mean muscularis thickness from a polygon-masked H&E crop.

Usage:
    python musc_thickness.py <image_tif> <pixel_width_um> <output_txt> <viz_png>

Arguments:
    image_tif        Polygon-masked H&E TIF saved by the Fiji macro.
                     Pixels outside the user's polygon are already white
                     (background), so no extra spatial constraint is needed.
    pixel_width_um   Pixel width in micrometres (from Fiji scale calibration)
    output_txt       Path where the result (single float, um) will be written
    viz_png          Path where the 4-panel verification PNG will be saved

Algorithm:
    1. H&E colour deconvolution (Ruifrok & Johnston 2001)
    2. Tissue mask: exclude near-white background (OD sum < 0.15)
    3. E-H feature thresholding:
         feature = Eosin - Haematoxylin
         muscle (pink)  : E-H >> 0   crypt (purple) : E-H << 0
         baseline threshold = midpoint of user-picked pink/purple reference
         E-H values (stored in musc_params.json); falls back to Otsu if not set.
         e_factor slider shifts the baseline: +scale per unit above 1.0.
    4. Morphological closing bridges muscle-cell nuclei inside the band.
    5. Keep only the largest connected component; EDT mean local thickness.

Dependencies (all included with Anaconda):
    numpy, scikit-image, scipy, matplotlib
"""

# ------------------------------------------------------------------
# AutoRun cleanup — runs under BOTH Jython (Fiji startup) and CPython.
#
# Fiji's AutoRun folder tries to execute every file it finds there,
# including .json files.  If musc_params.json is in AutoRun, Fiji's
# macro interpreter chokes on it at startup with:
#   "Statement cannot begin with ':' in line 2"
# This block deletes that file at startup so the error disappears on
# the next Fiji launch.  Uses only os module — safe for Jython/Python 2.
# ------------------------------------------------------------------
try:
    import os as _os
    try:
        _autorun_dir = _os.path.dirname(_os.path.abspath(__file__))
    except (NameError, TypeError):
        _autorun_dir = None
    if _autorun_dir:
        # 1. Remove stray musc_params.json directly in AutoRun
        _old_json = _os.path.join(_autorun_dir, "musc_params.json")
        if _os.path.exists(_old_json):
            _os.remove(_old_json)
        # 2. Remove GutMorphometry_config subdir that was accidentally created
        #    inside AutoRun (Fiji tries to execute it and shows "access denied")
        _wrong_config = _os.path.join(_autorun_dir, "GutMorphometry_config")
        _wrong_json   = _os.path.join(_wrong_config, "musc_params.json")
        if _os.path.exists(_wrong_json):
            _os.remove(_wrong_json)
        try:
            _os.rmdir(_wrong_config)   # only succeeds when empty
        except Exception:
            pass
except Exception:
    pass

# ------------------------------------------------------------------
# Jython / Python-2 guard
# Fiji's AutoRun folder runs all .py files at startup using Jython
# (Python 2).  Wrapping everything in this version check means
# Jython finds only a one-line no-op — no imports, no errors,
# no sys.exit() traceback.  Python 3 enters the block normally.
# ------------------------------------------------------------------
import sys as _sys

if _sys.version_info[0] >= 3:

    import sys
    import json
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage import io, filters
    from skimage.measure import label
    from skimage.morphology import remove_small_objects, binary_closing, disk, medial_axis
    from scipy.ndimage import distance_transform_edt, binary_fill_holes, gaussian_filter, convolve

    # ------------------------------------------------------------------
    # Standard H&E stain vectors -- Ruifrok & Johnston (2001)
    # ------------------------------------------------------------------
    _HE_MATRIX = np.array([
        [0.6500286, 0.7040924, 0.2860126],   # Haematoxylin
        [0.0704327, 0.9778137, 0.1972729],   # Eosin
        [0.7110272, 0.4232444, 0.5603568],   # Residual
    ], dtype=np.float64)

    _HE_INV = np.linalg.inv(_HE_MATRIX)

    def deconvolve_he(rgb):
        """
        Separate H&E stains via optical-density deconvolution.
        Returns (H_map, E_map) as float32 arrays.
        np.dot used instead of @ for broader tool compatibility.
        """
        rgb_f = rgb[:, :, :3].astype(np.float64)
        rgb_f = np.clip(rgb_f, 1, 255)
        od    = -np.log10(rgb_f / 255.0)

        stains = np.dot(od.reshape(-1, 3), _HE_INV.T)
        stains = np.clip(stains, 0, None)

        h = stains[:, 0].reshape(rgb.shape[:2]).astype(np.float32)
        e = stains[:, 1].reshape(rgb.shape[:2]).astype(np.float32)
        return h, e

    def get_tissue_mask(rgb):
        """
        True where tissue is present (not background).
        Background (including polygon exterior set to white by Fiji)
        has total optical density < 0.15 and is excluded.
        """
        rgb_f = np.clip(rgb[:, :, :3].astype(np.float64), 1, 255)
        od    = -np.log10(rgb_f / 255.0)
        return od.sum(axis=2) >= 0.15

    def _load_params():
        """
        Load tuning parameters from musc_params.json.

        Returns (e_factor, pink_rgb, purple_rgb).
          e_factor   : slider factor (default 1.0)
          pink_rgb   : [R, G, B] of user-picked pink (muscle) reference pixel,
                       or None if not set
          purple_rgb : [R, G, B] of user-picked purple (crypt) reference pixel,
                       or None if not set

        When both pink_rgb and purple_rgb are present, the mask is computed by
        nearest-centroid classification in optical-density space — independent
        of whether the target tissue has higher or lower E-H than background.
        Falls back to Otsu on E-H when references are absent.

        Params are stored at [Fiji root]/GutMorphometry_config/musc_params.json.
        """
        _autorun   = os.path.dirname(os.path.abspath(__file__))
        _fiji_root = os.path.dirname(os.path.dirname(os.path.dirname(_autorun)))
        params_path = os.path.join(_fiji_root, "GutMorphometry_config", "musc_params.json")

        if os.path.exists(params_path):
            try:
                with open(params_path) as f:
                    p = json.load(f)
                pink_rgb   = p.get("pink_rgb",   None)
                purple_rgb = p.get("purple_rgb", None)
                bg_rgb     = p.get("bg_rgb",     None)
                # Ignore old-format pink_eh/purple_eh keys (pre-OD-mode params)
                return float(p.get("e_factor", 1.0)), pink_rgb, purple_rgb, bg_rgb
            except Exception:
                pass
        return 1.0, None, None, None

    def muscularis_mask(rgb, pixel_width_um=1.0):
        """
        Return a boolean mask: True = muscularis.

        Thresholding on the E-H feature (Eosin minus Haematoxylin):

          muscle (pink)  : high E, low  H  ->  E-H >> 0
          crypt  (purple): low  E, high H  ->  E-H << 0
          background     : excluded by tissue mask before thresholding

        A single Otsu threshold on E-H within tissue gives a sharper
        pink/purple boundary than separate E-gate and H-gate thresholds,
        because both channels are used jointly in one discriminating
        feature rather than independently.

        e_factor (from musc_params.json, default 1.0):
          > 1  →  stricter (requires stronger pink-over-purple signal)
          < 1  →  more inclusive
        """
        e_factor, pink_rgb, purple_rgb, bg_rgb = _load_params()

        h, e   = deconvolve_he(rgb)
        tissue = get_tissue_mask(rgb)

        if pink_rgb is not None and purple_rgb is not None:
            # --- OD nearest-centroid mode ---
            # Classify each tissue pixel as muscle or crypt by comparing its
            # OD distance to each reference.  Works even when the target tissue
            # does not have the highest E-H in the image.
            def _rgb_to_od(rgb_pixel):
                px = np.clip(np.array(rgb_pixel[:3], dtype=np.float64), 1.0, 255.0)
                return -np.log10(px / 255.0)

            rgb_f    = np.clip(rgb[:, :, :3].astype(np.float64), 1.0, 255.0)
            od_img   = -np.log10(rgb_f / 255.0)          # (H, W, 3)
            pink_od  = _rgb_to_od(pink_rgb)
            purp_od  = _rgb_to_od(purple_rgb)
            d_pink   = np.sqrt(((od_img - pink_od)  ** 2).sum(axis=2))
            d_purple = np.sqrt(((od_img - purp_od)  ** 2).sum(axis=2))

            # Scale = 20% of OD separation between references
            od_sep   = float(np.linalg.norm(pink_od - purp_od))
            od_scale = max(0.02, od_sep * 0.20)
            margin   = (e_factor - 1.0) * od_scale
            mask = tissue & ((d_purple - d_pink) > margin)
            if bg_rgb is not None:
                # 3-class: also require closer to pink than background
                bg_od  = _rgb_to_od(bg_rgb)
                d_bg   = np.sqrt(((od_img - bg_od) ** 2).sum(axis=2))
                mask = mask & ((d_bg - d_pink) > margin)
        else:
            # --- Otsu fallback on E-H ---
            feature   = e - h
            feat_vals = feature[tissue]
            if feat_vals.size > 100:
                feat_base  = filters.threshold_otsu(feat_vals)
                feat_iqr   = float(np.percentile(feat_vals, 75) -
                                   np.percentile(feat_vals, 25))
                feat_scale = max(0.05, feat_iqr * 0.30)
                feat_thresh = feat_base + (e_factor - 1.0) * feat_scale
            else:
                feat_thresh = 0.0
            mask = tissue & (feature > feat_thresh)

        # Closing: bridge small H-rich nuclei that sit inside the muscle band.
        # Target bridge radius = 8 um (smooth-muscle nucleus ~4-8 um).
        r_bridge = max(3, min(20, int(round(8.0 / pixel_width_um))))
        mask = binary_closing(mask, disk(r_bridge))

        # Fill any residual holes within the bridged band
        mask = binary_fill_holes(mask)

        # Light final closing for smooth edges — done BEFORE size filtering
        # so that closing-merged regions are counted as one component.
        mask = binary_closing(mask, disk(3))

        # Remove fragments smaller than a (20 um)^2 area.
        # Runs after both closings so merged clusters count their full size.
        min_px = max(64, int(round(400.0 / (pixel_width_um ** 2))))
        mask = remove_small_objects(mask, min_size=min_px)

        # Keep only the single largest connected component.
        # Guarantees exactly one output region (muscularis propria).
        # Muscularis mucosae, loose smooth-muscle fibres, and artefacts
        # are all discarded here regardless of the e_factor setting.
        labeled = label(mask)
        if labeled.max() > 0:
            sizes    = np.bincount(labeled.ravel())
            sizes[0] = 0
            mask     = labeled == sizes.argmax()

        # Fill internal holes in the final component.
        mask = binary_fill_holes(mask)

        # Gentle boundary smoothing (sigma=2 px ≈ rounds jags < ~4 px wide).
        # Threshold at 0.5: equivalent to rounding corners without shrinking.
        _gs  = gaussian_filter(mask.astype(np.float32), sigma=2.0)
        mask = binary_fill_holes(_gs >= 0.5)

        return mask.astype(bool)

    def _prune_skel_spurs(skel, n_iter=10):
        """
        Strip degree-1 pixels (endpoints) from a skeleton iteratively.

        Each pass removes all pixels that have exactly one 8-connected skeleton
        neighbour.  Repeating n_iter times removes spurs up to n_iter pixels
        long — suppresses end-cap fanning and boundary-noise branches while
        leaving the main centerline intact.
        """
        kernel = np.ones((3, 3), dtype=np.uint8)
        kernel[1, 1] = 0
        skel = skel.copy()
        for _ in range(n_iter):
            nbrs = convolve(skel.astype(np.uint8), kernel, mode="constant", cval=0)
            endpoints = skel & (nbrs == 1)
            if not endpoints.any():
                break
            skel[endpoints] = False
        return skel

    def mean_local_thickness_um(mask, pixel_width_um):
        """
        Medial-axis EDT mean thickness in µm.

        Returns (mean_um, edt_pixel_array, pruned_skel).

        Algorithm:
          1. medial_axis(mask, return_distance=True): exact Euclidean medial
             axis — the locus of centres of maximal inscribed circles.
             `distance[skel]` gives the inscribed-circle radius at each axis
             pixel directly, without a separate EDT call.
          2. Prune short spurs via iterative endpoint removal (10 passes) to
             suppress end-cap fanning and boundary-noise branches.
          3. Mean of 2 × radius over the pruned axis = mean wall thickness.
             For a uniform ribbon of true width W every axis pixel has
             radius = W/2, so mean(2 × radius) = W  ✓
             Compare: mean(2 × EDT over ALL mask pixels) = W/2 for a uniform
             strip — the 50 % underestimate this replaces.
          4. A full EDT map is computed separately for the heatmap.
        """
        skel_raw, dist = medial_axis(mask, return_distance=True)
        skel = _prune_skel_spurs(skel_raw, n_iter=10)

        radii = dist[skel]
        if radii.size == 0:
            # Fallback: mask too thin for medial axis; use full-pixel EDT mean
            edt = distance_transform_edt(mask)
            fg  = edt[mask]
            mean_um = float(2.0 * np.mean(fg) * pixel_width_um) if fg.size > 0 else 0.0
            return mean_um, edt, skel_raw

        mean_um = float(2.0 * np.mean(radii) * pixel_width_um)

        # Full EDT for the heatmap visualization (medial_axis dist is 0 off-axis)
        edt = distance_transform_edt(mask)
        return mean_um, edt, skel

    def save_visualization(img_rgb, mask, edt_px, pixel_width_um,
                           mean_um, viz_path, skel=None):
        """
        4-panel verification PNG.

        Panel 1  Polygon-masked crop (original colours)
        Panel 2  Detected muscle mask overlaid in green; skeleton in yellow
        Panel 3  2*EDT map (local inscribed-circle diameter, µm)
                 Cyan dots = skeleton pixels where thickness is sampled
        Panel 4  E-H feature map with detected boundary contour (yellow)
                 Red = eosin-dominant (muscle); Blue = haematoxylin-dominant
                 (crypt/mucosa).  Yellow line = detected mask edge.
        """
        h, e    = deconvolve_he(img_rgb)
        feature = e - h

        fig, axes = plt.subplots(1, 4, figsize=(17, 4))

        # --- Panel 1: original crop ---
        axes[0].imshow(img_rgb[:, :, :3])
        axes[0].set_title("Masked crop (polygon ROI)", fontsize=9)
        axes[0].axis("off")

        # --- Panel 2: mask overlay + skeleton ---
        overlay = img_rgb[:, :, :3].astype(np.int16)
        overlay[mask, 1] = np.clip(overlay[mask, 1] + 50, 0, 255)
        overlay[mask, 0] = np.clip(overlay[mask, 0] - 20, 0, 255)
        overlay[mask, 2] = np.clip(overlay[mask, 2] - 20, 0, 255)
        axes[1].imshow(overlay.astype(np.uint8))
        if skel is not None and skel.any():
            ys, xs = np.where(skel)
            axes[1].scatter(xs, ys, c="yellow", s=0.3, linewidths=0)
        axes[1].set_title("Detected mask (green)\nskeleton (yellow)", fontsize=9)
        axes[1].axis("off")

        # --- Panel 3: 2*EDT map with skeleton sampling points ---
        thickness_map        = 2.0 * edt_px * pixel_width_um
        thickness_map[~mask] = 0.0
        vmax = float(np.percentile(thickness_map[mask], 95)) if mask.any() else 1.0
        im   = axes[2].imshow(thickness_map, cmap="hot", vmin=0, vmax=vmax)
        if skel is not None and skel.any():
            ys, xs = np.where(skel)
            axes[2].scatter(xs, ys, c="cyan", s=0.3, linewidths=0)
        axes[2].set_title(
            "2×EDT map (inscribed-circle diam, µm)\nmean along skeleton = {:.1f} µm".format(mean_um),
            fontsize=9)
        axes[2].axis("off")
        plt.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04, label="um")

        # --- Panel 4: E-H feature map + detected boundary ---
        abs98 = float(np.percentile(np.abs(feature), 98))
        axes[3].imshow(feature, cmap="RdBu_r", vmin=-abs98, vmax=abs98)
        if mask.any():
            axes[3].contour(mask, levels=[0.5], colors="yellow", linewidths=1.5)
        axes[3].set_title(
            "E-H feature  (red=muscle, blue=crypt)\nyellow = detected boundary",
            fontsize=9)
        axes[3].axis("off")

        plt.tight_layout()
        plt.savefig(viz_path, dpi=100, bbox_inches="tight")
        plt.close(fig)

    def save_result_image(img_rgb, mask, edt_px, pixel_width_um,
                          mean_um, result_path, label="", length_um=0.0, skel=None):
        """
        2-panel result image saved alongside the CSV.

        Left  : original crop with detected boundary contour (cyan) and
                skeleton (yellow dots).
        Right : 2*EDT heatmap (local inscribed-circle diameter, µm).
                Cyan dots = skeleton pixels used for mean thickness.

        This image is kept permanently — filename encodes the animal ID,
        region, section, loop position, and measurement index so it can
        be identified without opening the CSV.
        """
        thickness_map        = 2.0 * edt_px * pixel_width_um
        thickness_map[~mask] = np.nan

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.patch.set_facecolor("white")

        title_str = label.replace("_", "  ") if label else "muscularis"

        # --- Panel 1: original + boundary contour + skeleton ---
        axes[0].imshow(img_rgb[:, :, :3])
        if mask.any():
            axes[0].contour(mask, levels=[0.5], colors="cyan", linewidths=2.0)
        if skel is not None and skel.any():
            ys, xs = np.where(skel)
            axes[0].scatter(xs, ys, c="yellow", s=0.5, linewidths=0)
        axes[0].set_title(
            "Boundary (cyan) + skeleton (yellow)\n[{}]  mean={:.1f} µm  len={:.0f} µm".format(
                title_str, mean_um, length_um),
            fontsize=10)
        axes[0].axis("off")

        # --- Panel 2: 2*EDT heatmap + skeleton sampling points ---
        vmax = float(np.nanpercentile(thickness_map[mask], 95)) if mask.any() else 1.0
        im   = axes[1].imshow(thickness_map, cmap="hot", vmin=0, vmax=vmax)
        if skel is not None and skel.any():
            ys, xs = np.where(skel)
            axes[1].scatter(xs, ys, c="cyan", s=0.5, linewidths=0)
        axes[1].set_title(
            "2×EDT map (inscribed-circle diam, µm)\nmean along skeleton = {:.1f} µm".format(mean_um),
            fontsize=10)
        axes[1].axis("off")
        plt.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04, label="um")

        plt.tight_layout()
        plt.savefig(result_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

    def main():
        n = len(sys.argv)
        if n < 5:
            sys.stderr.write(
                "Usage: python musc_thickness.py "
                "<image_tif> <pixel_width_um> <output_txt> <viz_png>\n"
                "  or:  musc_thickness.py "
                "<image_tif> <pixel_width_um> <boundary_frac> "
                "<output_txt> <viz_png> [<label>] [<result_png>]\n"
            )
            sys.exit(1)

        image_path     = sys.argv[1]
        pixel_width_um = float(sys.argv[2])

        # Detect calling convention by argument count.
        # n=5 → old 4-arg format: image pixel_width output viz
        # n≥6 → new 5+-arg format: image pixel_width boundary_frac output viz [label] [result]
        if n == 5:
            output_path = sys.argv[3]
            viz_path    = sys.argv[4]
            label       = ""
            result_path = None
        else:
            # sys.argv[3] = boundary_frac (reserved for future spatial constraint)
            output_path = sys.argv[4]
            viz_path    = sys.argv[5]
            label       = sys.argv[6] if n > 6 else ""
            result_path = sys.argv[7] if n > 7 else None

        img = io.imread(image_path)
        if img.ndim < 3:
            sys.stderr.write("ERROR: Expected an RGB image.\n")
            sys.exit(1)

        mask = muscularis_mask(img, pixel_width_um)

        mean_um, edt_px, skel = mean_local_thickness_um(mask, pixel_width_um)

        # Effective length: total mask area / mean thickness.
        # This is NOT the arc length of the centerline — it is the length of a
        # rectangle with the same area and same mean wall thickness.
        # Named "effective_length" to distinguish from geometric arc length.
        area_um2             = float(mask.sum()) * (pixel_width_um ** 2)
        effective_length_um  = area_um2 / mean_um if mean_um > 0 else 0.0

        # Write both values; macro reads as "mean,effective_length"
        with open(output_path, "w") as f:
            f.write("{:.4f},{:.1f}".format(mean_um, effective_length_um))

        save_visualization(img, mask, edt_px, pixel_width_um, mean_um, viz_path,
                           skel=skel)

        if result_path:
            save_result_image(img, mask, edt_px, pixel_width_um,
                              mean_um, result_path, label=label,
                              length_um=effective_length_um, skel=skel)
            sys.stdout.write("Result image saved: {}\n".format(result_path))

        sys.stdout.write(
            "Mean thickness: {:.2f} um  |  Effective length: {:.1f} um\n".format(
                mean_um, effective_length_um))

    if __name__ == "__main__":
        main()
