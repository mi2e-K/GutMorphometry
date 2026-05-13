#!/usr/bin/env python3
"""
tune_muscularis.py  --  Interactive threshold tuning for muscularis detection.

Usage (called automatically from Fiji, or run manually):
    python tune_muscularis.py <image_tif> [<pixel_width_um>]

Opens an interactive matplotlib window showing:
  Left   : original polygon-masked crop  (click here to pick reference colours)
  Centre : detected muscle mask (green overlay) -- updates in real time
  Right  : E-H feature map with detected boundary (yellow contour)

Controls:
  "Pick Pink"   -- click this button, then click a pink (muscle) pixel in Panel 1
  "Pick Purple" -- click this button, then click a purple (crypt) pixel in Panel 1
  E-H factor slider -- fine-tune the boundary after reference colours are set

Algorithm:
  When BOTH reference colours are set (recommended):
      Nearest-centroid classification in optical-density (OD) space.
      Each tissue pixel is assigned to whichever reference (pink or purple)
      it is closer to in 3-D OD space.  This is independent of whether the
      target tissue has higher or lower E-H than the background -- it only
      requires that the two reference pixels have visually distinct colours.
      The slider shifts the decision boundary (positive = stricter / more pink).

  When no references are set (fallback):
      Single Otsu threshold on the E-H feature (E minus H) within tissue.

Click "Save & Close" to write values to musc_params.json.
musc_thickness.py reads this file automatically on every measurement.
"""

# ------------------------------------------------------------------
# Jython guard
# ------------------------------------------------------------------
import sys as _sys

if _sys.version_info[0] >= 3:

    import sys
    import json
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button
    from skimage import io, filters
    from skimage.measure import label
    from skimage.morphology import remove_small_objects, binary_closing, disk
    from scipy.ndimage import binary_fill_holes, gaussian_filter

    # ------------------------------------------------------------------
    # H&E deconvolution (identical to musc_thickness.py)
    # ------------------------------------------------------------------
    _HE_MATRIX = np.array([
        [0.6500286, 0.7040924, 0.2860126],
        [0.0704327, 0.9778137, 0.1972729],
        [0.7110272, 0.4232444, 0.5603568],
    ], dtype=np.float64)
    _HE_INV = np.linalg.inv(_HE_MATRIX)

    def _deconvolve_he(rgb):
        rgb_f  = np.clip(rgb[:, :, :3].astype(np.float64), 1, 255)
        od     = -np.log10(rgb_f / 255.0)
        stains = np.dot(od.reshape(-1, 3), _HE_INV.T)
        stains = np.clip(stains, 0, None)
        h = stains[:, 0].reshape(rgb.shape[:2]).astype(np.float32)
        e = stains[:, 1].reshape(rgb.shape[:2]).astype(np.float32)
        return h, e

    def _tissue_mask(rgb):
        rgb_f = np.clip(rgb[:, :, :3].astype(np.float64), 1, 255)
        od    = -np.log10(rgb_f / 255.0)
        return od.sum(axis=2) >= 0.15

    def _rgb_to_od(rgb_pixel):
        """3-element OD vector for a single RGB pixel (list or array)."""
        px = np.clip(np.array(rgb_pixel[:3], dtype=np.float64), 1.0, 255.0)
        return -np.log10(px / 255.0)

    # ------------------------------------------------------------------
    # Mask computation
    # ------------------------------------------------------------------
    def _compute_mask(rgb, pixel_width_um,
                      feat_base, feat_scale, e_factor,
                      pink_rgb=None, purple_rgb=None, bg_rgb=None):
        """
        Return a boolean mask selecting muscularis pixels.

        Pipeline is identical to musc_thickness.py so the calibration preview
        shows exactly the same boundary that will appear in the result image.

        Mode A3 – 3-class reference mode (pink, purple, AND background set):
            Nearest-centroid in OD space with 3 classes.
        Mode A2 – 2-class reference mode (pink and purple only):
            Nearest-centroid between pink and purple.
        Mode B – Otsu fallback (no references):
            E-H > (feat_base + (e_factor-1) * feat_scale)
        """
        h, e    = _deconvolve_he(rgb)
        tissue  = _tissue_mask(rgb)

        if pink_rgb is not None and purple_rgb is not None:
            # --- Mode A: nearest-centroid in OD space ---
            rgb_f     = np.clip(rgb[:, :, :3].astype(np.float64), 1.0, 255.0)
            od        = -np.log10(rgb_f / 255.0)
            pink_od   = _rgb_to_od(pink_rgb)
            purple_od = _rgb_to_od(purple_rgb)
            d_pink    = np.sqrt(((od - pink_od)   ** 2).sum(axis=2))
            d_purple  = np.sqrt(((od - purple_od) ** 2).sum(axis=2))
            margin = (e_factor - 1.0) * feat_scale
            mask = tissue & ((d_purple - d_pink) > margin)
            if bg_rgb is not None:
                bg_od  = _rgb_to_od(bg_rgb)
                d_bg   = np.sqrt(((od - bg_od) ** 2).sum(axis=2))
                mask = mask & ((d_bg - d_pink) > margin)
        else:
            # --- Mode B: E-H Otsu fallback ---
            feature     = e - h
            feat_thresh = feat_base + (e_factor - 1.0) * feat_scale
            mask = tissue & (feature > feat_thresh)

        # --- Full morphological pipeline (identical to musc_thickness.py) ---
        r_bridge = max(3, min(20, int(round(8.0 / pixel_width_um))))
        mask = binary_closing(mask, disk(r_bridge))
        mask = binary_fill_holes(mask)
        mask = binary_closing(mask, disk(3))
        min_px = max(64, int(round(400.0 / (pixel_width_um ** 2))))
        mask = remove_small_objects(mask, min_size=min_px)
        labeled = label(mask)
        if labeled.max() > 0:
            sizes    = np.bincount(labeled.ravel())
            sizes[0] = 0
            mask     = labeled == sizes.argmax()
        mask = binary_fill_holes(mask)
        _gs  = gaussian_filter(mask.astype(np.float32), sigma=2.0)
        mask = binary_fill_holes(_gs >= 0.5)

        return mask.astype(bool)

    def _make_overlay(img_rgb, mask):
        ov = img_rgb[:, :, :3].astype(np.int16)
        ov[mask, 1] = np.clip(ov[mask, 1] + 50, 0, 255)
        ov[mask, 0] = np.clip(ov[mask, 0] - 20, 0, 255)
        ov[mask, 2] = np.clip(ov[mask, 2] - 20, 0, 255)
        return ov.astype(np.uint8)

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------
    def main():
        if len(sys.argv) < 2 or len(sys.argv) > 3:
            sys.stderr.write(
                "Usage: python tune_muscularis.py <image_tif> [<pixel_width_um>]\n"
            )
            sys.exit(1)

        image_path     = sys.argv[1]
        pixel_width_um = float(sys.argv[2]) if len(sys.argv) == 3 else 1.0

        img = io.imread(image_path)
        if img.ndim < 3:
            sys.stderr.write("ERROR: Expected RGB image.\n")
            sys.exit(1)

        # ----------------------------------------------------------
        # Image-level statistics (fixed for this session)
        # ----------------------------------------------------------
        h_ch, e_ch  = _deconvolve_he(img)
        tissue      = _tissue_mask(img)
        feature     = e_ch - h_ch
        feat_vals   = feature[tissue]

        otsu_base  = float(filters.threshold_otsu(feat_vals)) if feat_vals.size > 100 else 0.0
        if feat_vals.size > 100:
            feat_iqr   = float(np.percentile(feat_vals, 75) - np.percentile(feat_vals, 25))
            feat_scale_eh = max(0.05, feat_iqr * 0.30)
        else:
            feat_scale_eh = 0.05
        abs98 = float(np.percentile(np.abs(feature), 98))

        # ----------------------------------------------------------
        # Params file
        # ----------------------------------------------------------
        _autorun    = os.path.dirname(os.path.abspath(__file__))
        _fiji_root  = os.path.dirname(os.path.dirname(os.path.dirname(_autorun)))
        _config_dir = os.path.join(_fiji_root, "GutMorphometry_config")
        os.makedirs(_config_dir, exist_ok=True)
        params_path = os.path.join(_config_dir, "musc_params.json")

        for _bad in [os.path.join(_autorun, "musc_params.json"),
                     os.path.join(_autorun, "GutMorphometry_config", "musc_params.json")]:
            try:
                if os.path.exists(_bad):
                    os.remove(_bad)
            except Exception:
                pass
        try:
            os.rmdir(os.path.join(_autorun, "GutMorphometry_config"))
        except Exception:
            pass

        # Load previously saved parameters
        e_init       = 1.0
        pink_rgb_0   = None
        purple_rgb_0 = None
        bg_rgb_0     = None
        if os.path.exists(params_path):
            try:
                with open(params_path) as f:
                    p = json.load(f)
                e_init       = float(p.get("e_factor", 1.0))
                pink_rgb_0   = p.get("pink_rgb",   None)
                purple_rgb_0 = p.get("purple_rgb", None)
                bg_rgb_0     = p.get("bg_rgb",     None)
            except Exception:
                pass

        # ----------------------------------------------------------
        # Mutable state shared across closures
        # ----------------------------------------------------------
        # ref_rgb[0] = pink, ref_rgb[1] = purple, ref_rgb[2] = background
        ref_rgb   = [pink_rgb_0, purple_rgb_0, bg_rgb_0]
        pick_mode = [None]   # "pink" | "purple" | "bg" | None

        def _has_refs():
            return ref_rgb[0] is not None and ref_rgb[1] is not None

        def _od_scale():
            """Slider scale in OD-distance units: 20% of pink-purple separation."""
            if not _has_refs():
                return feat_scale_eh
            d = np.linalg.norm(_rgb_to_od(ref_rgb[0]) - _rgb_to_od(ref_rgb[1]))
            return max(0.02, d * 0.20)

        state = {"e_factor": e_init}

        # ----------------------------------------------------------
        # Debug log
        # ----------------------------------------------------------
        import tempfile as _tmpmod
        _log_path = os.path.join(_tmpmod.gettempdir(), "tune_muscularis_debug.txt")
        _log_f    = open(_log_path, "w", buffering=1, encoding="utf-8")

        def _log(msg):
            print(msg)
            _log_f.write(msg + "\n")

        _log("tune_muscularis started")
        _log("  image      : {}".format(image_path))
        _log("  pixel_um   : {}".format(pixel_width_um))
        _log("  otsu_base  : {:.4f}".format(otsu_base))
        _log("  feat_scale : {:.4f}  (E-H Otsu mode)".format(feat_scale_eh))
        _log("  pink_rgb_0  : {}".format(pink_rgb_0))
        _log("  purple_rgb_0: {}".format(purple_rgb_0))
        _log("  bg_rgb_0    : {}".format(bg_rgb_0))
        _log("  tissue px  : {}".format(int(tissue.sum())))
        _log("---")

        # ----------------------------------------------------------
        # Initial mask
        # ----------------------------------------------------------
        mask0 = _compute_mask(img, pixel_width_um,
                              otsu_base, feat_scale_eh, e_init,
                              pink_rgb=ref_rgb[0], purple_rgb=ref_rgb[1],
                              bg_rgb=ref_rgb[2])
        _log("  initial mask px: {}".format(int(mask0.sum())))

        # ----------------------------------------------------------
        # Figure layout
        # ----------------------------------------------------------
        fig = plt.figure(figsize=(15, 8))
        fig.patch.set_facecolor("#1e1e1e")

        ax1 = fig.add_axes([0.01, 0.30, 0.31, 0.65])
        ax2 = fig.add_axes([0.34, 0.30, 0.31, 0.65])
        ax3 = fig.add_axes([0.67, 0.30, 0.31, 0.65])

        ax_se = fig.add_axes([0.12, 0.20, 0.76, 0.04])
        ax_sh = fig.add_axes([0.12, 0.13, 0.76, 0.04])

        ax_btn_pink   = fig.add_axes([0.01, 0.06, 0.13, 0.05])
        ax_btn_purple = fig.add_axes([0.16, 0.06, 0.13, 0.05])
        ax_btn_bg     = fig.add_axes([0.31, 0.06, 0.13, 0.05])
        ax_btn_clear  = fig.add_axes([0.46, 0.06, 0.08, 0.05])
        ax_ref_status = fig.add_axes([0.56, 0.06, 0.42, 0.05])

        ax_btn    = fig.add_axes([0.78, 0.005, 0.20, 0.05])
        ax_status = fig.add_axes([0.01, 0.005, 0.75, 0.05])

        for ax in [ax1, ax2, ax3, ax_ref_status, ax_status]:
            ax.axis("off")
            ax.set_facecolor("#1e1e1e")

        # Markers on ax1 showing where references were picked (pink, purple, bg)
        ref_markers = [None, None, None]

        # ----------------------------------------------------------
        # Panels
        # ----------------------------------------------------------
        ax1.imshow(img[:, :, :3])
        ax1.set_title("Panel 1 - Original  [click here after pressing a Pick button]",
                      fontsize=9, color="white")

        ov_img = ax2.imshow(_make_overlay(img, mask0))
        ax2.set_title("Detected mask  [preview - largest component only]",
                      fontsize=9, color="white")

        ax3.imshow(feature, cmap="RdBu_r", vmin=-abs98, vmax=abs98)
        contour_ref = [None]
        if mask0.any():
            contour_ref[0] = ax3.contour(mask0, levels=[0.5],
                                         colors="yellow", linewidths=1.5)
        ax3.set_title("E-H map  (red=eosin-dominant, blue=H-dominant)\nyellow = boundary",
                      fontsize=9, color="white")

        def _update_ax3_label():
            if _has_refs():
                p_od  = _rgb_to_od(ref_rgb[0])
                q_od  = _rgb_to_od(ref_rgb[1])
                sep   = float(np.linalg.norm(p_od - q_od))
                ax3.set_xlabel(
                    "Mode: OD nearest-centroid  |  OD separation = {:.3f}".format(sep),
                    color="#cccccc", fontsize=7)
            else:
                ax3.set_xlabel(
                    "Mode: E-H Otsu  |  baseline = {:.3f}".format(otsu_base),
                    color="#cccccc", fontsize=7)

        _update_ax3_label()

        # ----------------------------------------------------------
        # Slider
        # ----------------------------------------------------------
        slider_e = Slider(ax_se,
                          "Factor   (up = stricter/pink-only   down = more inclusive)",
                          -2.00, 3.00, valinit=e_init, valstep=0.05,
                          color="#e07b54")
        slider_e.label.set_color("white")
        slider_e.valtext.set_color("white")

        ax_sh.axis("off")
        ax_sh.set_facecolor("#1e1e1e")
        info_txt = ax_sh.text(
            0.5, 0.5,
            "OD mode: margin = (factor-1) x OD_scale  |  "
            "Otsu fallback: thresh = otsu + (factor-1) x scale  |  "
            "Use toolbar zoom to enlarge Panel 1 before picking reference pixels",
            transform=ax_sh.transAxes, fontsize=7.5,
            color="#999999", ha="center", va="center")

        # ----------------------------------------------------------
        # Reference status
        # ----------------------------------------------------------
        def _ref_status_str():
            p_rgb  = ref_rgb[0]
            q_rgb  = ref_rgb[1]
            b_rgb  = ref_rgb[2]
            if p_rgb is None and q_rgb is None:
                return "[Otsu mode - press Pick Pink then click on a pink pixel in Panel 1]", "#aaaaaa"
            if p_rgb is not None and q_rgb is not None:
                p_od  = _rgb_to_od(p_rgb)
                q_od  = _rgb_to_od(q_rgb)
                sep   = float(np.linalg.norm(p_od - q_od))
                if sep < 0.05:
                    col  = "#ff4444"
                    note = "  !! sep={:.3f} too small - re-pick !!".format(sep)
                elif sep < 0.15:
                    col  = "#ffaa44"
                    note = "  (sep={:.3f} marginal)".format(sep)
                else:
                    col  = "#44dd44"
                    note = "  (sep={:.3f} OK)".format(sep)
                bg_str = "  BG:{}".format(b_rgb) if b_rgb is not None else "  BG:(not set)"
                mode_str = "[3-class OD]" if b_rgb is not None else "[2-class OD]"
                return ("P:{} Q:{}{}{}  {}"
                        .format(p_rgb, q_rgb, bg_str, note, mode_str)), col
            elif p_rgb is not None:
                return "Pink set:{}   Purple:(not set)  [Otsu mode]".format(p_rgb), "#aaaaaa"
            else:
                return "Pink:(not set)   Purple set:{}  [Otsu mode]".format(q_rgb), "#aaaaaa"

        _rs0, _rc0 = _ref_status_str()
        ref_status_txt = ax_ref_status.text(
            0.0, 0.5, _rs0,
            transform=ax_ref_status.transAxes,
            fontsize=7.5, color=_rc0, va="center")

        btn_pink = Button(ax_btn_pink, "1. Pick Pink (muscle)",
                          color="#7a3a3a", hovercolor="#9a4a4a")
        btn_pink.label.set_color("white")
        btn_pink.label.set_fontsize(8)

        btn_purple = Button(ax_btn_purple, "2. Pick Purple (crypt)",
                            color="#3a2a5a", hovercolor="#4a3a7a")
        btn_purple.label.set_color("white")
        btn_purple.label.set_fontsize(8)

        btn_bg = Button(ax_btn_bg, "3. Pick Background",
                        color="#2a4a3a", hovercolor="#3a6a4a")
        btn_bg.label.set_color("white")
        btn_bg.label.set_fontsize(8)

        btn_clear = Button(ax_btn_clear, "Clear refs",
                           color="#4a4a2a", hovercolor="#6a6a3a")
        btn_clear.label.set_color("white")
        btn_clear.label.set_fontsize(8)

        # ----------------------------------------------------------
        # Status bar
        # ----------------------------------------------------------
        def _status_str(ef):
            if _has_refs():
                sc = _od_scale()
                margin = (ef - 1.0) * sc
                return ("Factor: {:.2f}   OD margin: {:.3f}   scale: {:.3f}"
                        "   [OD nearest-centroid mode]"
                        ).format(ef, margin, sc)
            else:
                thresh = otsu_base + (ef - 1.0) * feat_scale_eh
                return ("Factor: {:.2f}   Otsu baseline: {:.3f}   "
                        "scale: {:.3f}   threshold: {:.3f}   [Otsu mode]"
                        ).format(ef, otsu_base, feat_scale_eh, thresh)

        status_txt = ax_status.text(
            0.0, 0.5, _status_str(e_init),
            transform=ax_status.transAxes,
            fontsize=8, color="lightgrey", va="center")

        btn_save = Button(ax_btn, "Save & Close",
                          color="#2d6a2d", hovercolor="#3d8a3d")
        btn_save.label.set_color("white")

        # ----------------------------------------------------------
        # Core refresh
        # ----------------------------------------------------------
        def _refresh(ef):
            sc = _od_scale()
            new_mask = _compute_mask(img, pixel_width_um,
                                     otsu_base, sc, ef,
                                     pink_rgb=ref_rgb[0], purple_rgb=ref_rgb[1],
                                     bg_rgb=ref_rgb[2])
            _log("_refresh  ef={:.2f}  mode={}  px={}".format(
                ef, "OD" if _has_refs() else "Otsu", int(new_mask.sum())))

            ov_img.set_data(_make_overlay(img, new_mask))

            try:
                if contour_ref[0] is not None:
                    try:
                        contour_ref[0].remove()
                    except AttributeError:
                        for coll in contour_ref[0].collections:
                            coll.remove()
                    contour_ref[0] = None
                if new_mask.any():
                    contour_ref[0] = ax3.contour(new_mask, levels=[0.5],
                                                  colors="yellow", linewidths=1.5)
            except Exception as _exc:
                _log("  contour error: {}".format(_exc))

            status_txt.set_text(_status_str(ef))
            _rs, _rc = _ref_status_str()
            ref_status_txt.set_text(_rs)
            ref_status_txt.set_color(_rc)
            _update_ax3_label()
            fig.canvas.draw_idle()

        # ----------------------------------------------------------
        # Slider callback
        # ----------------------------------------------------------
        def _on_slider(_):
            ef = slider_e.val
            state["e_factor"] = ef
            _refresh(ef)

        # ----------------------------------------------------------
        # Reference colour picking
        # ----------------------------------------------------------
        def _activate_pink(_):
            pick_mode[0] = "pink"
            ax1.set_title(">>> Click on a PINK (muscle) pixel in Panel 1 <<<\n"
                          "Tip: use toolbar ZOOM first to enlarge the pink area",
                          fontsize=8, color="#ff9999")
            fig.canvas.draw_idle()
            _log("Pick mode: pink")

        def _activate_purple(_):
            pick_mode[0] = "purple"
            ax1.set_title(">>> Click on a PURPLE (crypt) pixel in Panel 1 <<<",
                          fontsize=9, color="#cc99ff")
            fig.canvas.draw_idle()
            _log("Pick mode: purple")

        def _activate_bg(_):
            pick_mode[0] = "bg"
            ax1.set_title(">>> Click on a BACKGROUND pixel in Panel 1 <<<\n"
                          "Pick from outside the tissue but inside the polygon",
                          fontsize=8, color="#66ffaa")
            fig.canvas.draw_idle()
            _log("Pick mode: bg")

        def _on_click(event):
            if pick_mode[0] is None:
                return
            if event.inaxes is not ax1:
                return
            col = int(round(event.xdata)) if event.xdata is not None else -1
            row = int(round(event.ydata)) if event.ydata is not None else -1
            if not (0 <= row < img.shape[0] and 0 <= col < img.shape[1]):
                return

            px_rgb = img[row, col, :3].tolist()   # [R, G, B]
            mode   = pick_mode[0]
            slot   = 0 if mode == "pink" else (1 if mode == "purple" else 2)
            pick_mode[0] = None

            ref_rgb[slot] = px_rgb
            _log("{} picked at ({},{}) RGB={} OD={}".format(
                mode, row, col, px_rgb,
                [round(v, 3) for v in _rgb_to_od(px_rgb).tolist()]))

            # Draw marker on ax1
            mcolor = "#ff6666" if slot == 0 else ("#bb88ff" if slot == 1 else "#66ffaa")
            if ref_markers[slot] is not None:
                try:
                    ref_markers[slot].remove()
                except Exception:
                    pass
            ref_markers[slot] = ax1.scatter(
                [col], [row], s=120, c=mcolor,
                marker="x", linewidths=2, zorder=10)

            ax1.set_title("Panel 1 - Original  [click here after pressing a Pick button]",
                          fontsize=9, color="white")

            # Log OD separation if both set
            if _has_refs():
                p_od = _rgb_to_od(ref_rgb[0])
                q_od = _rgb_to_od(ref_rgb[1])
                sep  = float(np.linalg.norm(p_od - q_od))
                _log("  OD separation = {:.4f}  OD_scale = {:.4f}".format(
                    sep, _od_scale()))
                if sep < 0.05:
                    _log("  WARNING: very small OD separation - try picking more distinct pixels")

            _refresh(state["e_factor"])

        # ----------------------------------------------------------
        # Clear references
        # ----------------------------------------------------------
        def _clear_refs(_):
            ref_rgb[0] = None
            ref_rgb[1] = None
            ref_rgb[2] = None
            for i, m in enumerate(ref_markers):
                if m is not None:
                    try:
                        m.remove()
                    except Exception:
                        pass
                ref_markers[i] = None
            ax1.set_title("Panel 1 - Original  [click here after pressing a Pick button]",
                          fontsize=9, color="white")
            _log("References cleared - back to Otsu mode")
            _refresh(state["e_factor"])

        # ----------------------------------------------------------
        # Save
        # ----------------------------------------------------------
        def _save(_):
            params = {"e_factor": state["e_factor"]}
            if ref_rgb[0] is not None:
                params["pink_rgb"] = ref_rgb[0]
            if ref_rgb[1] is not None:
                params["purple_rgb"] = ref_rgb[1]
            if ref_rgb[2] is not None:
                params["bg_rgb"] = ref_rgb[2]
            with open(params_path, "w") as f:
                json.dump(params, f, indent=2)
            _log("Saved: {}".format(params))
            _log_f.close()
            plt.close(fig)

        # ----------------------------------------------------------
        # Wire up events
        # ----------------------------------------------------------
        slider_e.on_changed(_on_slider)
        btn_pink.on_clicked(_activate_pink)
        btn_purple.on_clicked(_activate_purple)
        btn_bg.on_clicked(_activate_bg)
        btn_clear.on_clicked(_clear_refs)
        btn_save.on_clicked(_save)
        fig.canvas.mpl_connect("button_press_event", _on_click)

        fig.suptitle(
            "Muscularis Threshold Tuning   |   {}".format(
                os.path.basename(image_path)),
            fontsize=10, color="white")

        plt.show()

    if __name__ == "__main__":
        main()
