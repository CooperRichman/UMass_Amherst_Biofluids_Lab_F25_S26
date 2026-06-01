"""
Cooper Richman
Professor Sang Hyun Lee
University of Massachusetts Amherst Biofluids Lab
30 May 2026

The function of this script is to take in the three .csv's that are output from "Final Tracking Version.py" and graph them for analysis and figures.

CHANGE LOG
5/30/26: Reading and adding comments -C. Richman


"""

# Importing libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.widgets import Slider, Button
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os


# =========================================================
# CSV PATHS
# =========================================================

root = tk.Tk()
root.withdraw()

files = filedialog.askopenfilenames(
    title="Select ALL THREE CSV files",
    filetypes=[("CSV files", "*.csv")]
)

if len(files) != 3:
    print("Please select exactly 3 files.")
    exit()

# =========================================================
# PARSE FILES
# =========================================================
def classify_file(path):
    name = os.path.basename(path)

    if "full_paths" in name:
        return "full_paths", path
    elif "frame_summary" in name:
        return "frame_summary", path
    elif "tip_data" in name:
        return "tip_data", path
    else:
        return None, path

classified = {}

for f in files:
    key, path = classify_file(f)
    if key is None:
        print(f"Unexpected file: {path}")
        exit()
    classified[key] = path

# =========================================================
# CHECK ALL PRESENT
# =========================================================
required = ["full_paths", "frame_summary", "tip_data"]

for r in required:
    if r not in classified:
        print(f"Missing file type: {r}")
        exit()

# =========================================================
# CHECK PREFIX MATCH (XXX consistency)
# =========================================================
def get_prefix(filename):
    return filename.split("_")[0]

prefixes = [
    get_prefix(os.path.basename(classified[k]))
    for k in required
]

if len(set(prefixes)) != 1:
    print("ERROR: File prefixes do not match.")
    print("All files must start with the same XXX_ prefix.")
    exit()

# =========================================================
# ASSIGN FINAL PATHS
# =========================================================
csv1_path = classified["full_paths"]
csv2_path = classified["frame_summary"]
csv3_path = classified["tip_data"]

print("Files loaded successfully!")

# =========================================================
# LOAD CSV FILES
# =========================================================
df = pd.read_csv(csv1_path)

frame_metrics = pd.read_csv(csv2_path)

tip_metrics = pd.read_csv(csv3_path)

# =========================================================
# CLEAN DATA
# =========================================================
df = df.dropna(subset=["tip_id"])

df["tip_id"] = df["tip_id"].astype(str)

tip_metrics["tip_id"] = (
    tip_metrics["tip_id"]
    .astype(str)
)

# =========================================================
# SORT DATA
# =========================================================
df = df.sort_values(
    by=[
        "frame",
        "tip_id",
        "step_index"
    ]
)

tip_metrics = tip_metrics.sort_values(
    by=[
        "frame",
        "tip_id"
    ]
)

# =========================================================
# COORDINATE COLUMNS
# =========================================================
x_col = "pixel_x_um"
y_col = "pixel_y_um"

# =========================================================
# GET FRAMES / TIPS
# =========================================================
frames = sorted(
    df["frame"].unique()
)

tip_ids = sorted(
    df["tip_id"].unique()
)

# =========================================================
# GENERATE COLORS
# =========================================================
def generate_distinct_colors(n):

    hsv_colors = [
        (i / n, 0.85, 0.95)
        for i in range(n)
    ]

    rgb_colors = [
        mcolors.hsv_to_rgb(c)
        for c in hsv_colors
    ]

    return rgb_colors

unique_colors = generate_distinct_colors(
    len(tip_ids)
)

colors = {
    tip_id: unique_colors[i]
    for i, tip_id in enumerate(tip_ids)
}

# =========================================================
# GLOBAL LIMITS
# =========================================================
padding = 20

xmin = (
    df[x_col].min()
    - padding
)

xmax = (
    df[x_col].max()
    + padding
)

ymin = (
    df[y_col].min()
    - padding
)

ymax = (
    df[y_col].max()
    + padding
)

# =========================================================
# FIGURE
# =========================================================
fig, ax = plt.subplots(
    figsize=(16, 10)
)

# leave room for metrics
plt.subplots_adjust(
    left=0.06,
    right=0.70,
    bottom=0.18
)

# =========================================================
# FRAME SLIDER
# =========================================================
slider_ax = plt.axes(
    [0.15, 0.07, 0.40, 0.03]
)

frame_slider = Slider(
    ax=slider_ax,
    label="Frame",
    valmin=min(frames),
    valmax=max(frames),
    valinit=min(frames),
    valstep=frames
)

# =========================================================
# BUTTONS
# =========================================================
zoom_in_ax = plt.axes(
    [0.58, 0.065, 0.06, 0.04]
)

zoom_out_ax = plt.axes(
    [0.65, 0.065, 0.06, 0.04]
)

reset_ax = plt.axes(
    [0.72, 0.065, 0.09, 0.04]
)

zoom_in_button = Button(
    zoom_in_ax,
    "+"
)

zoom_out_button = Button(
    zoom_out_ax,
    "-"
)

reset_button = Button(
    reset_ax,
    "Reset"
)

# =========================================================
# STATE VARIABLES
# =========================================================
current_frame_index = 0

current_xlim = [xmin, xmax]

current_ylim = [ymax, ymin]

# =========================================================
# PAN VARIABLES
# =========================================================
is_panning = False

pan_start = None

start_xlim = None

start_ylim = None

# =========================================================
# METRICS PANEL
# =========================================================
metrics_text = fig.text(
    0.73,
    0.92,
    "",
    fontsize=10,
    family='monospace',
    verticalalignment='top'
)

# =========================================================
# DRAW FRAME
# =========================================================
def draw_frame(frame_num):

    global current_xlim
    global current_ylim

    ax.clear()

    frame_df = df[
        df["frame"] == frame_num
    ]

    # =====================================================
    # DRAW ALL PATHS
    # =====================================================
    for tip_id in tip_ids:

        tip_df = frame_df[
            frame_df["tip_id"] == tip_id
        ]

        if len(tip_df) == 0:
            continue

        x = tip_df[x_col].values
        y = tip_df[y_col].values

        # =================================================
        # DRAW PATH
        # =================================================
        ax.plot(
            x,
            y,
            '-',
            linewidth=1.2,
            color=colors[tip_id]
        )

        # =================================================
        # DRAW CURRENT TIP
        # =================================================
        ax.scatter(
            x[-1],
            y[-1],
            s=25,
            color=colors[tip_id],
            zorder=5
        )

        # =================================================
        # MOVING LABEL
        # =================================================
        ax.annotate(
            tip_id,
            (x[-1], y[-1]),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8,
            weight='bold',
            color=colors[tip_id]
        )

    # =====================================================
    # APPLY CURRENT VIEW
    # =====================================================
    ax.set_xlim(current_xlim)

    ax.set_ylim(current_ylim)

    ax.set_aspect('equal')

    # =====================================================
    # LABELS
    # =====================================================
    ax.set_title(
        f"Frame {frame_num}",
        fontsize=16
    )

    ax.set_xlabel(
        "X Distance (µm)",
        fontsize=12
    )

    ax.set_ylabel(
        "Y Distance (µm)",
        fontsize=12
    )

    # =====================================================
    # GRID
    # =====================================================
    ax.grid(
        alpha=0.3,
        linestyle='--'
    )

    # =====================================================
    # METRICS
    # =====================================================
    metrics_row = frame_metrics[
        frame_metrics["frame"] == frame_num
    ]

    if len(metrics_row) > 0:

        metrics_row = metrics_row.iloc[0]

        metrics_string = (
            f"FRAME METRICS\n"
            f"{'-'*30}\n\n"

            f"Tips:\n"
            f"{int(metrics_row['n_tips'])}\n\n"

            f"Faux Tips:\n"
            f"{int(metrics_row['n_faux_tips'])}\n\n"

            f"Junctions:\n"
            f"{int(metrics_row['n_junctions_normal'])}\n\n"

            f"Crossover Junctions:\n"
            f"{int(metrics_row['n_junctions_crossover'])}\n\n"

            f"Cancelled Junctions:\n"
            f"{int(metrics_row['n_junctions_removed'])}\n\n"

            f"Skeleton Area:\n"
            f"{metrics_row['skeleton_area_um2']:.2f} µm²"
        )

        metrics_text.set_text(
            metrics_string
        )

    else:

        metrics_text.set_text(
            "No metrics found"
        )

    plt.draw()

# =========================================================
# SLIDER UPDATE
# =========================================================
def update(val):

    global current_frame_index

    frame_num = int(
        frame_slider.val
    )

    current_frame_index = (
        frames.index(frame_num)
    )

    draw_frame(frame_num)

frame_slider.on_changed(update)

# =========================================================
# KEYBOARD CONTROLS
# =========================================================
def on_key(event):

    global current_frame_index

    if event.key == 'right':

        current_frame_index = min(
            current_frame_index + 1,
            len(frames) - 1
        )

    elif event.key == 'left':

        current_frame_index = max(
            current_frame_index - 1,
            0
        )

    else:
        return

    new_frame = frames[
        current_frame_index
    ]

    frame_slider.set_val(new_frame)

fig.canvas.mpl_connect(
    'key_press_event',
    on_key
)

# =========================================================
# ZOOM FUNCTION
# =========================================================
def zoom(factor):

    global current_xlim
    global current_ylim

    x_center = (
        current_xlim[0]
        + current_xlim[1]
    ) / 2

    y_center = (
        current_ylim[0]
        + current_ylim[1]
    ) / 2

    x_width = (
        current_xlim[1]
        - current_xlim[0]
    ) * factor

    y_width = (
        current_ylim[0]
        - current_ylim[1]
    ) * factor

    current_xlim = [
        x_center - x_width / 2,
        x_center + x_width / 2
    ]

    current_ylim = [
        y_center + y_width / 2,
        y_center - y_width / 2
    ]

    draw_frame(
        int(frame_slider.val)
    )

# =========================================================
# BUTTON CALLBACKS
# =========================================================
def zoom_in(event):

    zoom(0.8)

def zoom_out(event):

    zoom(1.25)

def reset_view(event):

    global current_xlim
    global current_ylim

    current_xlim = [
        xmin,
        xmax
    ]

    current_ylim = [
        ymax,
        ymin
    ]

    draw_frame(
        int(frame_slider.val)
    )

zoom_in_button.on_clicked(
    zoom_in
)

zoom_out_button.on_clicked(
    zoom_out
)

reset_button.on_clicked(
    reset_view
)

# =========================================================
# PAN FUNCTIONS
# =========================================================
def on_mouse_press(event):

    global is_panning
    global pan_start
    global start_xlim
    global start_ylim

    if event.inaxes != ax:
        return

    if event.button != 1:
        return

    is_panning = True

    pan_start = (
        event.xdata,
        event.ydata
    )

    start_xlim = ax.get_xlim()

    start_ylim = ax.get_ylim()

def on_mouse_release(event):

    global is_panning

    is_panning = False

def on_mouse_move(event):

    global current_xlim
    global current_ylim

    if not is_panning:
        return

    if event.inaxes != ax:
        return

    if event.xdata is None:
        return

    if event.ydata is None:
        return

    dx = (
        event.xdata
        - pan_start[0]
    )

    dy = (
        event.ydata
        - pan_start[1]
    )

    current_xlim = [
        start_xlim[0] - dx,
        start_xlim[1] - dx
    ]

    current_ylim = [
        start_ylim[0] - dy,
        start_ylim[1] - dy
    ]

    ax.set_xlim(current_xlim)

    ax.set_ylim(current_ylim)

    fig.canvas.draw_idle()

# =========================================================
# CONNECT PAN EVENTS
# =========================================================
fig.canvas.mpl_connect(
    'button_press_event',
    on_mouse_press
)

fig.canvas.mpl_connect(
    'button_release_event',
    on_mouse_release
)

fig.canvas.mpl_connect(
    'motion_notify_event',
    on_mouse_move
)

# =========================================================
# INITIAL DRAW
# =========================================================
draw_frame(frames[0])

plt.show()
