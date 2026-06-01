"""
Cooper Richman
Professor Sang Hyun Lee
University of Massachusetts Amherst Biofluids Lab
29 May 2026

***** PLEASE SEE RELEVANT DOCUMENTATION IN THE GITHUB TO FULLY UNDERSTAND THE ARCHITECTURE OF THIS PROGRAM *****
----------------------------------------------------------------------------------------------------------------


================================================================================
Mycelial Path Viewer & Export Tool
================================================================================

WHAT THIS SCRIPT DOES
--------------------------------------------------------------------------------
This script loads a time-lapse image stack of a skeleton in the form of a .tif file.
It detects:

  • TIPS                  — the free ends of branches
  • JUNCTIONS/CROSSOVERS  — places where branches meet / cross
  • ORIGIN                — the single "root" point of the whole network

It then tries to figure out which tip is growing from which junction by
finding the shortest skeleton path between them.

You review the junctions in PHASE 1 (mark crossovers, remove junctions etc.).
Then in PHASE 2 you see the pairings drawn as coloured paths.  You can
manually override any pairing that the algorithm got wrong.

Results are exported to three CSV files for downstream analysis.

================================================================================


CHANGE LOG
5/13/26: Polishing the first version of this script to be handed off to UMass Biofluids Lab -C. Richman
5/14/26: Making some more small edits to comments and layout -C. Richman
5/29/26: Ensuring that things are clear and understandable -C. Richman

"""


# ==============================================================================
# ███  SECTION 1 — IMPORTS  ███
# ==============================================================================
# There are lots of extra tools that are necessary to make a program like this
# possible. They are imported here. Note that most of these do not come with
# a default Python install, so they have been installed using pip. See more
# info here: https://pypi.org/project/pip/
# Python comes with a huge "standard library" of pre-written tools.  We also
# install extra packages (numpy, matplotlib, skimage, scipy) with pip.
# "import X" means "load this tool so we can use it later as X.something".
# "from X import Y" means "load just the Y part from tool X".
# ==============================================================================

import io           # Tools for reading/writing data in memory (used for clipboard saving of images)
import math         # Basic maths: sqrt, atan2, hypot, etc.
import os           # Operating-system tools: file paths, directory listing
import csv          # Read and write comma-separated-value (.csv) files
import heapq        # A "priority queue" — needed by Dijkstra's pathfinding
import numpy as np  # NumPy: fast arrays and matrix maths (images are big arrays)

import matplotlib                              # The main plotting library
matplotlib.use("TkAgg")                        # Use the Tkinter backend (desktop window)
import matplotlib.pyplot as plt                # The high-level plotting API
import matplotlib.backends.backend_tkagg       # Connects matplotlib to the Tkinter window

import tkinter as tk                           # Tkinter: Python's built-in GUI toolkit
from tkinter import scrolledtext, simpledialog, filedialog, messagebox
# scrolledtext  — a text box with a scrollbar
# simpledialog  — pop-up "ask a question" dialogs
# filedialog    — file-open / file-save dialogs
# messagebox    — pop-up info / warning / error boxes

from skimage.io import imread                  # Read image files (TIFF stacks, etc.)
from scipy.ndimage import convolve, binary_dilation
# convolve        — slide a small filter over an image (used to count neighbours)
# binary_dilation — "inflate" a binary mask by a few pixels (cosmetic rendering)

from collections import defaultdict
# defaultdict — a dictionary that creates a default value for missing keys
# (used to build the skeleton graph without crashing on missing edges)

# ==============================================================================
# ███  FILEPATH SELECTION  ███
# ==============================================================================
# You can either have a small popup to select your file or paste your filepath
# directly. Make sure to comment out the method that you are not using. This
# means that it should be EITHER loading the popup OR using your pasted
# filepath. It CANNOT be both. Use the # to comment out lines.
# ==============================================================================


# If you would like a file selection popup, leave this uncommented
# ------------------------------------------------------------------------------
# To get the filepath we must start a Tkinter window and then close it right after
root = tk.Tk()
root.withdraw()

# Prompt the user to pick their file, note that it can be .tiff or .tif
SKELETON_TIFF = filedialog.askopenfilename(
    title="Select skeletonized .tiff file",
    filetypes=[("TIFF files", "*.tif *.tiff")]
)

# If there was nothing selected, close the whole program
if not SKELETON_TIFF:
    print("No file selected.")
    exit()

# Fully and cleanly close this Tkinter window
root.destroy()
# ------------------------------------------------------------------------------

# If you would like the direct filepath paste, uncomment this AND comment
# out the lines above.
# ------------------------------------------------------------------------------

#SKELETON_TIFF = r"C:\Users\Documents\my_file.tif"

# ------------------------------------------------------------------------------

# ==============================================================================
# ███  SECTION 2 — COLOUR PALETTE  ███
# ==============================================================================
# These are HTML hex colour codes.
# Change them here and every part of the script that uses them will update.
# ==============================================================================

PALETTE = dict(
    JUNCTION_USED   = "#e040fb",   # Purple  — junction that has a tip paired to it
    JUNCTION_UNUSED = "#e040fb",   # Purple  — junction with no tip (same colour here, could maybe be changed)
    CROSSOVER       = "#ff5e00",   # Orange  — junction marked as a crossover
    TIP             = "#384ce4",   # Blue    — free tip of a branch
    ORIGIN          = "#ff1744",   # Red     — the root / origin point
)


# ==============================================================================
# ███  SECTION 3 — CONFIGURATION CONSTANTS  ███
# ==============================================================================
# These numbers control how the algorithm and GUI behave.
# You can tweak them without understanding the rest of the code.
# ==============================================================================









# How many micrometres (µm) each pixel represents in real space.
# Used to convert pixel distances into physical distances in the CSV exports.
PIXEL_SIZE_UM = 1.73

# ── Display defaults ──────────────────────────────────────────────────────────
invert_view         = True    # True = white background (easier to see skeleton)
show_overlay_master = True    # True = show paths, labels, markers on startup

# ── Anchor / origin detection ────────────────────────────────────────────────
# The "anchor" is the stable reference point used as the network origin.
# We find it by looking for pixels that are lit in at least this fraction
# of all frames — i.e. a pixel that is always part of the skeleton.
ANCHOR_MIN_FRAME_FRACTION = 0.50   # 50 % of frames must contain the pixel

# When snapping the anchor to the nearest skeleton pixel, search within
# this many pixels.
ANCHOR_SNAP_RADIUS = 20   # pixels

# ── Hover / click detection ───────────────────────────────────────────────────
# How close (in pixels) does the mouse have to be to "hover" over a path?
HOVER_RADIUS_PX = 6
# How close to be considered hovering over a junction centroid?
JUNCTION_HOVER_PX = 14
# How close to be considered "clicking" a junction?
CLICK_ZONE_RADIUS_PX = 15
# How close to be considered "clicking" a tip?
CLICK_TIP_RADIUS_PX = 10

# ── Junction ID tracking ──────────────────────────────────────────────────────
# When matching junctions across frames we look for junctions whose centres
# are within this many pixels of each other. Note that this can be changed 
# with a slider in the GUI as well.
DEFAULT_MATCH_RADIUS = 20   # pixels

# ── Zone state constants ──────────────────────────────────────────────────────
# Every junction has a "state" — an integer that says what kind it is.
# 0 = normal (use it in pairing)
# 1 = crossover ×  (ignore in pairing)
# 2 = crossover ∥  (ignore in pairing)
# 3 = removed      (ignore completely)
ZONE_STATE_NORMAL   = 0
ZONE_STATE_CROSS    = 1
ZONE_STATE_PARALLEL = 2
ZONE_STATE_REMOVED  = 3

# A special "zone ID" used to mark a pairing that goes to the ORIGIN rather
# than to a real junction.  We use -1 because real zone IDs start at 0.
ORIGIN_ZONE_ID = -1

# ==============================================================================
# ███  SECTION 4 — PHASE / OVERRIDE / FAUX-TIP STATE VARIABLES  ███
# ==============================================================================
# "Global state" — these variables are shared across the whole script and
# change as the user interacts with the GUI.
# ==============================================================================

# Which phase of the workflow are we in?
PHASE_REVIEW  = "review"    # Phase 1 — review and mark junctions
PHASE_ANALYZE = "analyze"   # Phase 2 — see/export tip–junction pairings

# Override mode — used when the user wants to manually re-route a tip.
OVERRIDE_MODE_NONE = "none"        # Not in override mode
OVERRIDE_MODE_TIP  = "tip_selected"  # A tip has been selected, waiting for target

# These hold the current override state.
override_mode = OVERRIDE_MODE_NONE
override_tip  = None   # The tip (pixel coordinate) selected for overriding

# Faux tip mode — used when the user wants to create a ghost copy of a tip.
FAUX_MODE_NONE     = "none"            # Not in faux-tip mode
FAUX_MODE_WAITING  = "waiting_tip"     # Waiting for the user to click a real tip

faux_mode = FAUX_MODE_NONE   # Current faux-tip mode state


# ==============================================================================
# ███  SECTION 5 — SMALL HELPER FUNCTIONS  ███
# ==============================================================================
# Tiny reusable utilities used throughout the rest of the script.
# ==============================================================================

def ipt(p):
    """
    Convert a point (which might be floats from numpy) to a plain Python tuple
    of integers.  We need integer tuples because we use them as dictionary keys
    and dictionary keys must be hashable (plain integers are, numpy floats are not).

    Example:  ipt((1.9, 3.0))  ->  (1, 3)
    """
    return (int(p[0]), int(p[1]))


def euclid(a, b):
    """
    Euclidean (straight-line) distance between two 2-D points a and b.
    Uses the Pythagorean theorem:  distance = sqrt((Δrow)² + (Δcol)²)
    """
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)


def hex_to_rgba(h, alpha=1.0):
    """
    Convert a HTML hex colour string like "#e040fb" into a tuple of four
    floats (red, green, blue, alpha) each in the range 0.0 – 1.0.
    matplotlib needs colours in this (R,G,B,A) format for pixel-level drawing.
    """
    h = h.lstrip('#')                                      # Remove the '#'
    r, g, b = (int(h[i:i+2], 16)/255.0 for i in (0, 2, 4)) # Parse 2 hex digits each
    return (r, g, b, alpha)


def zone_centroid_float(zone):
    """
    Find the centroid (centre of mass) of a junction ZONE.
    A zone is a set of (row, col) pixel coordinates.
    Returns a (float_row, float_col) tuple — i.e. an exact average position.
    """
    ys = [p[0] for p in zone]   # All row values
    xs = [p[1] for p in zone]   # All column values
    return (sum(ys)/len(ys), sum(xs)/len(xs))


def zone_centroid(zone):
    """
    Same as zone_centroid_float but rounds to the nearest integer pixel.
    Used when we need to draw something AT the centroid.
    """
    cr, cc = zone_centroid_float(zone)
    return (int(round(cr)), int(round(cc)))


# ==============================================================================
# ███  SECTION 6 — SKELETON-TO-GRAPH CONVERSION  ███
# ==============================================================================
# The skeleton is a binary image (True = skeleton pixel, False = background).
# For pathfinding we need a GRAPH — a data structure where each skeleton pixel
# is a NODE and each pair of adjacent skeleton pixels is an EDGE with a weight
# equal to the distance between them (1.0 for straight, √2 ≈ 1.414 for diagonal).
# ==============================================================================

def skeleton_to_graph(skel):
    """
    Convert a 2-D boolean skeleton image into a weighted adjacency graph.

    Returns a defaultdict where:
        graph[(row, col)]  ->  list of  ((neighbour_row, neighbour_col), weight)

    We check all 8 neighbours (4 cardinal + 4 diagonal directions).
    Cardinal edges have weight 1.0; diagonal edges have weight √2 ≈ 1.414.
    """
    # Build a Python set of all skeleton pixel coordinates for O(1) lookup.
    # np.where(skel) returns two arrays: [all_row_indices, all_col_indices].
    # zip pairs them up into (row, col) tuples.  ipt converts to plain ints.
    graph    = defaultdict(list)
    skel_set = {ipt(p) for p in zip(*np.where(skel))}

    # All 8 neighbour directions as (Δrow, Δcol, edge_weight)
    DIRS = [
        (-1,  0, 1.000),   # Up
        ( 1,  0, 1.000),   # Down
        ( 0, -1, 1.000),   # Left
        ( 0,  1, 1.000),   # Right
        (-1, -1, 1.414),   # Up-Left   (diagonal)
        (-1,  1, 1.414),   # Up-Right  (diagonal)
        ( 1, -1, 1.414),   # Down-Left (diagonal)
        ( 1,  1, 1.414),   # Down-Right(diagonal)
    ]

    for y, x in skel_set:
        for dy, dx, w in DIRS:
            nb = (y+dy, x+dx)   # Candidate neighbour pixel
            if nb in skel_set:   # Only add the edge if that neighbour is on the skeleton
                graph[(y,x)].append((nb, w))

    return graph


# ==============================================================================
# ███  SECTION 7 — TIP DETECTION  ███
# ==============================================================================
# A TIP is a skeleton pixel with exactly ONE skeleton neighbour — it's a dead
# end, like the outermost tip of a hypha.
# ==============================================================================

def detect_tips(skel):
    """
    Return a list of (row, col) integer tuples for every TIP pixel in the
    skeleton image 'skel'.

    Method:
      1. Use a 3×3 convolution kernel (a sliding window) to count, for each
         pixel, how many of its 8 neighbours are also skeleton pixels.
      2. A pixel with exactly 1 skeleton neighbour is a tip.
         (A pixel with 0 neighbours would be isolated; ≥3 is a junction.)
    """
    # k is a 3×3 array of 1s with the centre set to 0.
    # Convolving the skeleton with k at each pixel gives the number of
    # skeleton neighbours at that pixel.
    k   = np.ones((3,3), dtype=int)
    k[1,1] = 0                                    # Don't count the pixel itself
    deg = convolve(skel.astype(int), k, mode='constant') * skel
    # deg[r,c] = number of skeleton neighbours of pixel (r,c), or 0 if not on skel.
    # A tip has deg == 1.
    return [ipt(p) for p in np.argwhere((skel==1) & (deg==1))]


# ==============================================================================
# ███  SECTION 8 — JUNCTION ZONE DETECTION  ███
# ==============================================================================
# A JUNCTION is a skeleton pixel with 3 or more skeleton neighbours — a branch
# point where paths fork or merge.
#
# Nearby junction pixels are grouped into ZONES (clusters) because a real
# physical junction might span several pixels in the image.
# ==============================================================================

def detect_junction_zones(skel, graph):
    """
    Find all junction zones in the skeleton.

    Returns:
        zones         — list of sets; each set contains the (row,col) pixels
                        belonging to that zone.
        pixel_to_zone — dict mapping each junction pixel to its zone index.

    Algorithm:
      1. Find all pixels with ≥3 skeleton neighbours (raw junction pixels).
      2. Group touching raw junction pixels into clusters using flood-fill
         (breadth-first search starting from each unvisited junction pixel).
    """
    # Count skeleton neighbours using the same convolution trick as detect_tips.
    k   = np.ones((3,3), dtype=int)
    k[1,1] = 0
    deg = convolve(skel.astype(int), k, mode='constant') * skel

    # raw_set: pixels with 3 or more skeleton neighbours
    raw_set = {ipt(p) for p in np.argwhere((skel==1) & (deg>=3))}

    visited       = set()     # Track which pixels we've already processed
    zones         = []        # Final list of zone pixel-sets
    pixel_to_zone = {}        # Map: pixel -> zone index

    for seed in raw_set:
        if seed in visited:
            continue   # Already part of a previously found cluster — skip

        # BFS flood-fill: collect all pixels connected to 'seed' via raw_set
        cluster = set()
        queue   = [seed]
        visited.add(seed)

        while queue:
            cur = queue.pop()          # Take next pixel from queue
            cluster.add(cur)           # Add it to this cluster
            cy, cx = cur
            # Check all 8 neighbours
            for dy in range(-1, 2):
                for dx in range(-1, 2):
                    nb = (cy+dy, cx+dx)
                    if nb in raw_set and nb not in visited:
                        visited.add(nb)
                        queue.append(nb)

        # Record this cluster as a new zone
        zid = len(zones)
        zones.append(cluster)
        for px in cluster:
            pixel_to_zone[px] = zid   # Every pixel in this cluster maps to zid

    return zones, pixel_to_zone


# ==============================================================================
# ███  SECTION 9 — JUNCTION ID TRACKING ACROSS FRAMES  ███
# ==============================================================================
# Junctions don't have exactly matching positions across frames, the network 
# grows and shifts.  We assign persistent integer IDs by matching each new frame's
# junction centroids to the previous frame's centroids within a search radius.
# ==============================================================================

def assign_ids(zones, prev_centroids, next_id, match_radius):
    """
    Assign persistent integer IDs to junction zones by matching them to the
    zones from the previous frame.

    Parameters:
        zones          — list of zone pixel-sets for the current frame
        prev_centroids — list of (row, col, id) tuples from the PREVIOUS frame
        next_id        — the next unused integer ID
        match_radius   — max distance (pixels) to consider a match

    Returns:
        ids            — list of integer IDs, one per zone
        new_prev       — updated list of (row, col, id) to pass to the next frame
        next_id        — updated next unused ID
    """
    if not zones:
        return [], [], next_id   # Nothing to do if there are no zones

    # Calculate centroid of each zone as floating-point (row, col)
    centroids = [zone_centroid_float(z) for z in zones]
    ids       = [None] * len(zones)   # Start with no IDs assigned

    # Build a list of all (distance, new_zone_idx, prev_zone_idx, prev_id) candidates
    # where the distance is within match_radius.
    candidates = []
    for ni, (nr, nc) in enumerate(centroids):
        for pi, (pr, pc, pid) in enumerate(prev_centroids):
            d = math.hypot(nr - pr, nc - pc)
            if d <= match_radius:
                candidates.append((d, ni, pi, pid))

    # Sort by distance — best matches first
    candidates.sort()

    # Greedy matching: assign the closest pair, then exclude both from further matching
    used_new  = set()
    used_prev = set()
    for d, ni, pi, pid in candidates:
        if ni in used_new or pi in used_prev:
            continue   # One of these has already been matched
        ids[ni] = pid
        used_new.add(ni)
        used_prev.add(pi)

    # Any zone that wasn't matched gets a brand-new ID
    for ni in range(len(zones)):
        if ids[ni] is None:
            ids[ni] = next_id
            next_id += 1

    # Build the "previous centroids" list to pass into the next frame
    new_prev = [(centroids[i][0], centroids[i][1], ids[i]) for i in range(len(zones))]
    return ids, new_prev, next_id


def run_tracker(match_radius):
    """
    Run the junction-ID tracker across ALL frames in sequence.

    For each frame, calls assign_ids() passing in the previous frame's
    centroids.  Updates each frame's 'junction_ids' and 'centroids' fields.

    Returns the total number of unique junction IDs used.
    """
    prev_centroids = []   # Empty before the first frame
    next_id        = 1    # IDs start from 1

    for fd in frame_data:
        ids, prev_centroids, next_id = assign_ids(
            fd['zones'], prev_centroids, next_id, match_radius)
        fd['junction_ids'] = ids
        fd['centroids']    = [zone_centroid_float(z) for z in fd['zones']]

    return next_id - 1   # Total unique IDs assigned (they start at 1)


def build_tip_labels_from_pairings():
    """
    Build the 'tip_ids' lookup table for every frame.

    'tip_ids[t]' is a dict mapping each tip pixel to a human-readable label:
      - Tips paired to a JUNCTION get the label  "J<global_id>"
        e.g. "J3" if the junction's persistent ID is 3.
      - Tips paired to the ORIGIN get unique labels: "orig1", "orig2", …
        (unique within a frame; the first origin-paired tip is "orig1", etc.)
      - FAUX tips get "faux_N" labels derived from their ID.

    This function is called every time pairings change (after overrides,
    re-tracking, origin capacity changes, etc.).
    """
    # Iterate over each time point
    for t in range(num_frames):
        fd          = frame_data[t]
        j_ids       = fd['junction_ids']
        t_ids       = {}
        orig_counter = 0   # How many origin pairings have we seen so far?

        for seg in fd['pairings']:
            tip = seg['tip']

            if seg['zone_id'] == ORIGIN_ZONE_ID:
                # Origin pairing, give it a unique "orig1", "orig2", … label
                orig_counter += 1
                t_ids[tip] = f'orig{orig_counter}'

            else:
                # Junction pairing, look up the global persistent junction ID
                zid = seg['zone_id']
                if zid < len(j_ids):
                    t_ids[tip] = f"J{j_ids[zid]}"
                else:
                    t_ids[tip] = f"Z{zid}"   # Fallback if IDs are out of sync

        tip_ids[t] = t_ids


# ==============================================================================
# ███  SECTION 10 — ANCHOR SNAPPING  ███
# ==============================================================================
# "Snapping" means moving a point to the nearest skeleton pixel.
# We snap the anchor so it is always exactly on the skeleton.
# ==============================================================================

def snap_to_skeleton(pos, graph, radius=ANCHOR_SNAP_RADIUS):
    """
    Find the skeleton pixel nearest to 'pos' within 'radius' pixels.

    'graph' is the skeleton graph (keys are skeleton pixel coordinates).
    Returns the nearest skeleton pixel as (row, col), or None if none found.
    """
    py, px = int(pos[0]), int(pos[1])

    # If the exact pixel is already on the skeleton, no need to search
    if (py, px) in graph:
        return (py, px)
    
    # If it is not on the skeleton, then search
    best_d, best = radius+1, None
    for dy in range(-radius, radius+1):
        for dx in range(-radius, radius+1):
            c = (py+dy, px+dx)
            if c in graph:
                d = math.hypot(dy, dx)
                if d <= radius and d < best_d:
                    best_d, best = d, c
    return best   # Returns None if nothing found within radius


# ==============================================================================
# ███  SECTION 11 — ANCHOR DETECTION  ███
# ==============================================================================
# The anchor is the stable "root" pixel of the whole network.
# We find it by taking the centroid of pixels that are present in at least
# ANCHOR_MIN_FRAME_FRACTION of all frames, then snapping to the skeleton.
# ==============================================================================

def find_anchor(skeletons):
    """
    Find the global anchor pixel by averaging the "stable" part of the skeleton
    (pixels present in ≥50% of frames) across all frames.

    Returns a (row, col) integer tuple.
    """
    print("\nFinding anchor …")
    n = len(skeletons)

    # Sum all frames — heat[r,c] = number of frames that contain pixel (r,c)
    heat = sum(s.astype(float) for s in skeletons)

    # Threshold: keep only pixels present in at least 50% of frames
    thr  = ANCHOR_MIN_FRAME_FRACTION * n
    mask = heat >= thr

    # Fallback: if no pixel passes the threshold, use the first frame's skeleton
    first_skel = next((s for s in skeletons if s.any()), None)
    if not mask.any():
        ys, xs = np.where(first_skel) if first_skel is not None else ([0],[0])
    else:
        ys, xs = np.where(mask)

    # Take the mean position of the "stable" pixels as the centroid
    centroid = (int(ys.mean()), int(xs.mean()))

    # Snap that centroid to the actual skeleton of the first frame
    if first_skel is not None:
        g       = skeleton_to_graph(first_skel)
        snapped = snap_to_skeleton(centroid, g)
        anchor  = snapped if snapped else centroid
    else:
        anchor = centroid

    print(f"  centroid={centroid}  ->  anchor={anchor}")
    return anchor


def sticky_anchor_for_frame(prev_anchor, graph, skel):
    """
    Return the anchor pixel for a specific frame, keeping it "sticky"
    i.e. prefer the same pixel as the previous frame if it's still on the skeleton.

    If the previous anchor has moved off the skeleton, snap to the nearest pixel.
    If snapping fails too, fall back to the nearest skeleton pixel of all.
    """
    # Best case: anchor is already on the skeleton for this frame
    if prev_anchor in graph:
        return prev_anchor

    # Try snapping within the search radius
    snapped = snap_to_skeleton(prev_anchor, graph, radius=ANCHOR_SNAP_RADIUS)
    if snapped is not None:
        print(f"  Anchor {prev_anchor} left skeleton -> snapped to {snapped}")
        return snapped

    # Last resort: brute-force nearest skeleton pixel
    skel_pts = list(zip(*np.where(skel)))
    if not skel_pts:
        return prev_anchor   # No skeleton at all, just keep the old anchor

    py, px = prev_anchor
    best_d, best = math.inf, None
    for (sy, sx) in skel_pts:
        d = math.hypot(sy-py, sx-px)
        if d < best_d:
            best_d, best = d, (sy, sx)

    print(f"  Anchor {prev_anchor} left skeleton (no snap) -> fallback {best}")
    return best


# ==============================================================================
# ███  SECTION 12 — DIJKSTRA'S ALGORITHM (EDGE-EXCLUDING)  ███
# ==============================================================================
# Dijkstra's algorithm finds the SHORTEST PATH from one node to all other nodes
# in a weighted graph.  We use it to find the shortest skeleton path from each
# tip to each junction.
#
# The "edge-excluding" variant skips over edges that have already been used by
# a previously committed pairing, so two tips can't share the same path segment.
# ==============================================================================

def dijkstra_with_exclusions(graph, tip, pixel_to_zone, zone_states, excluded_edges):
    """
    Run Dijkstra from 'tip', finding the shortest path to every NORMAL junction
    zone, while avoiding edges in 'excluded_edges'.

    Parameters:
        graph           — the skeleton graph
        tip             — start pixel (row, col)
        pixel_to_zone   — dict mapping pixels to their zone index
        zone_states     — dict mapping zone index to its state (normal/cross/etc.)
        excluded_edges  — set of frozensets of 2 pixels whose edge is forbidden

    Returns:
        zone_hits — list of (zone_id, entry_pixel, dist, path) sorted by distance
        dist      — dict of shortest distances from tip to every pixel
        prev      — dict used to reconstruct paths (prev[pixel] = pixel before it)
    """
    # Priority queue entries are (distance, pixel)
    # heapq is a MIN-heap: the smallest distance is always at the front
    pq   = [(0.0, tip)]
    dist = {tip: 0.0}    # Distance from tip to each pixel (filled as we explore)
    prev = {}            # prev[pixel] = which pixel we came from (for path reconstruction)
    visited_zones = {}   # zone_id -> (entry_pixel, distance, path) — first time we hit each zone

    while pq:
        d, u = heapq.heappop(pq)   # Get the closest unvisited pixel

        # If we've already found a shorter path to u, skip this stale entry
        if d > dist.get(u, math.inf):
            continue

        # Check if this pixel belongs to a junction zone
        zid = pixel_to_zone.get(u)
        if zid is not None and zid not in visited_zones:
            state_val = zone_states.get(zid, ZONE_STATE_NORMAL)
            if state_val == ZONE_STATE_NORMAL:
                # First time reaching this NORMAL zone — record the entry point and path
                path, node = [], u
                while node in prev:
                    path.append(node)
                    node = prev[node]
                path.append(tip)
                path.reverse()   # Path goes from tip -> zone entry
                visited_zones[zid] = (u, d, path)

        # Explore neighbours
        for v, w in graph.get(u, []):
            edge = frozenset((u, v))
            if edge in excluded_edges:
                continue   # This edge is reserved — skip it
            nd = d + w
            if nd < dist.get(v, math.inf):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))

    # Convert visited_zones dict into a sorted list
    zone_hits = [(zid, ep, d, p) for zid,(ep,d,p) in visited_zones.items()]
    zone_hits.sort(key=lambda x: x[2])   # Sort by distance (closest first)
    return zone_hits, dist, prev


def path_to_pixel(prev_map, tip, target):
    """
    Reconstruct the path from 'tip' to 'target' using the 'prev_map' dictionary
    built by Dijkstra.

    Returns a list of pixels from tip -> target, or None if no path exists.
    """
    if target not in prev_map and target != tip:
        return None   # Target was never reached

    path, node = [], target
    while node in prev_map:
        path.append(node)
        node = prev_map[node]
    path.append(tip)
    path.reverse()   # We built it backwards; reverse to get tip->target
    return path


def edges_of_path(path):
    """
    Return the set of edges in a path as frozensets of pixel pairs.
    Used to populate the 'excluded_edges' set after a pairing is committed.

    Example: path = [(0,0), (0,1), (0,2)]
    -> { frozenset({(0,0),(0,1)}), frozenset({(0,1),(0,2)}) }
    """
    s = set()
    for i in range(len(path)-1):
        s.add(frozenset((path[i], path[i+1])))
    return s


# ==============================================================================
# ███  SECTION 13 — DIRECT DIJKSTRA (for manual overrides)  ███
# ==============================================================================
# A simpler version of Dijkstra that just finds the shortest path from A to B
# with no zone-exclusion logic.  Used when the user manually picks a target.
# ==============================================================================

def dijkstra_direct(graph, start, end):
    """
    Find the shortest path from 'start' to 'end' in 'graph'.
    No edge exclusions, just the plain shortest path.

    Returns (path, distance) where path is a list of pixels tip->end,
    or (None, inf) if no path exists.
    """
    pq   = [(0.0, start)]
    dist = {start: 0.0}
    prev = {}

    while pq:
        d, u = heapq.heappop(pq)
        if u == end:
            break   # Found the target — stop early
        if d > dist.get(u, math.inf):
            continue
        for v, w in graph.get(u, []):
            nd = d + w
            if nd < dist.get(v, math.inf):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))

    if end not in dist:
        return None, math.inf   # No path found

    # Reconstruct path
    path, node = [], end
    while node in prev:
        path.append(node)
        node = prev[node]
    path.append(start)
    path.reverse()
    return path, dist[end]


# ==============================================================================
# ███  SECTION 14 — GREEDY PAIRING ALGORITHM  ███
# ==============================================================================
# This is the core algorithm that matches tips to junctions.
#
# The strategy is "greedy by distance":
#   1. Pre-compute the nearest junction for every tip (ignoring exclusions).
#   2. Sort tips by their best-case distance (shortest-distance tip goes first).
#   3. For each tip (in order), find the nearest AVAILABLE junction and claim it.
#      "Available" = not already claimed by a closer tip, and not a crossover/removed.
#   4. After claiming a junction, mark all edges of the path as excluded so the
#      next tip can't use the same skeleton segments.
# ==============================================================================

def compute_pairings(tips, graph, zones, pixel_to_zone, zone_states,
                     origin_pixel, origin_capacity):
    """
    Compute greedy tip->junction pairings.

    Parameters:
        tips             — list of (row,col) tip pixels
        graph            — skeleton graph
        zones            — list of zone pixel-sets
        pixel_to_zone    — dict pixel->zone_id
        zone_states      — dict zone_id->state
        origin_pixel     — the origin/anchor pixel (or None)
        origin_capacity  — how many tips are allowed to pair to the origin

    Returns:
        pairings      — list of pairing dicts (see structure below)
        unpaired_tips — list of tips that could not be paired
    """
    excluded_edges = set()   # Edges "used up" by committed pairings
    claimed_zones  = {}      # zone_id -> tip pixel (which tip claimed this zone)
    origin_count   = 0       # How many tips have paired to the origin so far
    pairings       = []
    unpaired_tips  = []

    # ── STEP 1: Pre-sort tips by their unconstrained nearest-junction distance ──
    # This ensures the tip that is closest to ANY junction goes first and gets
    # first pick.  Longer tips fill in the gaps afterwards.
    prelim = []
    for tip in tips:
        if tip not in graph:
            continue   # Tip is not on the skeleton (shouldn't happen, but be safe)

        # Run Dijkstra with NO exclusions to get the "ideal" distances
        zh, dm, pm = dijkstra_with_exclusions(
            graph, tip, pixel_to_zone, zone_states, set())

        eligible = [(zid, ep, d, p) for zid, ep, d, p in zh
                    if zone_states.get(zid, ZONE_STATE_NORMAL) == ZONE_STATE_NORMAL]

        # Also check if origin is reachable
        origin_entry = None
        if origin_pixel and origin_pixel in dm:
            od = dm[origin_pixel]
            op = path_to_pixel(pm, tip, origin_pixel)
            if op:
                origin_entry = (ORIGIN_ZONE_ID, origin_pixel, od, op)

        merged = sorted(eligible + ([origin_entry] if origin_entry else []),
                        key=lambda x: x[2])

        if not merged:
            unpaired_tips.append(tip)
            continue

        prelim.append((merged[0][2], tip))   # (best_distance, tip)

    prelim.sort(key=lambda x: x[0])   # Sort tips by nearest-junction distance

    # ── STEP 2: Greedy assignment with edge exclusion ──────────────────────────
    for _, tip in prelim:
        if tip not in graph:
            unpaired_tips.append(tip)
            continue

        # Re-run Dijkstra WITH the current excluded_edges so already-used paths
        # are avoided.
        zh, dm, pm = dijkstra_with_exclusions(
            graph, tip, pixel_to_zone, zone_states, excluded_edges)

        eligible = [(zid, ep, d, p) for zid, ep, d, p in zh
                    if zone_states.get(zid, ZONE_STATE_NORMAL) == ZONE_STATE_NORMAL]

        origin_entry = None
        if origin_pixel and origin_pixel in dm:
            od = dm[origin_pixel]
            op = path_to_pixel(pm, tip, origin_pixel)
            if op:
                origin_entry = (ORIGIN_ZONE_ID, origin_pixel, od, op)

        merged = sorted(eligible + ([origin_entry] if origin_entry else []),
                        key=lambda x: x[2])

        paired = False
        for zid, zone_entry_px, d, path in merged:
            if zid == ORIGIN_ZONE_ID:
                # Try to claim the origin
                if origin_count < origin_capacity:
                    origin_count += 1
                    excluded_edges |= edges_of_path(path)
                    pairings.append({
                        'tip'          : tip,
                        'zone_id'      : ORIGIN_ZONE_ID,
                        'zone_entry_px': zone_entry_px,
                        'dist_px'      : d,
                        'path'         : path,
                        'manual'       : False,
                    })
                    paired = True
                    break
            else:
                # Try to claim this junction zone
                if zid not in claimed_zones:
                    claimed_zones[zid] = tip
                    excluded_edges |= edges_of_path(path)
                    pairings.append({
                        'tip'          : tip,
                        'zone_id'      : zid,
                        'zone_entry_px': zone_entry_px,
                        'dist_px'      : d,
                        'path'         : path,
                        'manual'       : False,
                    })
                    paired = True
                    break

        if not paired:
            unpaired_tips.append(tip)

    return pairings, unpaired_tips


# ==============================================================================
# ███  SECTION 15 — PER-FRAME BASE DATA  ███
# ==============================================================================
# This section sets up the data dictionary for a single frame and provides the
# refresh_pairings function that recomputes pairings whenever something changes.
# ==============================================================================

def compute_base_frame_data(skel, anchor):
    """
    Build the initial data dictionary for one frame.

    This is called ONCE at startup for each frame and stores the skeleton,
    graph, tips, junction zones, and origin pixel.  Pairings are computed
    separately by refresh_pairings().

    The returned dict has these keys:
        skel            — 2-D boolean numpy array (the skeleton image)
        graph           — skeleton adjacency graph
        tips            — list of tip (row,col) pixels
        zones           — list of zone pixel-sets
        pixel_to_zone   — dict pixel->zone_id
        origin_pixel    — the snapped anchor/origin pixel
        skeleton_area   — number of skeleton pixels (used in frame_summary.csv)
        junction_ids    — [] until run_tracker() fills it in
        centroids       — [] until run_tracker() fills it in
        pairings        — [] until refresh_pairings() fills it in
        unpaired_tips   — [] until refresh_pairings() fills it in
        manual_overrides— dict of tip->pairing_dict for manually overridden tips
        faux_tips       — dict of faux_pixel->faux_id_string for ghost tips
    """
    graph      = skeleton_to_graph(skel)
    tips       = detect_tips(skel)
    zones, p2z = detect_junction_zones(skel, graph)
    origin_px  = snap_to_skeleton(anchor, graph, radius=ANCHOR_SNAP_RADIUS)

    return {
        'skel'            : skel,
        'graph'           : graph,
        'tips'            : tips,
        'zones'           : zones,
        'pixel_to_zone'   : p2z,
        'origin_pixel'    : origin_px,
        'skeleton_area'   : int(skel.sum()),
        'junction_ids'    : [],
        'centroids'       : [],
        'pairings'        : [],
        'unpaired_tips'   : [],
        'manual_overrides': {},    # tip -> pairing dict
        'faux_tips'       : {},    # faux_pixel -> faux_id string (e.g. "faux_3")
    }


def refresh_pairings(fd, zone_states, origin_capacity):
    """
    Recompute pairings for a single frame dictionary 'fd'.

    1. Identify tips with manual overrides.
    2. Run the pairing algorithm ONLY on the NON-override tips.
        - The manual override paths' edges are pre-added to 'excluded_edges'
            inside compute_pairings so other tips route around them.
        - The manually overridden junctions are pre-marked as "claimed" so
            other tips can't steal them.
    3. Merge algorithm pairings + manual override pairings.

    This guarantees each junction is claimed by at most one tip.
    """
    manual = fd.get('manual_overrides', {})   # tip -> pairing dict

    # ── Collect all faux tip pixels so they participate in pairing ─────────────
    faux_tips = fd.get('faux_tips', {})        # faux_pixel -> faux_id_string

    # The "effective" tip list is: real tips  +  faux tip pixels
    real_tips  = fd['tips']
    faux_pixels = list(faux_tips.keys())
    all_tips    = real_tips + faux_pixels

    if not manual:
        # ── Simple case: no overrides — just run the algorithm on all tips ────
        pairings, unpaired = compute_pairings(
            all_tips, fd['graph'], fd['zones'],
            fd['pixel_to_zone'], zone_states,
            fd['origin_pixel'], origin_capacity)
    else:
        # Tips that have manual overrides — don't feed into the algorithm
        override_tip_set = set(manual.keys())

        # Tips that need automatic pairing (real + faux, minus override tips)
        auto_tips = [t for t in all_tips if t not in override_tip_set]

        # Pre-seed excluded_edges with edges from every manual override path
        # so that the auto-tips route AROUND the manually chosen paths.
        pre_excluded = set()
        for seg in manual.values():
            pre_excluded |= edges_of_path(seg['path'])

        # Run compute_pairings normally, then de-duplicate:
        # after the algorithm runs, if an algorithm pairing claims the same zone
        # as a manual override, drop the algorithm pairing.

        pairings, unpaired = compute_pairings(
            auto_tips, fd['graph'], fd['zones'],
            fd['pixel_to_zone'], zone_states,
            fd['origin_pixel'], origin_capacity)

        # De-duplicate: remove any algorithm pairing that conflicts with
        # a manual override's claimed junction 
        # Build the set of (zone_id) claimed by manual overrides.
        # Origin pairings use ORIGIN_ZONE_ID (-1), each is unique by tip.
        manual_claimed_zones = set()
        for seg in manual.values():
            if seg['zone_id'] != ORIGIN_ZONE_ID:
                manual_claimed_zones.add(seg['zone_id'])
            # Origin can be claimed multiple times up to origin_capacity,
            # so we don't add ORIGIN_ZONE_ID here... the capacity check in
            # compute_pairings handles it.

        # Keep only algorithm pairings whose target zone is NOT already taken
        pairings = [p for p in pairings
                    if p['zone_id'] not in manual_claimed_zones]

        # Similarly, make sure unpaired doesn't include override tips
        unpaired = [t for t in unpaired if t not in override_tip_set]

        # Add the manual overrides to the final pairings list
        for seg in manual.values():
            pairings.append(seg)

    fd['pairings']      = pairings
    fd['unpaired_tips'] = unpaired


# ==============================================================================
# ███  SECTION 16 — JUNCTION STATE PROPAGATION  ███
# ==============================================================================
# When you mark a junction as a crossover in one frame, that same junction
# (identified by its persistent global ID) is automatically marked the same
# way in all SUBSEQUENT frames.  This saves you from having to click the same
# junction in every frame.
# ==============================================================================

def propagate_zone_state_forward(frame_idx, local_zone_idx, new_state):
    """
    Propagate a zone state change from frame 'frame_idx' forward to all later
    frames that contain a junction with the same global ID.

    Parameters:
        frame_idx       — the frame where the change was made
        local_zone_idx  — the index of the zone within that frame's zone list
        new_state       — the new ZONE_STATE_* value to assign
    """
    fd    = frame_data[frame_idx]
    j_ids = fd['junction_ids']

    if local_zone_idx >= len(j_ids):
        return   # Safety check — index out of range

    target_jid = j_ids[local_zone_idx]   # The global persistent junction ID

    # Walk forward through all subsequent frames
    for t in range(frame_idx + 1, num_frames):
        fd_t   = frame_data[t]
        jids_t = fd_t['junction_ids']

        # Find the zone in this frame that has the matching global ID
        for zi, jid in enumerate(jids_t):
            if jid == target_jid:
                zone_states_all[t][zi] = new_state
                break   # Found it — move to the next frame


# ==============================================================================
# ███  SECTION 17 — LOAD IMAGE DATA  ███
# ==============================================================================
# Read the TIFF stack from disk and set up the per-frame data structures.
# ==============================================================================

print(f"Loading: {SKELETON_TIFF}")
raw = imread(SKELETON_TIFF)   # Read the image file into a numpy array
raw = raw > 0                 # Convert to boolean (True = skeleton pixel)

# If the image is 2-D (a single frame), add a "frame" dimension so the rest
# of the code can always assume 3-D (frames × rows × cols).
if raw.ndim == 2:
    raw = raw[np.newaxis, ...]

num_frames = raw.shape[0]
print(f"✓ {num_frames} frame(s) loaded")

# Split the 3-D array into a list of 2-D per-frame arrays
skeletons = [raw[t] for t in range(num_frames)]

# Find the global anchor (stable root pixel across all frames)
anchor_pt = find_anchor(skeletons)

# Compute per-frame "sticky" anchors — each frame's anchor starts at the
# previous frame's anchor and snaps to the nearest skeleton pixel if needed.
print("Computing per-frame sticky anchors …")
frame_anchors = []
prev_anchor   = anchor_pt
for t, skel in enumerate(skeletons):
    g   = skeleton_to_graph(skel)
    anc = sticky_anchor_for_frame(prev_anchor, g, skel)
    frame_anchors.append(anc)
    prev_anchor = anc
    print(f"  Frame {t}: anchor = {anc}")

# Build the per-frame data dictionaries
print("Pre-computing per-frame base data …")
frame_data      = []                                   # One dict per frame
zone_states_all = [dict() for _ in range(num_frames)]  # zone_states[t][zid] = state
origin_capacity = 1                                    # How many tips can pair to origin
tip_ids         = [dict() for _ in range(num_frames)]  # tip_ids[t][tip] = label
faux_tip_counter = 0                                   # Global counter to give faux tips unique IDs like "faux_1"

for t, skel in enumerate(skeletons):
    fd = compute_base_frame_data(skel, frame_anchors[t])
    frame_data.append(fd)
    print(f"  Frame {t}: {len(fd['tips'])} tips, {len(fd['zones'])} zones")

# Run the junction ID tracker to assign persistent IDs across frames
total_unique_ids = run_tracker(DEFAULT_MATCH_RADIUS)
print(f"✓ Done  ({total_unique_ids} unique junction IDs)\n")


# ==============================================================================
# ███  SECTION 18 — GUI STATE VARIABLES  ███
# ==============================================================================
# Variables that track the current state of the interactive GUI.
# These are changed by mouse/keyboard events and read by the render functions.
# ==============================================================================

current_frame   = 0              # Which frame is currently displayed (0-based)
current_phase   = PHASE_REVIEW   # Phase 1 or Phase 2
hovered_path    = None           # Index of the path the mouse is hovering over
hovered_junc    = None           # Index of the junction the mouse is hovering over
zoom_level      = 1.0            # 1.0 = no zoom; >1 = zoomed in
pan_offset      = [0.0, 0.0]     # How far the view has been panned (col, row)
is_panning      = False          # Is the user currently dragging to pan?
pan_start       = [0.0, 0.0]     # Mouse position when panning started
dilation_radius = 0              # How many pixels to "inflate" skeleton pixels cosmetically
match_radius    = DEFAULT_MATCH_RADIUS   # Junction ID matching radius (adjustable via slider)
flash_counter   = 0              # Countdown for the "✓ Re-tracked" flash message


# ==============================================================================
# ███  SECTION 19 — COLOUR HELPERS  ███
# ==============================================================================

# A list of distinct colours used to colour each path/tip pair.
# They cycle — path 0 gets index 0, path 10 gets index 0 again, etc.
PAIR_COLORS = ['#4fc3f7','#81c784','#ffb74d','#e57373','#ce93d8',
               '#80cbc4','#f9a825','#ff8a65','#90caf9','#a5d6a7']

# Map zone state integers to display colours
ZONE_STATE_COLORS = {
    ZONE_STATE_NORMAL  : '#e040fb',             # Purple
    ZONE_STATE_CROSS   : PALETTE['CROSSOVER'],    # Orange
    ZONE_STATE_PARALLEL: PALETTE['CROSSOVER'],    # Orange
    ZONE_STATE_REMOVED : '#555555',             # Dark grey
}

# Map zone state integers to matplotlib marker shapes
ZONE_STATE_MARKERS = {
    ZONE_STATE_NORMAL  : 'o',   # Circle
    ZONE_STATE_CROSS   : 'D',   # Diamond
    ZONE_STATE_PARALLEL: 's',   # Square
    ZONE_STATE_REMOVED : 'x',   # X
}


def pair_color(i):
    """Return the colour for the i-th pairing (cycles through PAIR_COLORS)."""
    return PAIR_COLORS[i % len(PAIR_COLORS)]


# ==============================================================================
# ███  SECTION 20 — DILATION HELPERS  ███
# ==============================================================================
# "Dilation" makes skeleton pixels appear thicker on screen.
# It has NO effect on measurements, it's purely cosmetic for easier viewing and exporting figures.
# ==============================================================================

def dilated_skel(skel, radius):
    """
    Return a version of 'skel' where every True pixel has been "inflated"
    by 'radius' pixels in all directions.  radius=0 means no dilation.
    """
    if radius <= 0:
        return skel
    struct = np.ones((2*radius+1, 2*radius+1), dtype=bool)   # Square structuring element
    return binary_dilation(skel, structure=struct)


def dilation_scale():
    """
    Return a scale factor for marker sizes based on the dilation radius.
    When the skeleton is drawn thicker, markers should be bigger to match.
    """
    return 1.0 + dilation_radius * 0.35


# ==============================================================================
# ███  SECTION 21 — BLOB RENDERER  ███
# ==============================================================================
# Renders a set of pixels as a coloured blob on a matplotlib axis.
# Used to draw junctions, tips, and origin markers as visible blobs.
# ==============================================================================

def render_blob(ax, skel_shape, pixels, hex_color, alpha, dilation_rad, zorder):
    """
    Draw a coloured blob on 'ax' at the given 'pixels' positions.

    Parameters:
        ax          — the matplotlib Axes object to draw on
        skel_shape  — (rows, cols) shape of the skeleton image
        pixels      — list/set of (row, col) coordinates to paint
        hex_color   — HTML colour string like "#ff1744"
        alpha       — opacity (0.0 = invisible, 1.0 = opaque)
        dilation_rad— inflate the blob by this many pixels
        zorder      — drawing order (higher = drawn on top)
    """
    # Create an empty image the same size as the skeleton
    mask = np.zeros(skel_shape, dtype=float)
    for (py, px_) in pixels:
        mask[py, px_] = 1.0   # Mark these pixels

    # Optionally dilate (inflate) the mask
    if dilation_rad > 0:
        struct = np.ones((2*dilation_rad+1, 2*dilation_rad+1), bool)
        mask   = binary_dilation(mask.astype(bool), structure=struct).astype(float)

    # Build an RGBA image where the mask pixels get the chosen colour
    r, g, b, _ = hex_to_rgba(hex_color)
    rgba = np.zeros((*skel_shape, 4), dtype=float)
    rgba[..., 0] = mask * r      # Red channel
    rgba[..., 1] = mask * g      # Green channel
    rgba[..., 2] = mask * b      # Blue channel
    rgba[..., 3] = mask * alpha  # Alpha (transparency) channel

    return ax.imshow(rgba, zorder=zorder)


# ==============================================================================
# ███  SECTION 22 — VIEW LIMITS  ███
# ==============================================================================
# Applies the current zoom level and pan offset to the matplotlib axis so the
# user sees the right part of the image.
# ==============================================================================

def apply_view(ax, skel_shape):
    """
    Set the x/y axis limits based on the current zoom_level and pan_offset.

    The image is centred at (w/2 + pan_offset[0],  h/2 + pan_offset[1]).
    At zoom=1 the full image is visible.  At zoom=2 only the central quarter
    of the image is visible.
    """
    h, w = skel_shape
    cx   = w/2 + pan_offset[0]     # Centre column (with pan)
    cy   = h/2 + pan_offset[1]     # Centre row    (with pan)
    hw   = (w/2) / zoom_level      # Half-width of the visible region
    hh   = (h/2) / zoom_level      # Half-height of the visible region
    ax.set_xlim(cx - hw, cx + hw)
    ax.set_ylim(cy + hh, cy - hh)  # y is flipped (image origin is top-left)


# ==============================================================================
# ███  SECTION 23 — IMBALANCE INDICATOR  ███
# ==============================================================================
# If the number of tips doesn't equal the number of active (normal) junctions,
# something is probably wrong and we warn the user.
# ==============================================================================

def imbalance_string(fd, zs):
    """
    Return a warning string if the tip count ≠ active junction count,
    or an empty string if they match.

    The idea is that ideally every tip pairs to exactly one junction and
    vice versa.  A mismatch means some tips will be unpaired or some
    junctions will be unclaimed.
    """
    n_tips   = len(fd['tips'])
    n_active = sum(1 for zid in range(len(fd['zones']))
                   if zs.get(zid, ZONE_STATE_NORMAL) == ZONE_STATE_NORMAL)
    if n_tips > n_active:
        return f"  ⚠ +{n_tips - n_active} TIPS"
    elif n_active > n_tips:
        return f"  ⚠ +{n_active - n_tips} JUNCTIONS"
    return ""


# ==============================================================================
# ███  SECTION 24 — ORIGIN CLICK DETECTION  ███
# ==============================================================================

def find_clicked_origin(xdata, ydata):
    """
    Return True if the mouse click at (xdata, ydata) is within
    CLICK_ZONE_RADIUS_PX pixels of the current frame's origin pixel.

    Used in override mode to detect "user clicked the origin".
    Note: matplotlib gives us (x=col, y=row) in data coordinates.
    """
    origin = frame_data[current_frame]['origin_pixel']
    if origin is None or xdata is None or ydata is None:
        return False
    oy, ox = origin   # origin is stored as (row, col)
    return math.hypot(xdata - ox, ydata - oy) <= CLICK_ZONE_RADIUS_PX


# ==============================================================================
# ███  SECTION 25 — RENDER: PHASE 1 (Junction Review)  ███
# ==============================================================================
# Draws the junction review view — shows the skeleton with coloured junction
# blobs and centroid markers.  No paths are shown here.
# ==============================================================================

# Pre-compute a circle for drawing the match-radius ring around each junction.
# We compute it once here to avoid recalculating 360 points every render frame.
_theta = np.linspace(0, 2 * math.pi, 120)
_ux    = np.cos(_theta)   # x-components of a unit circle
_uy    = np.sin(_theta)   # y-components of a unit circle


def render_review():
    """
    Draw Phase 1: the junction review view.

    Shows:
      - The (optionally dilated and inverted) skeleton
      - Each junction zone as a coloured blob
      - A centroid dot with the junction's persistent ID label
      - (Optionally) the match-radius circle around each junction
      - A hover tooltip listing distances to nearby junctions
    """
    show_circles = show_circles_var.get()   # Is the "show circles" checkbox ticked?
    ax.clear()                              # Wipe the previous frame's drawing

    fd    = frame_data[current_frame]
    skel  = fd['skel']
    zones = fd['zones']
    zs    = zone_states_all[current_frame]
    cents = fd['centroids']
    j_ids = fd['junction_ids']
    ds    = dilation_scale()   # Marker size scale factor

    # ── Background colour ──────────────────────────────────────────────────────
    if invert_view:
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

    # ── Skeleton ───────────────────────────────────────────────────────────────
    d_skel = dilated_skel(skel, dilation_radius)   # Cosmetically fattened skeleton
    if not invert_view:
        ax.imshow(d_skel, cmap='gray',   alpha=0.20, zorder=0)   # Dark background
    else:
        ax.imshow(d_skel, cmap='gray_r', alpha=1.00, zorder=0)   # White background

    # ── Junction blobs ─────────────────────────────────────────────────────────
    for zid, zone in enumerate(zones):
        state_val = zs.get(zid, ZONE_STATE_NORMAL)
        if state_val == ZONE_STATE_NORMAL:
            col = PALETTE['JUNCTION_UNUSED']; alpha = 0.85
        elif state_val in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL):
            col = PALETTE['CROSSOVER'];       alpha = 0.90
        else:
            col = '#555555';                  alpha = 0.60
        render_blob(ax, skel.shape, zone, col, alpha, dilation_radius, zorder=2)

    # ── Centroid markers and labels ────────────────────────────────────────────
    for zid, zone in enumerate(zones):
        state_val = zs.get(zid, ZONE_STATE_NORMAL)
        cy, cx    = zone_centroid(zone)

        if state_val == ZONE_STATE_NORMAL:
            col  = ZONE_STATE_COLORS[ZONE_STATE_NORMAL]
            mkr  = 'o'; sz = 28*ds; alph = 0.55
        else:
            col  = ZONE_STATE_COLORS[state_val]
            mkr  = ZONE_STATE_MARKERS[state_val]
            sz   = 55*ds; alph = 0.95

        ax.scatter(cx, cy, c=col, s=sz, marker=mkr,
                   edgecolors='white', linewidths=0.8, zorder=13, alpha=alph)

    # ── ID labels + hover tooltip ──────────────────────────────────────────────
    for i, (cr, cc) in enumerate(cents):
        jid    = j_ids[i] if i < len(j_ids) else '?'
        is_hov = (hovered_junc == i)   # Is the mouse hovering this junction?

        # Optionally draw the match-radius ring
        if show_circles:
            cx_pts = cc + match_radius * _ux
            cy_pts = cr + match_radius * _uy
            ax.plot(cx_pts, cy_pts,
                    color='#ffffff' if is_hov else '#ff0000',
                    linewidth=2.0   if is_hov else 1.8,
                    linestyle='--',
                    alpha=0.95      if is_hov else 0.45, zorder=3)

        # Draw the centroid dot
        ax.plot(cc, cr, 'o',
                markersize=15*ds if is_hov else 11*ds,
                color='#ffffff'  if is_hov else '#1a1a2e',
                alpha=0.95 if is_hov else 0.75, zorder=6)

        # Draw the ID number label
        ax.text(cc, cr, str(jid),
                fontsize=6, color='#1a1a2e' if is_hov else 'white',
                fontweight='bold', ha='center', va='center', zorder=7)

        # If hovering, show a tooltip with distances to nearby junctions
        if is_hov and len(cents) > 1:
            dists = []
            for j, (or_, oc) in enumerate(cents):
                if j == i: continue
                d_px = math.hypot(cr-or_, cc-oc)
                dists.append((d_px, d_px*PIXEL_SIZE_UM,
                               j_ids[j] if j < len(j_ids) else '?', or_, oc))
            dists.sort()
            inside  = [(d,u,jid2,r,c) for d,u,jid2,r,c in dists if d <= match_radius]
            outside = [(d,u,jid2,r,c) for d,u,jid2,r,c in dists if d >  match_radius]

            # Draw yellow lines to junctions within matching radius
            for d_px, _, _, or_, oc in inside:
                ax.plot([cc, oc], [cr, or_],
                        color='#ffd600', linewidth=1.0,
                        alpha=0.70, linestyle='-', zorder=7)

            # Build the tooltip text box
            lines = [f"J{jid}  (r={match_radius}px)"]
            if inside:
                lines.append("  within radius:")
                for d_px, d_um, jid2, _, _ in inside[:8]:
                    lines.append(f"   ● J{jid2:<4}  {d_px:6.1f}px  {d_um:5.1f}µm")
            if outside:
                lines.append("  outside:")
                for d_px, d_um, jid2, _, _ in outside[:4]:
                    lines.append(f"   ○ J{jid2:<4}  {d_px:6.1f}px  {d_um:5.1f}µm")
            ax.text(cc + match_radius + 6, cr, "\n".join(lines),
                    fontsize=6.5, color='#ffd600', fontfamily='monospace',
                    va='center', ha='left',
                    bbox=dict(boxstyle='round,pad=0.35', facecolor='#0d0d1e',
                              edgecolor='#ffd600', linewidth=0.8, alpha=0.92),
                    zorder=15)

    # ── Title bar ──────────────────────────────────────────────────────────────
    zs_vals    = list(zs.values())
    n_cross    = sum(1 for s in zs_vals if s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL))
    n_rem      = sum(1 for s in zs_vals if s == ZONE_STATE_REMOVED)
    n_norm     = len(zones) - n_cross - n_rem
    all_marked = sum(1 for zs_ in zone_states_all
                     for s in zs_.values()
                     if s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL, ZONE_STATE_REMOVED))
    imb = imbalance_string(fd, zs)

    ax.set_title(
        f"PHASE 1 — Junction Review  │  "
        f"Frame {current_frame}/{num_frames-1}  │  "
        f"{len(zones)} junctions  │  "
        f"normal {n_norm}  cross/para {n_cross}  removed {n_rem}  │  "
        f"r={match_radius}px  │  "
       # f"Total marked across all frames: {all_marked}" # Commented out because it seems to not count properly
        + imb,
        fontsize=9, fontweight='bold', color='#1a3a6e')

    ax.axis('off')
    apply_view(ax, skel.shape)
    canvas.draw()
    canvas.flush_events()
    update_table_review()


def update_table_review():
    """Update the right-side text panel for Phase 1."""
    table.delete('1.0', tk.END)   # Clear existing text
    fd    = frame_data[current_frame]
    zones = fd['zones']
    zs    = zone_states_all[current_frame]
    j_ids = fd['junction_ids']
    cents = fd['centroids']

    n_cross = sum(1 for s in zs.values() if s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL))
    n_rem   = sum(1 for s in zs.values() if s == ZONE_STATE_REMOVED)
    n_norm  = len(zones) - n_cross - n_rem
    imb = imbalance_string(fd, zs)

    table.insert(tk.END,
        f"PHASE 1 — Junction Review\n"
        f"{'='*38}\n"
        f"Frame {current_frame} / {num_frames-1}\n"
        f"  Tips      : {len(fd['tips'])}\n"
        f"  Junctions : {len(zones)}{imb}\n"
        f"    normal  : {n_norm}\n"
        f"    crossover: {n_cross}\n"
        f"    removed : {n_rem}\n"
        f"\n"
        f"L-click a junction to cycle:\n"
        f"  normal -> cross× -> cross∥ -> removed\n"
        f"\n"
        f"{'─'*38}\n"
        f"Junction states this frame:\n")

    state_str_map = {0:'normal',1:'cross×',2:'cross∥',3:'removed'}
    for i, jid in enumerate(j_ids):
        cr, cc    = cents[i]
        state_val = zs.get(i, ZONE_STATE_NORMAL)
        sstr      = state_str_map.get(state_val, '?')
        table.insert(tk.END,
            f"  J{jid:>3}  ({int(cr):>4},{int(cc):>4})  [{sstr}]\n")

    table.insert(tk.END, f"\n{'='*38}\nAll-frames summary:\n")
    for t in range(num_frames):
        fd_t = frame_data[t]
        zs_t = zone_states_all[t]
        nc   = sum(1 for s in zs_t.values() if s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL))
        nr   = sum(1 for s in zs_t.values() if s == ZONE_STATE_REMOVED)
        nn   = len(fd_t['zones']) - nc - nr
        mark = "◀" if t == current_frame else " "
        table.insert(tk.END,
            f"  {mark} F{t:>2}: {nn}nrm {nc}cross {nr}rem\n")


# ==============================================================================
# ███  SECTION 26 — RENDER: PHASE 2 (Analysis)  ███
# ==============================================================================
# Draws the full analysis view — skeleton, junctions, tips, origin, AND the
# coloured pairing paths with labels.
# ==============================================================================

def render_analyze():
    """
    Draw Phase 2: the full analysis / pairing view.

    Shows everything from Phase 1, plus:
      - Coloured lines for each tip->junction pairing path
      - Star markers at each tip, coloured to match their path
      - Label badges showing the tip ID (e.g. "J3", "orig1", "faux_2")
      - Origin marker with pairing count
      - Override mode banner if active
      - Faux tip mode banner if active
    """
    show_overlay = show_overlay_var.get()   # Show paths / markers?
    show_circles = show_circles_var.get()   # Show match-radius rings?

    ax.clear()
    fd       = frame_data[current_frame]
    skel     = fd['skel']
    pairings = fd['pairings']
    zones    = fd['zones']
    zs       = zone_states_all[current_frame]
    origin   = fd['origin_pixel']
    t_ids    = tip_ids[current_frame]
    j_ids    = fd['junction_ids']
    cents    = fd['centroids']
    faux_tips = fd.get('faux_tips', {})   # faux_pixel -> faux_id string

    ds     = dilation_scale()
    d_skel = dilated_skel(skel, dilation_radius)

    paired_zone_ids      = {p['zone_id'] for p in pairings}
    origin_pairing_count = sum(1 for p in pairings if p['zone_id'] == ORIGIN_ZONE_ID)

    # ── Background ─────────────────────────────────────────────────────────────
    if invert_view:
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

    if not invert_view:
        ax.imshow(d_skel, cmap='gray',   alpha=0.20, zorder=0)
    else:
        ax.imshow(d_skel, cmap='gray_r', alpha=1.00, zorder=0)

    # ── Junction blobs ─────────────────────────────────────────────────────────
    for zid, zone in enumerate(zones):
        state_val = zs.get(zid, ZONE_STATE_NORMAL)
        if state_val == ZONE_STATE_NORMAL:
            col   = PALETTE['JUNCTION_USED'] if zid in paired_zone_ids else PALETTE['JUNCTION_UNUSED']
            alpha = 0.85
        elif state_val in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL):
            col = PALETTE['CROSSOVER']; alpha = 0.90
        else:
            col = '#555555';            alpha = 0.60
        render_blob(ax, skel.shape, zone, col, alpha, dilation_radius, zorder=2)

    # ── Tips — colour real tips and faux tips differently ─────────────────────
    # In override mode the selected tip glows green; others are normal
    if override_mode == OVERRIDE_MODE_TIP and override_tip is not None:
        other_real  = [t for t in fd['tips']           if t != override_tip]
        other_faux  = [t for t in faux_tips.keys()     if t != override_tip]
        if other_real:
            render_blob(ax, skel.shape, other_real, PALETTE['TIP'],    0.9, dilation_radius, zorder=3)
        if other_faux:
            render_blob(ax, skel.shape, other_faux, '#ffd600',          0.9, dilation_radius, zorder=3)
        render_blob(ax, skel.shape, [override_tip], '#00ff88', 1.0, max(dilation_radius,3), zorder=4)
    elif faux_mode == FAUX_MODE_WAITING:
        # In faux-tip mode, draw all real tips in teal to signal "click one of these"
        render_blob(ax, skel.shape, fd['tips'],        '#00e5ff',   0.9, dilation_radius, zorder=3)
        render_blob(ax, skel.shape, list(faux_tips.keys()), '#ffd600', 0.9, dilation_radius, zorder=3)
    else:
        # Normal mode: real tips in blue, faux tips in gold
        render_blob(ax, skel.shape, fd['tips'],              PALETTE['TIP'], 0.9, dilation_radius, zorder=3)
        if faux_tips:
            render_blob(ax, skel.shape, list(faux_tips.keys()), '#ffd600',    0.9, dilation_radius, zorder=3)

    # ── Origin marker ──────────────────────────────────────────────────────────
    if origin:
        render_blob(ax, skel.shape, [origin], PALETTE['ORIGIN'], 0.95, dilation_radius, zorder=4)

    # ── Junction centroid dots + match-radius circles ─────────────────────────
    for i, (cr, cc) in enumerate(cents):
        jid    = j_ids[i] if i < len(j_ids) else '?'
        is_hov = (hovered_junc == i)

        if show_circles:
            cx_pts = cc + match_radius * _ux
            cy_pts = cr + match_radius * _uy
            ax.plot(cx_pts, cy_pts,
                    color='#ffffff' if is_hov else '#ff0000',
                    linewidth=2.0   if is_hov else 1.8,
                    linestyle='--',
                    alpha=0.95      if is_hov else 0.45, zorder=3)

        # In override tip-selected mode, colour junction dots green = valid target
        if override_mode == OVERRIDE_MODE_TIP:
            dot_color = '#00cc66'
        else:
            dot_color = '#ffffff' if is_hov else '#1a1a2e'

        ax.plot(cc, cr, 'o',
                markersize=15*ds if is_hov else 11*ds,
                color=dot_color,
                alpha=0.95 if is_hov else 0.75, zorder=6)
        ax.text(cc, cr, str(jid),
                fontsize=6, color='#1a1a2e' if is_hov else 'white',
                fontweight='bold', ha='center', va='center', zorder=7)

        # Hover tooltip (same as Phase 1)
        if is_hov and len(cents) > 1:
            dists = []

            for j, (or_, oc) in enumerate(cents):
                if j == i: continue
                d_px = math.hypot(cr-or_, cc-oc)
                dists.append((d_px, d_px*PIXEL_SIZE_UM,
                               j_ids[j] if j < len(j_ids) else '?', or_, oc))
                
            dists.sort()
            inside  = [(d,u,jid2,r,c) for d,u,jid2,r,c in dists if d <= match_radius]
            outside = [(d,u,jid2,r,c) for d,u,jid2,r,c in dists if d >  match_radius]

            for d_px, _, _, or_, oc in inside:
                ax.plot([cc, oc], [cr, or_],
                        color='#ffd600', linewidth=1.0,
                        alpha=0.70, linestyle='-', zorder=7)
                
            lines = [f"J{jid}  (r={match_radius}px)"]

            if inside:
                lines.append("  within radius:")
                for d_px, d_um, jid2, _, _ in inside[:8]:
                    lines.append(f"   ● J{jid2:<4}  {d_px:6.1f}px  {d_um:5.1f}µm")

            if outside:
                lines.append("  outside:")
                for d_px, d_um, jid2, _, _ in outside[:4]:
                    lines.append(f"   ○ J{jid2:<4}  {d_px:6.1f}px  {d_um:5.1f}µm")
            ax.text(cc + match_radius + 6, cr, "\n".join(lines),
                    fontsize=6.5, color='#ffd600', fontfamily='monospace',
                    va='center', ha='left',
                    bbox=dict(boxstyle='round,pad=0.35', facecolor='#0d0d1e',
                              edgecolor='#ffd600', linewidth=0.8, alpha=0.92),
                    zorder=15)

    # ── OVERLAY (paths, stars, labels, origin marker) ─────────────────────────
    if show_overlay:

        # Junction state markers (circle / diamond / square / X)
        for zid, zone in enumerate(zones):
            state_val = zs.get(zid, ZONE_STATE_NORMAL)
            cy, cx    = zone_centroid(zone)
            if state_val == ZONE_STATE_NORMAL:
                col = ZONE_STATE_COLORS[ZONE_STATE_NORMAL]
                mkr = 'o'; sz = 28*ds; alph = 0.55
            else:
                col = ZONE_STATE_COLORS[state_val]
                mkr = ZONE_STATE_MARKERS[state_val]
                sz  = 55*ds; alph = 0.95
            ax.scatter(cx, cy, c=col, s=sz, marker=mkr,
                       edgecolors='white', linewidths=0.8, zorder=13, alpha=alph)

        # Origin icon
        if origin:
            oy, ox = origin
            if override_mode == OVERRIDE_MODE_TIP:
                # In override mode, origin glows green as a valid target
                ax.scatter(ox, oy, c='#00ff88', s=220*ds, marker='P',
                           edgecolors='white', linewidths=2.0*ds, zorder=22, alpha=0.85)
                ax.text(ox+4, oy-4, "← click to route here",
                        fontsize=7, color='#00ff88', zorder=23,
                        bbox=dict(boxstyle='round,pad=0.15', fc='#0d2a1a', alpha=0.8))
            else:
                # Normal: show pairing count e.g. "1/2" if capacity is 2
                col = '#ff6e40' if origin_pairing_count > 0 else '#ff1744'
                ax.scatter(ox, oy, c=col, s=120*ds, marker='P',
                           edgecolors='white', linewidths=1.4*ds, zorder=22)
                ax.text(ox+4, oy-4, f"{origin_pairing_count}/{origin_capacity}",
                        fontsize=7, color='white', zorder=23,
                        bbox=dict(boxstyle='round,pad=0.15', fc='#b71c1c', alpha=0.7))

        # Pairing paths — draw each as a coloured line along the skeleton
        for i, seg in enumerate(pairings):
            arr = np.array(seg['path'])
            hov = (hovered_path == i)
            col = pair_color(i)
            lw  = (3.0 if hov else 1.4) * ds
            al  = 1.0 if hov else 0.60
            z   = 8   if hov else 5
            ls  = '--' if seg.get('manual') else '-'   # Dashed = manual override
            ax.plot(arr[:,1], arr[:,0], ls, color=col,
                    linewidth=lw, alpha=al, zorder=z, solid_capstyle='round')

        # Tip star markers and ID labels
        paired_tips = {seg['tip']: i for i, seg in enumerate(pairings)}
        all_displayed_tips = list(fd['tips']) + list(faux_tips.keys())

        for tip in all_displayed_tips:
            i       = paired_tips.get(tip)
            is_faux = tip in faux_tips

            # Pick marker colour: paired -> matching path colour, unpaired -> blue/gold
            if i is not None:
                col = pair_color(i)
            elif is_faux:
                col = '#ffd600'   # Gold for unpaired faux tips
            else:
                col = PALETTE['TIP']

            sz = (90 if (hovered_path == i) else 45) * ds
            # Faux tips get a different marker shape (triangle) so they look distinct
            marker = '^' if is_faux else '*'
            ax.scatter(tip[1], tip[0], c=col, s=sz, marker=marker,
                       edgecolors='white', linewidths=0.9*ds, zorder=10)

            # Build the label text
            tid = t_ids.get(tip)
            if tid:
                label_text = f"[F] {tid}" if is_faux else tid
                ax.text(tip[1]+4, tip[0]-4, label_text,
                        fontsize=8, color='white', fontweight='bold', zorder=16,
                        bbox=dict(boxstyle='round,pad=0.2', fc=col, alpha=0.80,
                                  edgecolor='white', linewidth=0.6))

    # ── Mode banners ───────────────────────────────────────────────────────────
    h_shape, w_shape = skel.shape
    cx_b = w_shape/2 + pan_offset[0]
    cy_b = h_shape/2 + pan_offset[1]
    hh_b = (h_shape/2) / zoom_level

    # Override mode banner
    if override_mode != OVERRIDE_MODE_NONE:
        if override_mode == OVERRIDE_MODE_TIP:
            msg = f"OVERRIDE: tip {override_tip} selected — click a junction OR the origin ✛"
        else:
            msg = "OVERRIDE: click a tip to select it"
        ax.text(cx_b, cy_b - hh_b + hh_b*0.08, msg,
                fontsize=11, fontweight='bold',
                color='#00ff88', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#0d2a1a',
                          edgecolor='#00ff88', linewidth=1.5, alpha=0.95),
                zorder=30)

    # Faux tip mode banner
    if faux_mode == FAUX_MODE_WAITING:
        ax.text(cx_b, cy_b - hh_b + hh_b*0.16,
                "FAUX TIP MODE: click a real tip to create a ghost copy",
                fontsize=11, fontweight='bold',
                color='#ffd600', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#2a1a00',
                          edgecolor='#ffd600', linewidth=1.5, alpha=0.95),
                zorder=30)

    # "Re-tracked" success flash banner
    if flash_counter > 0:
        ax.text(cx_b, cy_b - hh_b + hh_b*0.08,
                f"✓  Re-tracked  r={match_radius}px  │  {total_unique_ids} unique IDs",
                fontsize=11, fontweight='bold',
                color='#00ff88', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#0d2a1a',
                          edgecolor='#00ff88', linewidth=1.5, alpha=0.93),
                zorder=30)

    imb = imbalance_string(fd, zs)

    # ── Title bar ──────────────────────────────────────────────────────────────
    if show_overlay and hovered_path is not None and hovered_path < len(pairings):
        seg     = pairings[hovered_path]
        dist_um = seg['dist_px'] * PIXEL_SIZE_UM
        if seg['zone_id'] == ORIGIN_ZONE_ID:
            zlab = f"origin ({t_ids.get(seg['tip'], 'orig?')})"
        else:
            zlab = f"Zone {seg['zone_id']}"
        tid     = t_ids.get(seg['tip'], '')
        tlabel  = f"Tip {seg['tip']}" + (f"  [{tid}]" if tid else "")
        man_str = "  [MANUAL]" if seg.get('manual') else ""
        ax.set_title(
            f"PHASE 2  │  Frame {current_frame}/{num_frames-1}  │  "
            f"{tlabel} -> {zlab}  │  "
            f"{seg['dist_px']:.1f} px  /  {dist_um:.1f} µm"
            + man_str + imb,
            fontsize=10, fontweight='bold')
    else:
        n_cross  = sum(1 for s in zone_states_all[current_frame].values()
                       if s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL))
        n_faux   = len(fd.get('faux_tips', {}))
        ovr_str  = "ON" if show_overlay else "OFF"
        circ_str = "ON" if show_circles else "OFF"
        ax.set_title(
            f"PHASE 2 — Analysis  │  Frame {current_frame}/{num_frames-1}  │  "
            f"{len(zones)} zones  │  {len(pairings)} pairings  │  "
            f"{len(fd['unpaired_tips'])} unpaired  │  "
            f"{n_faux} faux  │  "
            f"orig cap {origin_capacity}  │  r={match_radius}px  │  "
            f"overlay {ovr_str}  │  circles {circ_str}  │  dil {dilation_radius}px"
            + imb,
            fontsize=9, fontweight='bold')

    ax.axis('off')
    apply_view(ax, skel.shape)
    canvas.draw()
    canvas.flush_events()
    update_table_analyze()


def update_table_analyze():
    """Update the right-side text panel for Phase 2."""
    table.delete('1.0', tk.END)
    fd        = frame_data[current_frame]
    pairings  = fd['pairings']
    zs        = zone_states_all[current_frame]
    t_ids     = tip_ids[current_frame]
    j_ids     = fd['junction_ids']
    cents     = fd['centroids']
    faux_tips = fd.get('faux_tips', {})

    n_norm  = sum(1 for s in zs.values() if s == ZONE_STATE_NORMAL)
    n_cross = sum(1 for s in zs.values() if s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL))
    n_rem   = sum(1 for s in zs.values() if s == ZONE_STATE_REMOVED)
    n_norm_total = len(fd['zones']) - len(zs) + n_norm

    show_overlay = show_overlay_var.get()
    show_circles = show_circles_var.get()
    imb = imbalance_string(fd, zs)

    # Build status string for active special modes
    override_str = ""
    if override_mode == OVERRIDE_MODE_TIP:
        override_str = f"\n⚡ OVERRIDE ACTIVE\n  Tip: {override_tip}\n  -> click a junction OR origin\n"
    elif override_mode != OVERRIDE_MODE_NONE:
        override_str = "\n⚡ OVERRIDE ACTIVE\n  -> click a tip\n"

    faux_str = ""
    if faux_mode == FAUX_MODE_WAITING:
        faux_str = "\n👻 FAUX TIP MODE\n  -> click a real tip\n"

    table.insert(tk.END,
        f"PHASE 2 — Analysis\n"
        f"{'='*38}\n"
        f"Frame {current_frame}\n"
        f"Tips: {len(fd['tips'])}  Faux: {len(faux_tips)}  Paired: {len(pairings)}\n"
        f"Junctions: {n_norm_total}  Crossovers: {n_cross}  Removed: {n_rem}\n"
        + (f"⚠ {imb.strip()}\n" if imb else "")
        + f"Skel area: {fd['skeleton_area']} px²\n"
        f"Anchor: {frame_anchors[current_frame]}\n"
        f"Origin cap: {origin_capacity}\n"
        f"Match radius: {match_radius}px\n"
        f"Overlay: {'ON' if show_overlay else 'OFF'}  "
        f"Circles: {'ON' if show_circles else 'OFF'}\n"
        + override_str
        + faux_str
        + f"{'='*38}\n\nJunction IDs:\n")

    for i, jid in enumerate(j_ids):
        cr, cc    = cents[i]
        state_val = zs.get(i, ZONE_STATE_NORMAL)
        sstr      = {0:'normal',1:'cross×',2:'cross∥',3:'removed'}.get(state_val,'?')
        table.insert(tk.END, f"  J{jid:>3}  ({int(cr):>4},{int(cc):>4})  [{sstr}]\n")

    table.insert(tk.END, f"\n{'='*38}\n\nPairings:\n")
    for i, seg in enumerate(pairings):
        marker  = "► " if i == hovered_path else "  "
        dist_um = seg['dist_px'] * PIXEL_SIZE_UM
        if seg['zone_id'] == ORIGIN_ZONE_ID:
            zlab = f"origin ({t_ids.get(seg['tip'], '?')})"
        else:
            zlab = f"Zone {seg['zone_id']}"
        tid     = t_ids.get(seg['tip'], '—')
        man_str = " [M]"  if seg.get('manual') else ""
        fau_str = " [F]"  if seg.get('is_faux') else ""
        table.insert(tk.END,
            f"{marker}[{tid}]{man_str}{fau_str} Tip {seg['tip']}\n"
            f"  -> {zlab}\n"
            f"  {seg['dist_px']:.1f} px / {dist_um:.1f} µm\n"
            f"{'─'*38}\n")

    if fd['unpaired_tips']:
        table.insert(tk.END, "\nUnpaired tips:\n")
        for tip in fd['unpaired_tips']:
            tid  = t_ids.get(tip, '—')
            fstr = " [F]" if tip in faux_tips else ""
            table.insert(tk.END, f"  [{tid}]{fstr} {tip}\n")

    n_manual = len(fd.get('manual_overrides', {}))
    if n_manual:
        table.insert(tk.END, f"\n[M] = {n_manual} manual override(s)\n")
    if faux_tips:
        table.insert(tk.END, f"\nFaux tips this frame:\n")
        for fpx, fid in faux_tips.items():
            table.insert(tk.END, f"  {fid} @ {fpx}\n")


# ==============================================================================
# ███  SECTION 27 — UNIFIED RENDER DISPATCHER  ███
# ==============================================================================

def render():
    """
    Call the correct render function based on which phase we're in.
    This is the ONE function to call whenever you want to redraw the canvas.
    """
    if current_phase == PHASE_REVIEW:
        render_review()
    else:
        render_analyze()


# ==============================================================================
# ███  SECTION 28 — PHASE TRANSITIONS  ███
# ==============================================================================

def enter_phase_review():
    """
    Switch to Phase 1 (Junction Review).
    Resets hover state and override/faux modes.
    Junction states and manual overrides are PRESERVED so you don't lose work.
    """
    global current_phase, hovered_path, hovered_junc, override_mode, override_tip, faux_mode
    current_phase = PHASE_REVIEW
    hovered_path  = hovered_junc = None
    override_mode = OVERRIDE_MODE_NONE
    override_tip  = None
    faux_mode     = FAUX_MODE_NONE
    _sync_panel_visibility()
    render()


def enter_phase_analyze():
    """
    Switch to Phase 2 (Analysis).
    Runs the junction tracker and computes all pairings.
    """
    global current_phase, total_unique_ids, hovered_path, hovered_junc
    hovered_path = hovered_junc = None

    # Re-run the ID tracker with the current match radius
    total_unique_ids = run_tracker(match_radius)
    print(f"✓ Tracker run: {total_unique_ids} unique junction IDs")

    # Compute pairings for every frame
    for t, fd in enumerate(frame_data):
        refresh_pairings(fd, zone_states_all[t], origin_capacity)
        print(f"  Frame {t}: {len(fd['pairings'])} pairings, "
              f"{len(fd['unpaired_tips'])} unpaired")

    # Build the human-readable tip labels
    build_tip_labels_from_pairings()

    current_phase = PHASE_ANALYZE
    _sync_panel_visibility()
    render()


def _sync_panel_visibility():
    """
    Show/hide right-panel widgets depending on the current phase.
    Phase 1 widgets (e.g. copy button) are hidden in Phase 2 and vice versa.
    """
    if current_phase == PHASE_REVIEW:
        finalize_btn.config(state='normal')
        back_btn.config(state='disabled')
        phase_header.config(text="PHASE 1 — Junction Review", bg='#1a3a6e')
        for w in analyze_only_widgets:
            w.pack_forget()
        for w in review_only_widgets:
            w.pack(**review_only_pack[w])
    else:
        finalize_btn.config(state='disabled')
        back_btn.config(state='normal')
        phase_header.config(text="PHASE 2 — Analysis", bg='#1b5e20')
        for w in analyze_only_widgets:
            w.pack(**analyze_only_pack[w])
        for w in review_only_widgets:
            w.pack_forget()


# ==============================================================================
# ███  SECTION 29 — HOVER DETECTION  ███
# ==============================================================================

def point_to_seg_dist(mouse, a, b):
    """
    Perpendicular distance from point 'mouse' to the line segment from 'a' to 'b'.

    This is used to decide whether the mouse is "hovering over" a path segment.
    The mouse and segment points are all (row, col) arrays.
    """
    ab = b - a
    ln = float(np.dot(ab, ab))      # Squared length of segment
    if ln < 1e-9:
        return float(np.linalg.norm(mouse - a))   # Degenerate segment (zero length)
    # Project mouse onto the segment, clamping to [0,1]
    t  = np.clip(np.dot(mouse-a, ab)/ln, 0.0, 1.0)
    return float(np.linalg.norm(mouse - (a + t*ab)))


def find_hovered_path(xdata, ydata):
    """
    Return the index of the pairing path currently under the mouse cursor,
    or None if no path is close enough.

    'xdata' and 'ydata' are matplotlib data coordinates (x=col, y=row).
    """
    if xdata is None or ydata is None or current_phase != PHASE_ANALYZE:
        return None
    mouse    = np.array([ydata, xdata], dtype=float)   # Convert to (row, col)
    pairings = frame_data[current_frame]['pairings']
    best_d, best_i = HOVER_RADIUS_PX, None

    for i, seg in enumerate(pairings):
        path = seg['path']
        for j in range(len(path)-1):
            a = np.array(path[j],   dtype=float)
            b = np.array(path[j+1], dtype=float)
            d = point_to_seg_dist(mouse, a, b)
            if d < best_d:
                best_d, best_i = d, i
    return best_i


def find_hovered_junc(xdata, ydata):
    """
    Return the index of the junction centroid currently under the mouse cursor,
    or None if no junction is close enough.
    """
    if xdata is None or ydata is None:
        return None
    cents   = frame_data[current_frame]['centroids']
    best_d, best_i = JUNCTION_HOVER_PX, None
    for i, (cr, cc) in enumerate(cents):
        d = math.hypot(ydata - cr, xdata - cc)   # ydata = row, xdata = col
        if d < best_d:
            best_d, best_i = d, i
    return best_i


# ==============================================================================
# ███  SECTION 30 — ZONE AND TIP CLICK DETECTION  ███
# ==============================================================================

def find_clicked_zone(xdata, ydata):
    """
    Return the zone index of the junction centroid nearest to a mouse click,
    or None if none is close enough.
    """
    if xdata is None or ydata is None:
        return None
    zones = frame_data[current_frame]['zones']
    best_d, best_zid = CLICK_ZONE_RADIUS_PX, None
    for zid, zone in enumerate(zones):
        cy, cx = zone_centroid(zone)
        d = math.hypot(xdata - cx, ydata - cy)
        if d < best_d:
            best_d, best_zid = d, zid
    return best_zid


def find_clicked_tip(xdata, ydata):
    """
    Return the (row, col) of the real tip pixel nearest to a mouse click,
    or None if none is close enough.

    NOTE: This searches ONLY real tips.  Faux tips are stored separately
    and searched by find_clicked_faux_or_real_tip() when needed.
    """
    if xdata is None or ydata is None:
        return None
    tips   = frame_data[current_frame]['tips']
    best_d, best_tip = CLICK_TIP_RADIUS_PX, None
    for tip in tips:
        d = math.hypot(xdata - tip[1], ydata - tip[0])
        if d < best_d:
            best_d, best_tip = d, tip
    return best_tip


def find_clicked_any_tip(xdata, ydata):
    """
    Return the (row, col) of the nearest tip (real OR faux) near the click,
    or None if none found.

    Used in override mode so you can override the path of a faux tip too.
    """
    if xdata is None or ydata is None:
        return None
    fd    = frame_data[current_frame]
    tips  = list(fd['tips']) + list(fd.get('faux_tips', {}).keys())
    best_d, best_tip = CLICK_TIP_RADIUS_PX, None
    for tip in tips:
        d = math.hypot(xdata - tip[1], ydata - tip[0])
        if d < best_d:
            best_d, best_tip = d, tip
    return best_tip


# ==============================================================================
# ███  SECTION 31 — MANUAL OVERRIDE LOGIC  ███
# ==============================================================================
# The "Override Path" feature lets you manually choose which junction a tip
# connects to, instead of relying on the shortest-path algorithm.
# ==============================================================================

def start_override_mode():
    """Enter override mode: next click selects a tip."""
    global override_mode, override_tip
    if current_phase != PHASE_ANALYZE:
        return
    override_mode = "waiting_tip"
    override_tip  = None
    override_btn.config(text="Cancel Override", bg='#b71c1c')
    render()


def cancel_override_mode():
    """Exit override mode without making any change."""
    global override_mode, override_tip
    override_mode = OVERRIDE_MODE_NONE
    override_tip  = None
    override_btn.config(text="Override Path (O)", bg='#f57f17')
    render()


def toggle_override_mode():
    """Toggle between override mode on and off."""
    if override_mode == OVERRIDE_MODE_NONE:
        start_override_mode()
    else:
        cancel_override_mode()


def clear_override_for_tip(tip):
    """
    Remove the manual override for a specific tip and re-run automatic pairing.
    The tip will then be handled by the greedy algorithm again.
    """
    fd = frame_data[current_frame]
    if tip in fd.get('manual_overrides', {}):
        del fd['manual_overrides'][tip]
        refresh_pairings(fd, zone_states_all[current_frame], origin_capacity)
        build_tip_labels_from_pairings()
        render()


def apply_manual_override(tip, zone_idx):
    """
    Route 'tip' to a specific junction zone (or the origin) via Dijkstra,
    store the result as a manual override, and re-run pairing for ALL other tips.

    Parameters:
        tip      — (row, col) tip pixel
        zone_idx — junction zone index, OR ORIGIN_ZONE_ID (-1) for the origin
    """
    global override_mode, override_tip
    fd     = frame_data[current_frame]
    graph  = fd['graph']
    zones  = fd['zones']
    j_ids  = fd['junction_ids']
    zs     = zone_states_all[current_frame]

    # ── Route to origin ────────────────────────────────────────────────────────
    if zone_idx == ORIGIN_ZONE_ID:
        origin_px = fd['origin_pixel']
        if origin_px is None:
            messagebox.showwarning("Override", "No origin pixel found for this frame.")
            cancel_override_mode()
            return
        path, dist = dijkstra_direct(graph, tip, origin_px)
        if path is None:
            messagebox.showwarning("Override",
                f"No path found from tip {tip} to origin {origin_px}.\n"
                "They may be disconnected in the skeleton.")
            cancel_override_mode()
            return

        seg = {
            'tip'          : tip,
            'zone_id'      : ORIGIN_ZONE_ID,
            'zone_entry_px': origin_px,
            'dist_px'      : dist,
            'path'         : path,
            'manual'       : True,
            'is_faux'      : tip in fd.get('faux_tips', {}),
        }
        fd['manual_overrides'][tip] = seg

        # Re-run all non-override pairings to prevent double-claiming
        refresh_pairings(fd, zs, origin_capacity)
        build_tip_labels_from_pairings()

        override_mode = OVERRIDE_MODE_NONE
        override_tip  = None
        override_btn.config(text="Override Path (O)", bg='#f57f17')
        render()
        print(f"  Manual override: tip {tip} -> origin {origin_px}, dist={dist:.1f}px")
        return

    # ── Route to junction ──────────────────────────────────────────────────────
    if zone_idx >= len(zones):
        cancel_override_mode()
        return

    # Find the closest path to ANY pixel in the target zone
    zone_pixels = list(zones[zone_idx])
    best_path, best_dist = None, math.inf
    for zp in zone_pixels:
        if zp not in graph:
            continue
        path, d = dijkstra_direct(graph, tip, zp)
        if path is not None and d < best_dist:
            best_dist = d
            best_path = path

    if best_path is None:
        messagebox.showwarning("Override",
            f"No path found from tip {tip} to junction {zone_idx}.\n"
            "They may be disconnected in the skeleton.")
        cancel_override_mode()
        return

    jid = j_ids[zone_idx] if zone_idx < len(j_ids) else zone_idx
    seg = {
        'tip'          : tip,
        'zone_id'      : zone_idx,
        'zone_entry_px': best_path[-1],
        'dist_px'      : best_dist,
        'path'         : best_path,
        'manual'       : True,
        'is_faux'      : tip in fd.get('faux_tips', {}),
    }
    fd['manual_overrides'][tip] = seg

    # Re-run all non-override pairings to prevent double-claiming
    refresh_pairings(fd, zs, origin_capacity)
    build_tip_labels_from_pairings()

    override_mode = OVERRIDE_MODE_NONE
    override_tip  = None
    override_btn.config(text="Override Path (O)", bg='#f57f17')
    render()
    print(f"  Manual override: tip {tip} -> J{jid} (zone {zone_idx}), "
          f"dist={best_dist:.1f}px")


# ==============================================================================
# ███  SECTION 32 — FAUX TIP LOGIC  ███
# ==============================================================================
#
# WHAT IS A FAUX TIP?
# ────────────────────
# A "faux tip" (ghost tip) is an invisible synthetic tip that lives at the
# SAME pixel location as a real tip but has a completely separate identity.
# It participates in the pairing algorithm just like a real tip, so a second
# path can depart from the same pixel to a different junction.
#
# USE CASE:
# Two branches grow on top of each other for several frames, so only ONE tip
# is visible.  Later they split.  To avoid a data gap, you add a faux tip
# to the real one during the "merged" frames, then use manual override to
# route the faux tip to the second junction.  Both paths exist simultaneously
# in the data.
#
# In the CSV export:
#   - Real tips -> normal rows
#   - Faux tips -> rows with is_faux = 1 and tip_id starting with "faux_"
#
# HOW TO USE:
#   1. Click "Add Faux Tip (F)" or press F.
#   2. Click any real tip on the canvas.
#   3. A ghost copy is created.  It appears as a gold triangle marker.
#   4. (Optional) Use Override Path to route the faux tip to a specific junction.
#   5. Use "Clear Faux Tips" to remove all faux tips for this frame.
# ==============================================================================

def start_faux_mode():
    """Enter faux-tip mode: the next tip click will create a ghost copy."""
    global faux_mode
    if current_phase != PHASE_ANALYZE:
        return
    faux_mode = FAUX_MODE_WAITING
    faux_btn.config(text="Cancel Faux (F)", bg='#b71c1c')
    render()


def cancel_faux_mode():
    """Exit faux-tip mode without creating anything."""
    global faux_mode
    faux_mode = FAUX_MODE_NONE
    faux_btn.config(text="Add Faux Tip (F)", bg='#7b1fa2')
    render()


def toggle_faux_mode():
    """Toggle between faux-tip mode on and off."""
    if faux_mode == FAUX_MODE_NONE:
        start_faux_mode()
    else:
        cancel_faux_mode()


def apply_faux_tip(real_tip):
    """
    Create a faux (ghost) tip at the same pixel as 'real_tip' for the current
    frame.

    The faux tip gets a unique ID like "faux_3" using the global counter.
    It is added to fd['faux_tips'] and then the frame's pairings are refreshed
    so the new tip participates in the greedy algorithm.

    Parameters:
        real_tip — (row, col) of the real tip to clone
    """
    global faux_tip_counter, faux_mode

    fd = frame_data[current_frame]

    # Assign a unique ID to this faux tip
    faux_tip_counter += 1
    faux_id = f"faux_{faux_tip_counter}"

    # Store the faux tip.  The pixel IS the same as the real tip.
    # We use the same (row, col) tuple as the "identity" of the faux tip.
    # But we need a KEY in fd['faux_tips'] that is distinct from the real tip key,
    # even if they're at the same pixel.
    # To solve this, we use the faux_id string as the value, and the same pixel as a key.
    # When two entries share the same pixel, refresh_pairings() treats them as
    # separate tips because it iterates over all_tips = real_tips + faux_pixels,
    # and the same pixel can appear twice (once from real_tips, once from faux_tips).
    # Dijkstra from the same start pixel will find paths equally, that's fine
    # because excluded_edges will force them to different junctions.
    fd['faux_tips'][real_tip] = faux_id

    # Also add to tip_ids so it gets a label
    if real_tip not in tip_ids[current_frame]:
        tip_ids[current_frame][real_tip] = faux_id
    # Note: build_tip_labels_from_pairings() will overwrite this if the faux tip
    # gets successfully paired,  which is what we want.

    # Refresh pairings so the new faux tip participates
    refresh_pairings(fd, zone_states_all[current_frame], origin_capacity)
    build_tip_labels_from_pairings()

    # Exit faux mode
    faux_mode = FAUX_MODE_NONE
    faux_btn.config(text="Add Faux Tip (F)", bg='#7b1fa2')
    render()
    print(f"  Faux tip created: {faux_id} @ {real_tip}")


def clear_faux_tips():
    """
    Remove ALL faux tips from the current frame and refresh pairings.
    Also removes any manual overrides that were associated with faux tips.
    """
    fd = frame_data[current_frame]

    faux_pixels = set(fd.get('faux_tips', {}).keys())

    # Remove manual overrides for faux tips
    for fp in faux_pixels:
        fd['manual_overrides'].pop(fp, None)

    # Clear the faux tips dictionary
    fd['faux_tips'].clear()

    # Refresh pairings without the faux tips
    refresh_pairings(fd, zone_states_all[current_frame], origin_capacity)
    build_tip_labels_from_pairings()
    render()
    print(f"  Cleared {len(faux_pixels)} faux tip(s) from frame {current_frame}")


# ==============================================================================
# ███  SECTION 33 — CLIPBOARD COPY  ███
# ==============================================================================

def copy_figure_to_clipboard():
    """
    Copy the current matplotlib figure to the Windows clipboard as a bitmap.

    Requires the win32clipboard and Pillow packages.
    If they're not available, saves to a temporary PNG file instead.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    buf.seek(0)
    png_data = buf.read()
    buf.close()

    try:
        import win32clipboard, win32con
        from PIL import Image
        img     = Image.open(io.BytesIO(png_data))
        bmp_buf = io.BytesIO()
        img.convert('RGB').save(bmp_buf, format='BMP')
        bmp_data = bmp_buf.getvalue()[14:]   # BMP data without the 14-byte file header
        win32clipboard.OpenClipboard()
        win32clipboard.EmptyClipboard()
        win32clipboard.SetClipboardData(win32con.CF_DIB, bmp_data)
        win32clipboard.CloseClipboard()
        print("✓ Figure copied to clipboard (Windows)")
        return
    except ImportError:
        pass   # win32clipboard or Pillow not installed, fall back to file

    # Fallback: save to a PNG file next to this script
    tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_clipboard_export.png")
    with open(tmp, 'wb') as f_:
        f_.write(png_data)
    messagebox.showinfo("Clipboard not available",
        f"win32clipboard / Pillow not found.\nFigure saved as:\n{tmp}")


# ==============================================================================
# ███  SECTION 34 — CSV EXPORT  ███
# ==============================================================================
# Exports three CSV files:
#   tip_data.csv      — one row per tip per frame with distance and speed
#   full_paths.csv    — one row per skeleton pixel per path per frame
#                       (path_type = 'geodesic' or 'untraced')
#   frame_summary.csv — one row per frame with aggregate stats
# ==============================================================================

def collect_untraced_pixels(fd):
    """
    Find skeleton pixels that are NOT part of any pairing path in this frame.

    These are loops, dead-end spurs, disconnected strands, etc., parts of the
    skeleton that the tip-to-junction routing algorithm never traversed.

    Returns a set of (row, col) pixel tuples.
    """
    skel     = fd['skel']
    pairings = fd['pairings']

    # All skeleton pixels in this frame
    all_skel_pixels = set(zip(*np.where(skel)))

    # All pixels covered by at least one pairing path
    covered = set()
    for seg in pairings:
        covered.update(seg['path'])

    return all_skel_pixels - covered   # Uncovered = untraced


def export_csv():
    """
    Export all analysis data to three CSV files chosen by the user via a
    file-save dialog.

    The user picks a base filename; the three files are automatically named:
        <base>_tip_data.csv
        <base>_frame_summary.csv
        <base>_full_paths.csv
    """
    if current_phase != PHASE_ANALYZE:
        messagebox.showwarning("Export", "Please finalize analysis first (Phase 2).")
        return

    out_path = filedialog.asksaveasfilename(
        title="Save CSV as …",
        defaultextension=".csv",
        filetypes=[("CSV files","*.csv"), ("All files","*.*")])
    if not out_path:
        return   # User cancelled

    base      = out_path.rsplit('.', 1)[0]
    tip_csv   = base + "_tip_data.csv"
    frame_csv = base + "_frame_summary.csv"
    path_csv  = base + "_full_paths.csv"

    # ── tip_data.csv ────────────────────────────────────────────────────────────
    # One row per tip per frame.  Includes:
    #   - Position, paired junction, path length, speed, direction
    #   - is_manual_override: 1 if this pairing was set manually
    #   - is_faux: 1 if this is a faux/ghost tip
    with open(tip_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['tip_id','frame','tip_x_px','tip_y_px',
                    'paired_junction','junction_global_id',
                    'path_length_px','path_length_um',
                    'speed_px_per_frame','speed_um_per_frame','direction_deg',
                    'is_manual_override','is_faux'])

        label_prev = {}   # label -> (previous_frame, prev_x, prev_y) for speed calc

        for t in range(num_frames):
            fd         = frame_data[t]
            j_ids      = fd['junction_ids']
            t_ids      = tip_ids[t]
            faux_tips  = fd.get('faux_tips', {})
            paired_map = {seg['tip']: seg for seg in fd['pairings']}

            # Export real tips AND faux tips
            all_tips_this_frame = list(fd['tips']) + list(faux_tips.keys())

            for tip in all_tips_this_frame:
                label   = t_ids.get(tip)
                if not label:
                    continue   # No label = not tracked, skip

                tx, ty   = tip[1], tip[0]          # x = column, y = row
                is_faux  = int(tip in faux_tips)
                seg      = paired_map.get(tip)

                if seg is None:
                    # Unpaired tip
                    junc_label = 'unpaired'; junc_gid = ''; path_px = ''; path_um = ''
                    is_manual  = False
                elif seg['zone_id'] == ORIGIN_ZONE_ID:
                    junc_label = label      # e.g. "orig1"
                    junc_gid   = 'origin'
                    path_px    = f"{seg['dist_px']:.2f}"
                    path_um    = f"{seg['dist_px']*PIXEL_SIZE_UM:.2f}"
                    is_manual  = seg.get('manual', False)
                else:
                    junc_label = str(seg['zone_id'])
                    junc_gid   = str(j_ids[seg['zone_id']]) if seg['zone_id'] < len(j_ids) else ''
                    path_px    = f"{seg['dist_px']:.2f}"
                    path_um    = f"{seg['dist_px']*PIXEL_SIZE_UM:.2f}"
                    is_manual  = seg.get('manual', False)

                # Calculate speed (pixels per frame) from previous labelled position
                speed_px = speed_um = direction = ''
                if label in label_prev:
                    pf, px_, py_ = label_prev[label]
                    if pf == t - 1:   # Only calculate speed for consecutive frames
                        dx = tx - px_; dy = ty - py_
                        sp = math.hypot(dx, dy)
                        speed_px  = f"{sp:.4f}"
                        speed_um  = f"{sp*PIXEL_SIZE_UM:.4f}"
                        direction = f"{math.degrees(math.atan2(-dy, dx)):.2f}"
                label_prev[label] = (t, tx, ty)

                w.writerow([label, t, tx, ty, junc_label, junc_gid,
                            path_px, path_um, speed_px, speed_um, direction,
                            int(is_manual), is_faux])

    # ── full_paths.csv ──────────────────────────────────────────────────────────
    # One row per skeleton pixel per pairing path per frame.
    # Also includes untraced pixels (loops, spurs, etc.) with path_type='untraced'.
    with open(path_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['tip_id','frame','step_index',
                    'pixel_x','pixel_y',
                    'pixel_x_um','pixel_y_um',
                    'path_type',
                    'is_manual_override',
                    'is_faux'])

        for t in range(num_frames):
            fd        = frame_data[t]
            t_ids     = tip_ids[t]
            faux_tips = fd.get('faux_tips', {})

            # Geodesic (measured) paths
            for seg in fd['pairings']:
                tip       = seg['tip']
                label     = t_ids.get(tip, '')
                is_manual = int(seg.get('manual', False))
                is_faux   = int(tip in faux_tips)

                for step, (py, px) in enumerate(seg['path']):
                    w.writerow([label, t, step,
                                px, py,
                                f"{px*PIXEL_SIZE_UM:.4f}",
                                f"{py*PIXEL_SIZE_UM:.4f}",
                                'geodesic',
                                is_manual,
                                is_faux])

            # Untraced pixels (parts of skeleton not on any pairing path)
            untraced = collect_untraced_pixels(fd)
            for (py, px) in sorted(untraced):
                w.writerow(['', t, -1,
                            px, py,
                            f"{px*PIXEL_SIZE_UM:.4f}",
                            f"{py*PIXEL_SIZE_UM:.4f}",
                            'untraced',
                            0, 0])

    # ── frame_summary.csv ───────────────────────────────────────────────────────
    # One row per frame with aggregate counts.
    with open(frame_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['frame','n_tips','n_faux_tips','n_junctions_normal',
                    'n_junctions_crossover','n_junctions_removed',
                    'skeleton_area_px2','skeleton_area_um2',
                    'n_geodesic_pixels','n_untraced_pixels'])

        for t in range(num_frames):
            fd    = frame_data[t]
            zs    = zone_states_all[t]
            n_zones = len(fd['zones'])
            n_norm = n_zones - len(zs); n_cross = 0; n_rem = 0
            for s in zs.values():
                if   s == ZONE_STATE_NORMAL:                       n_norm  += 1
                elif s in (ZONE_STATE_CROSS, ZONE_STATE_PARALLEL): n_cross += 1
                elif s == ZONE_STATE_REMOVED:                      n_rem   += 1
            area_px2  = fd['skeleton_area']
            area_um2  = area_px2 * (PIXEL_SIZE_UM ** 2)
            n_geo     = sum(len(seg['path']) for seg in fd['pairings'])
            n_untr    = len(collect_untraced_pixels(fd))
            n_faux    = len(fd.get('faux_tips', {}))
            w.writerow([t, len(fd['tips']), n_faux, n_norm, n_cross, n_rem,
                        area_px2, f"{area_um2:.3f}", n_geo, n_untr])

    messagebox.showinfo("Export complete",
        f"Tip data       ->  {os.path.basename(tip_csv)}\n"
        f"Full paths     ->  {os.path.basename(path_csv)}\n"
        f"Frame summary  ->  {os.path.basename(frame_csv)}\n\n"
        "path_type='geodesic'  -> measured growth path\n"
        "path_type='untraced'  -> loop/spur/disconnected strand\n\n"
        "Origin branches labelled orig1, orig2, …\n"
        "junction_global_id='origin' for all origin rows.")


# ==============================================================================
# ███  SECTION 35 — FLASH NOTIFICATION  ███
# ==============================================================================
# A brief "✓ Re-tracked" banner shown on the canvas for ~1.5 seconds after
# the junction tracker is re-run.
# ==============================================================================

def start_flash(n_unique):
    """Start the flash notification showing the re-track result."""
    global total_unique_ids, flash_counter
    total_unique_ids = n_unique
    flash_counter    = 6   # Will tick down 6×250ms = 1.5 seconds
    tick_flash()


def tick_flash():
    """Decrement the flash counter and redraw.  Stops when counter hits 0."""
    global flash_counter
    if flash_counter > 0:
        flash_counter -= 1
        render()
        root.after(250, tick_flash)   # Schedule next tick in 250 milliseconds
    else:
        render()   # Final redraw with no banner


# ==============================================================================
# ███  SECTION 36 — GUI LAYOUT  ███
# ==============================================================================
# Build the Tkinter window with the matplotlib canvas on the left and a
# scrollable control panel on the right.
# ==============================================================================

root = tk.Tk()
root.title("Mycelial Path Viewer & Export Tool")
root.geometry("1900x1050")
root.minsize(1400, 850)

# Tkinter BooleanVar objects let checkboxes automatically update when changed.
show_overlay_var = tk.BooleanVar(value=show_overlay_master)
show_circles_var = tk.BooleanVar(value=True)

# ── Matplotlib figure embedded in the Tkinter window ─────────────────────────
fig, ax  = plt.subplots(figsize=(10, 8))
canvas   = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(side='left', fill='both', expand=True)

# ── Right panel: a scrollable frame ──────────────────────────────────────────
# We need a scrollable panel because all the controls might not fit vertically.
right_outer = tk.Frame(root, width=440)
right_outer.pack(side='right', fill='y')
right_outer.pack_propagate(False)   # Don't let child widgets resize this frame

right_vscroll = tk.Scrollbar(right_outer, orient='vertical')
right_vscroll.pack(side='right', fill='y')

right_canvas = tk.Canvas(right_outer, width=420,
                          yscrollcommand=right_vscroll.set,
                          highlightthickness=0)
right_canvas.pack(side='left', fill='both', expand=True)
right_vscroll.config(command=right_canvas.yview)

right = tk.Frame(right_canvas)   # The actual content frame inside the canvas
right_canvas.create_window((0, 0), window=right, anchor='nw', width=420)

def _on_right_configure(event):
    """Update the scroll region whenever the right panel's contents resize."""
    right_canvas.configure(scrollregion=right_canvas.bbox('all'))

right.bind('<Configure>', _on_right_configure)

def _on_right_scroll(event):
    """Forward mouse-wheel events to the right panel's scrollbar."""
    if event.num == 4:
        right_canvas.yview_scroll(-1, 'units')
    elif event.num == 5:
        right_canvas.yview_scroll(1, 'units')
    else:
        right_canvas.yview_scroll(int(-event.delta / 120), 'units')

# Bind scroll events on both the canvas and the inner frame
right_canvas.bind('<MouseWheel>', _on_right_scroll)
right_canvas.bind('<Button-4>',   _on_right_scroll)
right_canvas.bind('<Button-5>',   _on_right_scroll)
right.bind('<MouseWheel>', _on_right_scroll)
right.bind('<Button-4>',   _on_right_scroll)
right.bind('<Button-5>',   _on_right_scroll)

# ── Phase header label ────────────────────────────────────────────────────────
phase_header = tk.Label(right,
    text="PHASE 1 — Junction Review",
    font=("Arial", 11, "bold"), fg='white', bg='#1a3a6e',
    pady=6, padx=6)
phase_header.pack(fill='x', pady=(4,2))

# ── Phase transition buttons ──────────────────────────────────────────────────
phase_btn_frame = tk.Frame(right)
phase_btn_frame.pack(fill='x', pady=4, padx=6)

finalize_btn = tk.Button(
    phase_btn_frame,
    text="Finalize & Track ▶",
    command=enter_phase_analyze,
    bg='#1b5e20', fg='white',
    activebackground='#2e7d32',
    font=("Arial", 10, "bold"),
    relief='raised', bd=2)
finalize_btn.pack(fill='x', pady=(0,2))

back_btn = tk.Button(
    phase_btn_frame,
    text="◀ Back to Review",
    command=enter_phase_review,
    bg='#4a148c', fg='white',
    activebackground='#6a1b9a',
    font=("Arial", 10, "bold"),
    state='disabled',
    relief='raised', bd=2)
back_btn.pack(fill='x')

# ── Quick Help text box ───────────────────────────────────────────────────────
leg = tk.LabelFrame(right, text="Quick Help", font=("Arial", 8))
leg.pack(pady=4, padx=4, fill="x")
tk.Label(leg,
    text=(
        "PHASE 1:\n"
        "  L-click junction -> cycle state\n"
        "  (normal->cross×->cross∥->removed)\n"
        "  State change propagates forward.\n"
        "  Then click 'Finalize & Track ▶'\n\n"
        "PHASE 2:\n"
        "  Hover path -> highlight + distance\n"
        "  Override Path (O): click tip ->\n"
        "    click junction OR origin ✛\n"
        "  Add Faux Tip (F): click real tip ->\n"
        "    creates ghost copy at same pixel\n"
        "    gold triangle = faux tip\n"
        "    use Override to route it elsewhere\n"
        "  R -> refresh   O -> override\n"
        "  F -> faux tip  Ctrl+C -> copy\n"
        "  Scroll -> zoom   Shift+drag -> pan\n\n"
        "Origins labelled orig1, orig2, …\n"
        "Untraced skel pixels in full_paths.csv\n"
        "  path_type='untraced' (loops/spurs)\n"
        "⚠ = tips ≠ active junctions\n"
        "[M] = manual override  [F] = faux tip"
    ),
    justify="left", anchor="w", font=("Courier", 7)
).pack(anchor="w", padx=4, pady=2)

# ── Match radius slider ───────────────────────────────────────────────────────
rad_frame = tk.LabelFrame(right, text="Junction ID matching", font=("Arial", 8))
rad_frame.pack(pady=6, padx=6, fill='x')

tk.Label(rad_frame, text="Match radius (px)", font=("Courier", 8)).pack(anchor='w', padx=4)

def _on_slider_move(v):
    """Called whenever the match-radius slider moves.  Updates match_radius and redraws."""
    global match_radius
    match_radius = int(float(v))
    render()

match_slider = tk.Scale(rad_frame, from_=3, to=120, resolution=1,
                         orient='horizontal', length=175, font=("Courier", 7),
                         command=_on_slider_move)
match_slider.set(DEFAULT_MATCH_RADIUS)
match_slider.pack(padx=4)

def do_retrack():
    """Re-run the junction tracker with the current match radius and flash the result."""
    n_unique = run_tracker(match_radius)
    if current_phase == PHASE_ANALYZE:
        build_tip_labels_from_pairings()
    start_flash(n_unique)

tk.Button(rad_frame, text="↺  Re-track with this radius",
          command=do_retrack,
          bg='#2d5a3d', fg='white', activebackground='#3d8a5d',
          font=("Arial", 9, "bold"), relief='raised', bd=2
          ).pack(pady=5, padx=4, fill='x')

tk.Checkbutton(rad_frame, text="Show radius circles",
               variable=show_circles_var,
               command=lambda: render(),
               font=("Courier", 8)).pack(anchor='w', padx=4, pady=(0,4))

# ── Dilation slider ───────────────────────────────────────────────────────────
tk.Label(right, text="Dilation (cosmetic):", font=("Arial", 9, "bold")).pack(pady=(10,0))
dil_frame = tk.Frame(right)
dil_frame.pack()
dil_var   = tk.IntVar(value=0)

def on_dilation_change(val):
    """Called whenever the dilation slider moves.  Updates dilation_radius and redraws."""
    global dilation_radius
    dilation_radius = int(float(val))
    render()

dil_slider = tk.Scale(dil_frame, from_=0, to=8, orient='horizontal',
                       variable=dil_var, command=on_dilation_change,
                       length=160, showvalue=True, font=("Courier", 8),
                       label="radius (px)")
dil_slider.pack(padx=4, pady=2)

# ── Frame navigation buttons ──────────────────────────────────────────────────
tk.Label(right, text="Frame Navigation:", font=("Arial", 9, "bold")).pack(pady=(8,0))
nav = tk.Frame(right)
nav.pack(pady=4)

def go_prev():
    """Go to the previous frame."""
    global current_frame, hovered_path, hovered_junc
    current_frame = max(0, current_frame-1)
    hovered_path = hovered_junc = None
    render()

def go_next():
    """Go to the next frame."""
    global current_frame, hovered_path, hovered_junc
    current_frame = min(num_frames-1, current_frame+1)
    hovered_path = hovered_junc = None
    render()

tk.Button(nav, text="◀ Prev", command=go_prev, width=8).pack(side='left', padx=3)
tk.Button(nav, text="Next ▶", command=go_next, width=8).pack(side='left', padx=3)

# ── Zoom buttons ──────────────────────────────────────────────────────────────
tk.Label(right, text="Zoom:", font=("Arial", 9, "bold")).pack(pady=(6,0))
zf = tk.Frame(right)
zf.pack(pady=4)

def zoom_in():
    """Zoom in by 30%."""
    global zoom_level; zoom_level *= 1.3; render()

def zoom_out():
    """Zoom out by 30% (minimum zoom = 0.3×)."""
    global zoom_level; zoom_level = max(0.3, zoom_level/1.3); render()

def reset_vw():
    """Reset zoom and pan to the default (full image)."""
    global zoom_level, pan_offset
    zoom_level = 1.0; pan_offset = [0.,0.]; render()

tk.Button(zf, text="＋", command=zoom_in,  width=4).pack(side='left', padx=2)
tk.Button(zf, text="－", command=zoom_out, width=4).pack(side='left', padx=2)
tk.Button(zf, text="⟳",  command=reset_vw, width=4).pack(side='left', padx=2)

# ── Info table ────────────────────────────────────────────────────────────────
tk.Label(right, text="Info:", font=("Arial", 9, "bold")).pack(pady=(6,0))
table = scrolledtext.ScrolledText(right, width=42, height=14, font=("Courier", 8))
table.pack(pady=5, padx=4)


# ==============================================================================
# ███  SECTION 37 — ANALYZE-ONLY WIDGETS  ███
# ==============================================================================
# These widgets are only visible in Phase 2. They are hidden in Phase 1 via
# _sync_panel_visibility().
# ==============================================================================

separator1 = tk.Frame(right, height=2, bg='#cccccc')

# ── Override Path button ──────────────────────────────────────────────────────
override_btn = tk.Button(right, text="Override Path (O)",
                          command=toggle_override_mode,
                          width=22, bg='#f57f17', fg='white',
                          activebackground='#e65100',
                          font=("Arial", 9, "bold"),
                          relief='raised', bd=2)

# ── Add Faux Tip button  ──────────────────────────────────────────────────────
faux_btn = tk.Button(right, text="Add Faux Tip (F)",
                      command=toggle_faux_mode,
                      width=22, bg='#7b1fa2', fg='white',
                      activebackground='#4a0072',
                      font=("Arial", 9, "bold"),
                      relief='raised', bd=2)

# ── Clear Faux Tips button ────────────────────────────────────────────────────
clear_faux_btn = tk.Button(right, text="Clear Faux Tips (this frame)",
                            command=clear_faux_tips,
                            width=22, bg='#4a148c', fg='white',
                            activebackground='#2d0057',
                            font=("Arial", 8),
                            relief='raised', bd=1)

# ── Clear Override button ─────────────────────────────────────────────────────
def clear_hovered_override():
    """
    Remove a manual override from the current frame.
    If there is only one, removes it automatically.
    If there are multiple, asks the user which tip to clear.
    """
    fd     = frame_data[current_frame]
    manual = fd.get('manual_overrides', {})
    if not manual:
        messagebox.showinfo("Clear Override", "No manual overrides on this frame.")
        return
    if len(manual) == 1:
        tip = list(manual.keys())[0]
    else:
        # Ask user to pick which tip to clear
        options = [f"{tip_ids[current_frame].get(t,'?')} {t}" for t in manual]
        choice  = simpledialog.askstring(
            "Clear Override",
            "Multiple overrides. Enter tip coords (y, x) to clear:\n" +
            "\n".join(options),
            parent=root)
        if not choice:
            return
        tip = None
        for t in manual:
            if str(t) in choice or tip_ids[current_frame].get(t,'') in choice:
                tip = t
                break
        if tip is None:
            messagebox.showwarning("Clear Override", "Tip not found.")
            return
    clear_override_for_tip(tip)

clear_override_btn = tk.Button(right, text="Clear Override (this frame)",
                                command=clear_hovered_override,
                                width=22, bg='#6d4c41', fg='white',
                                activebackground='#4e342e',
                                font=("Arial", 8),
                                relief='raised', bd=1)

# ── Origin capacity control ───────────────────────────────────────────────────
origin_cap_frame = tk.LabelFrame(right, text="Origin connections", font=("Arial", 8))
cap_frame2 = tk.Frame(origin_cap_frame)
cap_frame2.pack()
cap_var    = tk.IntVar(value=origin_capacity)
cap_label2 = tk.Label(cap_frame2, textvariable=cap_var, width=4,
                       font=("Courier", 12, "bold"), relief='sunken')
cap_label2.pack(side='left', padx=4)

def cap_inc():
    """Increase the number of tips allowed to pair to the origin by 1."""
    global origin_capacity
    origin_capacity += 1
    cap_var.set(origin_capacity)
    _recompute_all_and_render()

def cap_dec():
    """Decrease the number of tips allowed to pair to the origin by 1 (minimum 0)."""
    global origin_capacity
    origin_capacity = max(0, origin_capacity-1)
    cap_var.set(origin_capacity)
    _recompute_all_and_render()

def _recompute_all_and_render():
    """Re-run pairing for all frames and redraw.  Called when origin capacity changes."""
    global hovered_path
    hovered_path = None
    for t, fd in enumerate(frame_data):
        refresh_pairings(fd, zone_states_all[t], origin_capacity)
    build_tip_labels_from_pairings()
    render()

tk.Button(cap_frame2, text="＋", command=cap_inc, width=3).pack(side='left', padx=2)
tk.Button(cap_frame2, text="－", command=cap_dec, width=3).pack(side='left', padx=2)

# ── Refresh / Reset / Copy / Export buttons ───────────────────────────────────
refresh_btn_w   = tk.Button(right, text="Refresh Frame  (R)",
                             command=lambda: _do_refresh(), width=22, bg='#c8e6c9')
reset_zones_btn = tk.Button(right, text="Reset Zone States (this frame)",
                             command=lambda: _reset_zone_states(), width=22, bg='#ffe0b2')
copy_btn_w      = tk.Button(right, text="Copy View  (Ctrl+C)",
                             command=copy_figure_to_clipboard, width=22, bg='#e8eaf6')
export_btn_w    = tk.Button(right, text="Export CSV …",
                             command=export_csv, width=22, bg='#bbdefb',
                             font=("Arial", 9, "bold"))

def _do_refresh():
    """Refresh pairings for the current frame and redraw."""
    global hovered_path
    hovered_path = None
    refresh_pairings(frame_data[current_frame],
                     zone_states_all[current_frame], origin_capacity)
    build_tip_labels_from_pairings()
    render()

def _reset_zone_states():
    """Reset all zone states and manual overrides for the current frame, then redraw."""
    global hovered_path
    zone_states_all[current_frame].clear()
    frame_data[current_frame]['manual_overrides'].clear()
    hovered_path = None
    refresh_pairings(frame_data[current_frame],
                     zone_states_all[current_frame], origin_capacity)
    build_tip_labels_from_pairings()
    render()

# ── Widget show/hide configuration ────────────────────────────────────────────
# analyze_only_widgets: shown ONLY in Phase 2
# review_only_widgets:  shown ONLY in Phase 1
# The dicts map each widget to the keyword arguments for .pack()

analyze_only_widgets = [
    separator1,
    override_btn,
    faux_btn,
    clear_faux_btn,
    clear_override_btn,
    origin_cap_frame,
    refresh_btn_w,
    reset_zones_btn,
    copy_btn_w,
    export_btn_w,
]
analyze_only_pack = {
    separator1        : dict(fill='x', pady=6, padx=6),
    override_btn      : dict(pady=3, padx=6, fill='x'),
    faux_btn          : dict(pady=3, padx=6, fill='x'),
    clear_faux_btn    : dict(pady=2, padx=6, fill='x'),
    clear_override_btn: dict(pady=2, padx=6, fill='x'),
    origin_cap_frame  : dict(pady=6, padx=6, fill='x'),
    refresh_btn_w     : dict(pady=3),
    reset_zones_btn   : dict(pady=2),
    copy_btn_w        : dict(pady=2),
    export_btn_w      : dict(pady=6),
}

copy_review_btn  = tk.Button(right, text="Copy View  (Ctrl+C)",
                              command=copy_figure_to_clipboard,
                              width=22, bg='#e8eaf6')
review_only_widgets = [copy_review_btn]
review_only_pack    = {copy_review_btn: dict(pady=2)}

copy_review_btn.pack(pady=2)   # Visible immediately (Phase 1 is the start)

tk.Label(right, text="", height=2).pack()   # Spacer at the bottom


# ==============================================================================
# ███  SECTION 38 — EVENT HANDLERS  ███
# ==============================================================================
# These functions are called by matplotlib when the user does something with
# the mouse or keyboard on the canvas.
# ==============================================================================

def on_scroll(event):
    """
    Mouse-wheel scroll over the canvas: zoom in or out.
    Only responds to scroll events that are INSIDE the axes (the image area).
    """
    global zoom_level
    if event.inaxes == ax:
        if event.button == 'up':
            zoom_level *= 1.12    # Zoom in 12%
        elif event.button == 'down':
            zoom_level = max(0.3, zoom_level / 1.12)   # Zoom out 12% (floor = 0.3)
        render()


def on_motion(event):
    """
    Mouse movement over the canvas: update hover state and re-render if needed.
    Also handles panning (shift+drag or middle-button drag).
    """
    global hovered_path, hovered_junc, pan_offset, pan_start

    # If mouse is outside the image axes, clear hover
    if event.inaxes != ax:
        if hovered_path is not None or hovered_junc is not None:
            hovered_path = hovered_junc = None
            render()
        return

    # If panning, update the pan offset based on how far the mouse moved
    if is_panning:
        if event.xdata is not None:
            pan_offset[0] -= event.xdata - pan_start[0]
            pan_offset[1] -= event.ydata - pan_start[1]
            pan_start[0]   = event.xdata
            pan_start[1]   = event.ydata
            render()
        return

    # Check if mouse is near a junction centroid (takes priority over path hover)
    new_junc = find_hovered_junc(event.xdata, event.ydata)

    # Only check path hover if no junction is being hovered and overlay is visible
    new_path = (find_hovered_path(event.xdata, event.ydata)
                if show_overlay_var.get() and new_junc is None else None)

    # Only redraw if the hover state actually changed (avoid unnecessary redraws)
    if new_junc != hovered_junc or new_path != hovered_path:
        hovered_junc = new_junc
        hovered_path = new_path
        render()


def on_press(event):
    """
    Mouse button press on the canvas.

    Handles:
      - Middle-button / Shift+Left-button -> start panning
      - Left-click in override mode       -> select tip, then target
      - Left-click in faux-tip mode       -> create ghost tip
      - Left-click normally               -> cycle junction state
      - Right-click in Phase 2            -> rename tip label (or cancel override)
    """
    global is_panning, pan_start, override_mode, override_tip, faux_mode

    if event.inaxes != ax or event.xdata is None:
        return   # Click was outside the image therefore ignore it

    # ── Start panning ─────────────────────────────────────────────────────────
    if event.button == 2 or (event.button == 1 and event.key == 'shift'):
        is_panning   = True
        pan_start[0] = event.xdata
        pan_start[1] = event.ydata
        return

    # ── Left-click in Phase 2 ─────────────────────────────────────────────────
    if event.button == 1 and current_phase == PHASE_ANALYZE:

        # ── Faux-tip mode: waiting for user to click a real tip ───────────────
        if faux_mode == FAUX_MODE_WAITING:
            tip = find_clicked_tip(event.xdata, event.ydata)   # Real tips only
            if tip is not None:
                apply_faux_tip(tip)   # Create the ghost copy
            else:
                cancel_faux_mode()    # Clicked somewhere else so cancel
            return

        # ── Override mode: waiting for tip selection ──────────────────────────
        if override_mode == "waiting_tip":
            tip = find_clicked_any_tip(event.xdata, event.ydata)   # Real + faux
            if tip is not None:
                override_tip  = tip
                override_mode = OVERRIDE_MODE_TIP
                render()
            return

        # ── Override mode: tip selected, waiting for target click ─────────────
        if override_mode == OVERRIDE_MODE_TIP:
            # Check origin first, as it's a valid override target
            if find_clicked_origin(event.xdata, event.ydata):
                apply_manual_override(override_tip, ORIGIN_ZONE_ID)
                return
            # Otherwise try a junction centroid
            zid = find_clicked_zone(event.xdata, event.ydata)
            if zid is not None and override_tip is not None:
                apply_manual_override(override_tip, zid)
            return

    # ── Normal left-click: cycle junction state ───────────────────────────────
    if event.button == 1:
        zid = find_clicked_zone(event.xdata, event.ydata)
        if zid is not None:
            zs  = zone_states_all[current_frame]
            cur = zs.get(zid, ZONE_STATE_NORMAL)
            nxt = (cur + 1) % 4   # Cycle: 0->1->2->3->0->…
            zs[zid] = nxt

            # Propagate the new state forward to all subsequent frames
            if current_phase in (PHASE_REVIEW, PHASE_ANALYZE):
                propagate_zone_state_forward(current_frame, zid, nxt)

            # Re-run pairings if in Phase 2 (a state change may affect pairing)
            if current_phase == PHASE_ANALYZE:
                refresh_pairings(frame_data[current_frame], zs, origin_capacity)
                build_tip_labels_from_pairings()
            render()

    # ── Right-click in Phase 2: cancel override or rename tip label ───────────
    if event.button == 3 and current_phase == PHASE_ANALYZE:
        if override_mode != OVERRIDE_MODE_NONE:
            cancel_override_mode()   # Right-click = cancel override
            return
        if faux_mode != FAUX_MODE_NONE:
            cancel_faux_mode()       # Right-click = cancel faux mode
            return
        
        # Otherwise: allow manual renaming of a tip's display label
        tip = find_clicked_any_tip(event.xdata, event.ydata)
        if tip is not None:
            t_ids   = tip_ids[current_frame]
            current = t_ids.get(tip, '')
            new_id  = simpledialog.askstring(
                "Tip ID",
                f"Override label for tip at {tip}\n(auto-label = junction ID):",
                initialvalue=current,
                parent=root)
            if new_id is not None:
                if new_id.strip():
                    t_ids[tip] = new_id.strip()
                else:
                    t_ids.pop(tip, None)   # Empty string = remove label
            render()


def on_release(event):
    """Mouse button released: stop panning."""
    global is_panning
    is_panning = False


def on_key(event):
    """
    Keyboard shortcuts:
      <-/->      — navigate frames
      R          — refresh pairings (Phase 2 only)
      O          — toggle override mode (Phase 2 only)
      F          — toggle faux-tip mode (Phase 2 only)
      Ctrl+C     — copy view to clipboard
      Escape     — cancel override or faux mode
    """
    global current_frame, hovered_path, hovered_junc

    if event.key == 'right':
        current_frame = min(num_frames-1, current_frame+1)
        hovered_path = hovered_junc = None
        render()
    elif event.key == 'left':
        current_frame = max(0, current_frame-1)
        hovered_path = hovered_junc = None
        render()
    elif event.key == 'r' and current_phase == PHASE_ANALYZE:
        _do_refresh()
    elif event.key == 'o' and current_phase == PHASE_ANALYZE:
        toggle_override_mode()
    elif event.key == 'f' and current_phase == PHASE_ANALYZE:
        toggle_faux_mode()
    elif event.key == 'ctrl+c':
        copy_figure_to_clipboard()
    elif event.key == 'escape':
        # Cancel whichever special mode is active
        if override_mode != OVERRIDE_MODE_NONE:
            cancel_override_mode()
        elif faux_mode != FAUX_MODE_NONE:
            cancel_faux_mode()


# ── Connect event handlers to the matplotlib canvas ──────────────────────────
fig.canvas.mpl_connect('scroll_event',         on_scroll)
fig.canvas.mpl_connect('motion_notify_event',  on_motion)
fig.canvas.mpl_connect('button_press_event',   on_press)
fig.canvas.mpl_connect('button_release_event', on_release)
fig.canvas.mpl_connect('key_press_event',      on_key)


def on_close(event=None):
    """Clean up when the user closes the window."""
    plt.close('all')
    root.quit()
    root.destroy()

root.protocol("WM_DELETE_WINDOW", on_close)


# ==============================================================================
# ███  SECTION 39 — LAUNCH  ███
# ==============================================================================
# Print a summary of controls to the terminal and start the GUI event loop.
# ==============================================================================

print("=" * 60)
print("Mycelial Path Viewer & Export Tool")
print("=" * 60)
print("PHASE 1 — Junction Review:")
print("  L-click junction -> cycle state (normal->cross×->cross∥->removed)")
print("  State change propagates to all subsequent frames with same J-ID.")
print("  Navigate all frames, mark crossovers/junk junctions.")
print("  'Finalize & Track ▶' -> runs tracker + pairing -> Phase 2")
print("")
print("PHASE 2 — Analysis:")
print("    No more two tips sharing the same junction after override.")
print("")
print("  FAUX TIPS:")
print("    'Add Faux Tip (F)' or press F:")
print("    -> Click a real tip -> creates a ghost copy (gold triangle)")
print("    -> The faux tip gets its own ID (faux_1, faux_2, …)")
print("    -> Use Override Path to route it to a different junction")
print("    -> Both the real and faux tip export as separate rows in CSVs")
print("    -> is_faux=1 column lets you filter faux tips in analysis")
print("    'Clear Faux Tips (this frame)' removes all faux tips")
print("")
print("  Override Path (O): click tip -> click junction OR origin ✛")
print("  Origin highlighted green in override mode as valid target")
print("  Hover path -> highlight + distance in title bar")
print("  R -> refresh pairings for current frame")
print("  ↺ Re-track -> rerun junction ID assignment")
print("  Export CSV -> tip_data.csv + full_paths.csv + frame_summary.csv")
print("    full_paths.csv: path_type='geodesic' or 'untraced'")
print("    Untraced = loops, spurs, disconnected strands (tip_id='')")
print("  '◀ Back to Review' -> return to Phase 1 (states preserved)")
print("")
print("Both phases:")
print("  ⚠ in title = number of tips ≠ number of active junctions")
print("  Scroll (over image) -> zoom   |   Shift+drag / middle-drag -> pan")
print("  Ctrl+C -> copy view   |   Escape -> cancel override / faux mode")
print("=" * 60 + "\n")

# Initial render, draw the first frame in Phase 1
render()

# Hand control to the Tkinter event loop.
# This call blocks until the user closes the window.
root.mainloop()
