# Maze Solver

A browser tool that solves maze images. Upload one or more mazes, press **Solve**, and the route from the entrance to the exit is drawn on the image (default line color `#BCBCBC`).

## Use it

Open `index.html` in a browser. No install or server is needed.

1. Click the upload box (or drop / paste images). You can add several images and switch between them with the thumbnails.
2. Press **Solve** (or **Solve all**).
3. Change the line color or thickness at any time; the route redraws instantly.
4. Press **Download PNG** to save the solved maze at full size.

If the entrance and exit aren't guessed correctly, press **Set start & end** and click the two openings on the maze.

## How it works

`solver.js` works on raw pixels, so it handles hexagonal, square or round mazes drawn as dark lines on a light background:

1. Threshold the image into walls and open space (Otsu threshold, adjustable under *Image tuning*).
2. Thicken the walls slightly so anti-aliased lines don't leak.
3. Take the convex hull of the walls as the maze area; open gaps on its edge are the entrance and exit.
4. Run Dijkstra between the two openings with a cost that keeps the route centered in the corridors, then smooth it.

`sample.js` holds the sample maze (`samples/hex-maze-1.webp`) that loads on start.
