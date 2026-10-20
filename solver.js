/*
 * Maze solver core. Works on raw RGBA pixels, so it handles any maze drawn as
 * dark lines on a light background (hexagonal, square, circular...).
 *
 *   solveMaze({ width, height, data }, options) -> result
 *
 * Steps:
 *   1. Threshold the image into wall / open pixels (Otsu by default).
 *   2. Thicken walls a little so anti-aliased lines have no leaks.
 *   3. Take the convex hull of the walls as the maze area.
 *   4. Openings = open pixels on the hull edge that reach the maze interior.
 *   5. Dijkstra from one opening to the other, with a cost that keeps the
 *      route in the middle of the corridors.
 */
(function (root) {
  "use strict";

  function toGray(img) {
    const { width: w, height: h, data } = img;
    const g = new Uint8ClampedArray(w * h);
    for (let i = 0, p = 0; i < w * h; i++, p += 4) {
      const a = data[p + 3] / 255;
      const lum = 0.299 * data[p] + 0.587 * data[p + 1] + 0.114 * data[p + 2];
      g[i] = a * lum + (1 - a) * 255; // transparent counts as white paper
    }
    return g;
  }

  function otsu(gray) {
    const hist = new Float64Array(256);
    for (let i = 0; i < gray.length; i++) hist[gray[i]]++;
    const total = gray.length;
    let sum = 0;
    for (let t = 0; t < 256; t++) sum += t * hist[t];
    let sumB = 0, wB = 0, best = 0, thr = 128;
    for (let t = 0; t < 256; t++) {
      wB += hist[t];
      if (!wB) continue;
      const wF = total - wB;
      if (!wF) break;
      sumB += t * hist[t];
      const mB = sumB / wB, mF = (sum - sumB) / wF;
      const between = wB * wF * (mB - mF) * (mB - mF);
      if (between > best) { best = between; thr = t; }
    }
    return thr;
  }

  // Square dilation (op = 1) or erosion (op = 0) of a binary mask, separable.
  function morph(mask, w, h, r, op) {
    if (r <= 0) return mask.slice();
    const tmp = new Uint8Array(w * h), out = new Uint8Array(w * h);
    const hit = op ? 1 : 0;
    for (let y = 0; y < h; y++) {
      const row = y * w;
      for (let x = 0; x < w; x++) {
        let v = 1 - hit;
        for (let k = Math.max(0, x - r), e = Math.min(w - 1, x + r); k <= e; k++) {
          if (mask[row + k] === hit) { v = hit; break; }
        }
        // Erosion treats outside the image as empty.
        if (!op && (x - r < 0 || x + r > w - 1)) v = 0;
        tmp[row + x] = v;
      }
    }
    for (let x = 0; x < w; x++) {
      for (let y = 0; y < h; y++) {
        let v = 1 - hit;
        for (let k = Math.max(0, y - r), e = Math.min(h - 1, y + r); k <= e; k++) {
          if (tmp[k * w + x] === hit) { v = hit; break; }
        }
        if (!op && (y - r < 0 || y + r > h - 1)) v = 0;
        out[y * w + x] = v;
      }
    }
    return out;
  }

  // Labels 8-connected components of mask; returns { labels, sizes }.
  function components(mask, w, h, conn8) {
    const labels = new Int32Array(w * h).fill(-1);
    const sizes = [];
    const stack = new Int32Array(w * h);
    for (let s = 0; s < w * h; s++) {
      if (!mask[s] || labels[s] !== -1) continue;
      const id = sizes.length;
      let top = 0, n = 0;
      stack[top++] = s;
      labels[s] = id;
      while (top) {
        const i = stack[--top];
        n++;
        const x = i % w, y = (i / w) | 0;
        for (let dy = -1; dy <= 1; dy++) {
          const ny = y + dy;
          if (ny < 0 || ny >= h) continue;
          for (let dx = -1; dx <= 1; dx++) {
            if (!dx && !dy) continue;
            if (!conn8 && dx && dy) continue;
            const nx = x + dx;
            if (nx < 0 || nx >= w) continue;
            const j = ny * w + nx;
            if (mask[j] && labels[j] === -1) { labels[j] = id; stack[top++] = j; }
          }
        }
      }
      sizes.push(n);
    }
    return { labels, sizes };
  }

  function convexHull(pts) {
    pts.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const lower = [], upper = [];
    for (const p of pts) {
      while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) lower.pop();
      lower.push(p);
    }
    for (let i = pts.length - 1; i >= 0; i--) {
      const p = pts[i];
      while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) upper.pop();
      upper.push(p);
    }
    upper.pop(); lower.pop();
    return lower.concat(upper);
  }

  // Fills the convex polygon into a mask (one span per row).
  function fillConvex(poly, w, h) {
    const minX = new Float64Array(h).fill(Infinity), maxX = new Float64Array(h).fill(-Infinity);
    for (let i = 0; i < poly.length; i++) {
      const [x0, y0] = poly[i], [x1, y1] = poly[(i + 1) % poly.length];
      const ya = Math.min(y0, y1), yb = Math.max(y0, y1);
      for (let y = Math.max(0, Math.ceil(ya)); y <= Math.min(h - 1, Math.floor(yb)); y++) {
        const x = y1 === y0 ? Math.min(x0, x1) : x0 + (x1 - x0) * (y - y0) / (y1 - y0);
        const x2 = y1 === y0 ? Math.max(x0, x1) : x;
        if (x < minX[y]) minX[y] = x;
        if (x2 > maxX[y]) maxX[y] = x2;
      }
    }
    const mask = new Uint8Array(w * h);
    for (let y = 0; y < h; y++) {
      if (minX[y] > maxX[y]) continue;
      const a = Math.max(0, Math.ceil(minX[y] - 0.5)), b = Math.min(w - 1, Math.floor(maxX[y] + 0.5));
      mask.fill(1, y * w + a, y * w + b + 1);
    }
    return mask;
  }

  // Chamfer (3-4) distance to the nearest blocked pixel, in pixels.
  function distanceMap(free, w, h) {
    const INF = 1e9;
    const d = new Float32Array(w * h);
    for (let i = 0; i < w * h; i++) d[i] = free[i] ? INF : 0;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        const i = y * w + x;
        if (!d[i]) continue;
        let v = d[i];
        if (x > 0) v = Math.min(v, d[i - 1] + 3); else v = Math.min(v, 3);
        if (y > 0) {
          v = Math.min(v, d[i - w] + 3);
          if (x > 0) v = Math.min(v, d[i - w - 1] + 4);
          if (x < w - 1) v = Math.min(v, d[i - w + 1] + 4);
        } else v = Math.min(v, 3);
        d[i] = v;
      }
    }
    for (let y = h - 1; y >= 0; y--) {
      for (let x = w - 1; x >= 0; x--) {
        const i = y * w + x;
        if (!d[i]) continue;
        let v = d[i];
        if (x < w - 1) v = Math.min(v, d[i + 1] + 3); else v = Math.min(v, 3);
        if (y < h - 1) {
          v = Math.min(v, d[i + w] + 3);
          if (x < w - 1) v = Math.min(v, d[i + w + 1] + 4);
          if (x > 0) v = Math.min(v, d[i + w - 1] + 4);
        } else v = Math.min(v, 3);
        d[i] = v;
      }
    }
    for (let i = 0; i < w * h; i++) d[i] /= 3;
    return d;
  }

  // Binary min-heap of (key, index) pairs.
  function Heap(cap) {
    this.k = new Float64Array(cap);
    this.v = new Int32Array(cap);
    this.n = 0;
  }
  Heap.prototype.push = function (key, val) {
    if (this.n === this.k.length) {
      const k = new Float64Array(this.n * 2), v = new Int32Array(this.n * 2);
      k.set(this.k); v.set(this.v); this.k = k; this.v = v;
    }
    let i = this.n++;
    const K = this.k, V = this.v;
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (K[p] <= key) break;
      K[i] = K[p]; V[i] = V[p]; i = p;
    }
    K[i] = key; V[i] = val;
  };
  Heap.prototype.pop = function () {
    const K = this.k, V = this.v;
    const topV = V[0];
    this.topKey = K[0];
    const key = K[--this.n], val = V[this.n];
    let i = 0;
    for (;;) {
      let c = 2 * i + 1;
      if (c >= this.n) break;
      if (c + 1 < this.n && K[c + 1] < K[c]) c++;
      if (K[c] >= key) break;
      K[i] = K[c]; V[i] = V[c]; i = c;
    }
    K[i] = key; V[i] = val;
    return topV;
  };

  const SQ2 = Math.SQRT2;
  const NB = [[1, 0, 1], [-1, 0, 1], [0, 1, 1], [0, -1, 1], [1, 1, SQ2], [1, -1, SQ2], [-1, 1, SQ2], [-1, -1, SQ2]];

  // Dijkstra from a set of source pixels to any target pixel.
  function shortestPath(free, dist, halfWidth, w, h, sources, isTarget) {
    const cost = new Float32Array(w * h).fill(Infinity);
    const prev = new Int32Array(w * h).fill(-1);
    const done = new Uint8Array(w * h);
    const heap = new Heap(1 << 16);
    const c2 = halfWidth * halfWidth;
    for (const s of sources) { cost[s] = 0; heap.push(0, s); }
    let end = -1;
    while (heap.n) {
      const i = heap.pop();
      if (done[i]) continue;
      done[i] = 1;
      if (isTarget(i)) { end = i; break; }
      const x = i % w, y = (i / w) | 0, ci = cost[i];
      for (const [dx, dy, len] of NB) {
        const nx = x + dx, ny = y + dy;
        if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
        const j = ny * w + nx;
        if (!free[j] || done[j]) continue;
        const d = dist[j];
        // Cheap in the middle of a corridor, expensive near walls.
        const nc = ci + len * (1 + (4 * c2) / (d * d));
        if (nc < cost[j]) { cost[j] = nc; prev[j] = i; heap.push(nc, j); }
      }
    }
    if (end < 0) return null;
    const path = [];
    for (let i = end; i !== -1; i = prev[i]) path.push([i % w, (i / w) | 0]);
    return path.reverse();
  }

  function smooth(path, win) {
    if (path.length < 3 || win < 1) return path.map(p => [p[0], p[1]]);
    const out = [];
    const n = path.length;
    for (let i = 0; i < n; i++) {
      const r = Math.min(win, i, n - 1 - i);
      let sx = 0, sy = 0;
      for (let k = i - r; k <= i + r; k++) { sx += path[k][0]; sy += path[k][1]; }
      out.push([sx / (2 * r + 1), sy / (2 * r + 1)]);
    }
    return out;
  }

  // Ramer–Douglas–Peucker, keeps the drawn polyline light.
  function simplify(pts, eps) {
    if (pts.length < 3) return pts;
    const keep = new Uint8Array(pts.length);
    keep[0] = keep[pts.length - 1] = 1;
    const stack = [[0, pts.length - 1]];
    while (stack.length) {
      const [a, b] = stack.pop();
      const [ax, ay] = pts[a], [bx, by] = pts[b];
      const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy) || 1;
      let best = -1, idx = -1;
      for (let i = a + 1; i < b; i++) {
        const d = Math.abs(dy * pts[i][0] - dx * pts[i][1] + bx * ay - by * ax) / L;
        if (d > best) { best = d; idx = i; }
      }
      if (best > eps) { keep[idx] = 1; stack.push([a, idx], [idx, b]); }
    }
    return pts.filter((_, i) => keep[i]);
  }

  function nearestFree(free, w, h, px, py) {
    px = Math.round(px); py = Math.round(py);
    for (let r = 0; r < Math.max(w, h); r++) {
      let best = -1, bd = Infinity;
      for (let y = py - r; y <= py + r; y++) {
        if (y < 0 || y >= h) continue;
        for (let x = px - r; x <= px + r; x++) {
          if (x < 0 || x >= w) continue;
          if (Math.max(Math.abs(x - px), Math.abs(y - py)) !== r) continue;
          const i = y * w + x;
          if (free[i]) {
            const d = (x - px) ** 2 + (y - py) ** 2;
            if (d < bd) { bd = d; best = i; }
          }
        }
      }
      if (best >= 0) return best;
    }
    return -1;
  }

  /**
   * options:
   *   threshold  0-255 or "auto" (default "auto")
   *   thicken    wall thickening in px (default 1)
   *   start, end optional [x, y] points; when both are set they replace the
   *              detected openings
   */
  function solveMaze(img, options) {
    const o = Object.assign({ threshold: "auto", thicken: 1 }, options || {});
    const { width: w, height: h } = img;
    const N = w * h;
    const t0 = Date.now();

    const gray = toGray(img);
    const thr = o.threshold === "auto" ? otsu(gray) : +o.threshold;
    let wall = new Uint8Array(N), dark = 0;
    for (let i = 0; i < N; i++) if (gray[i] <= thr) { wall[i] = 1; dark++; }
    // Light lines on a dark background: flip so walls are the minority.
    const inverted = dark > N / 2;
    if (inverted) for (let i = 0; i < N; i++) wall[i] ^= 1;

    // Drop specks: keep wall components at least 2% the size of the largest.
    const wc = components(wall, w, h, true);
    const biggest = Math.max(0, ...wc.sizes);
    if (!biggest) throw new Error("No maze lines found. Try adjusting the line threshold.");
    const minKeep = biggest * 0.02;
    for (let i = 0; i < N; i++) if (wall[i] && wc.sizes[wc.labels[i]] < minKeep) wall[i] = 0;

    const thick = morph(wall, w, h, Math.max(0, o.thicken | 0), 1);

    // Maze area: convex hull of the wall pixels (row extremes are enough).
    const pts = [];
    for (let y = 0; y < h; y++) {
      let a = -1, b = -1;
      for (let x = 0; x < w; x++) if (wall[y * w + x]) { if (a < 0) a = x; b = x; }
      if (a >= 0) { pts.push([a, y]); if (b !== a) pts.push([b, y]); }
    }
    const hull = convexHull(pts);
    const inside = fillConvex(hull, w, h);

    const free = new Uint8Array(N);
    for (let i = 0; i < N; i++) free[i] = inside[i] && !thick[i] ? 1 : 0;

    const dist = distanceMap(free, w, h);
    const fc = components(free, w, h, false);
    let main = 0;
    for (let k = 1; k < fc.sizes.length; k++) if (fc.sizes[k] > fc.sizes[main]) main = k;

    // Typical corridor half-width: high percentile of distance inside the maze.
    const samples = [];
    for (let i = 0; i < N; i += 7) if (fc.labels[i] === main) samples.push(dist[i]);
    samples.sort((a, b) => a - b);
    const halfWidth = Math.max(1.5, samples[Math.floor(samples.length * 0.9)] || 3);

    let sources, isTarget, openings = [];
    if (o.start && o.end) {
      const s = nearestFree(free, w, h, o.start[0], o.start[1]);
      const e = nearestFree(free, w, h, o.end[0], o.end[1]);
      if (s < 0 || e < 0) throw new Error("Couldn't place the start or end point inside the maze.");
      sources = [s];
      isTarget = i => i === e;
      openings = [[s % w, (s / w) | 0], [e % w, (e / w) | 0]];
    } else {
      // Rim = band of the hull edge; open rim pixels in the main corridor
      // network are the entrance and exit.
      const rimWidth = Math.max(2, Math.round(halfWidth * 0.25));
      const core = morph(inside, w, h, rimWidth, 0);
      const rim = new Uint8Array(N);
      for (let i = 0; i < N; i++) rim[i] = inside[i] && !core[i] && free[i] && fc.labels[i] === main ? 1 : 0;
      const rc = components(rim, w, h, true);
      const clusters = rc.sizes.map((size, id) => ({ id, size, sx: 0, sy: 0, px: [] }));
      for (let i = 0; i < N; i++) {
        if (!rim[i]) continue;
        const c = clusters[rc.labels[i]];
        c.sx += i % w; c.sy += (i / w) | 0; c.px.push(i);
      }
      const minSize = Math.max(3, halfWidth * 0.5);
      const cand = clusters.filter(c => c.size >= minSize);
      for (const c of cand) { c.cx = c.sx / c.size; c.cy = c.sy / c.size; }
      if (cand.length < 2) {
        throw new Error("Couldn't find both the entrance and the exit. Pick them by hand: click Set start & end, then click the two openings.");
      }
      // Entrance and exit: the two openings farthest apart.
      let A = cand[0], B = cand[1], bestD = -1;
      for (let i = 0; i < cand.length; i++) {
        for (let j = i + 1; j < cand.length; j++) {
          const d = (cand[i].cx - cand[j].cx) ** 2 + (cand[i].cy - cand[j].cy) ** 2;
          if (d > bestD) { bestD = d; A = cand[i]; B = cand[j]; }
        }
      }
      // Start from the upper/left one so the route reads naturally.
      if (A.cy + A.cx > B.cy + B.cx) [A, B] = [B, A];
      const targetId = B.id;
      sources = A.px;
      isTarget = i => rim[i] && rc.labels[i] === targetId;
      openings = [[A.cx, A.cy], [B.cx, B.cy]];
    }

    const raw = shortestPath(free, dist, halfWidth, w, h, sources, isTarget);
    if (!raw) throw new Error("No route connects the entrance and exit. Try lowering the wall thickening.");

    const path = simplify(smooth(raw, Math.max(1, Math.round(halfWidth * 0.6))), 0.5);
    return {
      path,
      rawLength: raw.length,
      openings,
      halfWidth,
      threshold: thr,
      inverted,
      hull,
      ms: Date.now() - t0,
    };
  }

  const api = { solveMaze };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.MazeSolver = api;
})(typeof self !== "undefined" ? self : this);
