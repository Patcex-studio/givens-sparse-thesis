# Assets

This directory contains visual assets for the thesis.

## Planned Charts

### 1. Support Growth Plot
- **Description**: A line chart showing the maximum support size (number of non-zero amplitudes) as a function of the number of Givens rotations `M`.
- **Expected shape**: Linear growth from 1 to `1+M`.
- **Comparison**: Overlay with exponential growth `2^n` for full-state simulation.

### 2. Complexity Comparison
- **Description**: Log-scale plot comparing memory requirements:
  - Full state simulation: `O(2^n)`
  - Givens sparse simulation: `O(M)`
- **Insight**: Even for moderate `n`, the sparse method is dramatically more efficient.

## Adding Charts

To generate these plots, use Python with Matplotlib or generate SVG diagrams manually. Place the final images as:
- `support_growth.png`
- `complexity_plot.png`