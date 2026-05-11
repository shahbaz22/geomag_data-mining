# Geomagnetic Data Mining

Network analysis of ULF geomagnetic pulsations (Pc1-Pc5) using ground magnetometer data from [SuperMAG](https://supermag.jhuapl.edu/).

## Motivation

Extreme space weather events pose a significant risk to modern infrastructure, with potential economic impacts estimated in the trillions of dollars. Geomagnetically induced currents (GICs) can damage power grids, disrupt satellite operations, degrade GPS accuracy, and compromise communications systems. Understanding and characterising these events is critical for forecasting and mitigation.

This project uses **spatiotemporal pattern analysis** to characterise extreme space weather events by studying how ULF (ultra-low frequency) geomagnetic pulsations propagate across networks of ground-based magnetometer stations. By identifying coherent oscillation patterns, community structures, and temporal evolution in station networks, we can better understand the large-scale dynamics of geomagnetic disturbances and work towards improved early warning systems.

## Approach

We build time-evolving networks from cross-correlation analysis of geomagnetic pulsation signals across high-latitude magnetometer stations. Nodes represent stations; edges represent statistically significant coherence with associated time lags. The spatiotemporal structure of these networks reveals how energy propagates through the magnetosphere-ionosphere system during extreme events.

### Pipeline

```
SuperMAG CSV data
    -> Butterworth bandpass filter (isolate Pc1-Pc5 bands)
    -> Windowed cross-correlation (station pairs)
    -> Threshold significant correlations
    -> Build NetworkX graphs (directed + undirected)
    -> Analyse: community detection, degree distributions, temporal snapshots
    -> Visualise on geographic maps (Cartopy)
```

## Setup

```bash
pip install numpy pandas scipy matplotlib networkx dynetx python-igraph cartopy PyAstronomy tqdm
```

## Usage

```bash
# Main network construction
python networkx_pcmodel.py

# Parallel Pc power calculation
python parallel_pc_power.py

# Network with surrogate testing
python surrogate_pc_net.py

# Read and analyse results
python read_networks.py

# Geographic visualisation with communities
python Drawnet_community.py
```

## Data

Download magnetometer time-series from [SuperMAG](https://supermag.jhuapl.edu/) as CSV. Station metadata is in `supermag-stations.csv`.

## Key Scripts

| Script | Purpose |
|--------|---------|
| `networkx_pcmodel.py` | Main network construction pipeline |
| `parallel_pc_power.py` | Parallel Pc power spectral analysis |
| `surrogate_pc_net.py` | Network construction with surrogate validation |
| `Drawnet_community.py` | Geographic network plot with community detection |
| `read_networks.py` | Network reading and property analysis |
| `ULF_bp.py` | ULF bandpass filtering |
| `NetDistr.py` | Network distribution analysis class |

## Author

Shahbaz Chaudhary ([@shahbaz22](https://github.com/shahbaz22))
