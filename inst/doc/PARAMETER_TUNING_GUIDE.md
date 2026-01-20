# OncoImplexus Parameter Tuning Guide

This guide provides recommendations for tuning detection parameters in OncoImplexus to optimize results for different sample types and data quality levels.

## Quick Reference

### Default Parameters for Standard WGS Data

```r
result <- detect_chromoanagenesis(
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    
    # Chromothripsis defaults
    min_cluster_size = 10,
    min_sv_calls = 10,
    
    # Chromoplexy defaults  
    min_chromosomes = 3,
    min_translocations = 3,
    max_cn_change = 1,
    
    # Chromoanasynthesis defaults
    cn_gradient_threshold = 0.5,
    min_sv_density = 3
)
```

---

## Chromothripsis Detection Parameters

### `min_cluster_size` (default: 10)

Controls the minimum number of clustered SVs required to consider a region as chromothripsis.

| Sample Type | Recommended Value | Notes |
|-------------|------------------|-------|
| High-quality WGS (>60x) | 10 | Standard threshold |
| Low-coverage WGS (30-40x) | 7-8 | Lower to compensate for missed calls |
| Panel sequencing | 5-7 | Limited coverage of breakpoints |

### `min_sv_calls` (default: 10)

Minimum total SV calls in sample. Samples below this threshold may not have enough data for reliable detection.

### `min_oscillating_cn` (default: 7)

Minimum oscillating copy number segments. Classic chromothripsis shows 2-state CN oscillations.

**Adjustment guidance:**
- Increase to 10+ for samples with known high CN noise
- Decrease to 5-6 for low-purity samples

---

## Chromoplexy Detection Parameters

### `min_chromosomes` (default: 3)

Minimum chromosomes involved in a translocation chain.

| Setting | Use Case |
|---------|----------|
| 2 | Capture simple balanced translocations |
| 3 | Standard chromoplexy (recommended) |
| 4+ | High-confidence complex events only |

### `min_translocations` (default: 3)

Minimum inter-chromosomal translocations required.

### `max_cn_change` (default: 1)

Maximum copy number change allowed between adjacent segments. Chromoplexy typically shows CN-neutral rearrangements.

- **1**: Strict CN neutrality (classic chromoplexy)
- **2**: Allow some CN variation (tumors with additional events)

### `adjacency_distance` (default: 10Mb)

Maximum genomic distance for considering breakpoints as adjacent.

### `max_path_search` (default: 10000)

Maximum paths to explore during chain detection. Increase for complex samples (>100 translocations).

---

## Chromoanasynthesis Detection Parameters

### `cn_gradient_threshold` (default: 0.5)

Threshold for detecting copy number gradients (typical of replication-based mechanisms).

- **Higher (0.7-1.0)**: More stringent, captures clear gradients only
- **Lower (0.3-0.5)**: More sensitive, may include false positives

### `min_sv_density` (default: 3)

Minimum SVs per Mb in a region.

### `template_switch_distance` (default: 2000 bp)

Maximum distance between breakpoints to consider as template switching (FoSTeS/MMBIR signature).

---

## Sample-Specific Recommendations

### High-Purity Tumor Samples (>80%)

```r
# Standard settings work well
result <- detect_chromoanagenesis(sv_data, cnv_data)
```

### Low-Purity Samples (<40%)

```r
result <- detect_chromoanagenesis(
    sv_data, cnv_data,
    min_cluster_size = 7,
    min_oscillating_cn = 5,
    max_cn_change = 2  # Allow more CN variation
)
```

### FFPE Samples

FFPE samples may have artifacts. Consider:
- Filtering SVs with low support reads
- Increasing `min_cluster_size` to 12-15
- Using stricter confidence thresholds

### Pediatric Tumors

Some pediatric cancers have simpler genomes:
```r
result <- detect_chromoanagenesis(
    sv_data, cnv_data,
    min_cluster_size = 8,
    min_chromosomes = 2
)
```

---

## Performance Tuning

### Large SV Datasets (>500 SVs)

```r
result <- detect_chromoanagenesis(
    sv_data, cnv_data,
    max_path_search = 5000,     # Limit search space
    max_neighbors = 30          # Reduce graph complexity
)
```

### Speed vs Sensitivity Trade-off

| Priority | Settings |
|----------|----------|
| Speed | `max_path_search = 1000`, `max_neighbors = 20` |
| Balanced | Default settings |
| Maximum sensitivity | `max_path_search = 20000`, `max_neighbors = 100` |

---

## Validation Recommendations

1. **Visual inspection**: Always visualize high-confidence events with `plot_chromothripsis()` or circos plots
2. **Cross-validation**: Compare results with different parameter sets
3. **Literature comparison**: Check if detected events match known patterns for tumor type

---

## Troubleshooting

### Too Many False Positives

- Increase `min_cluster_size`
- Increase statistical significance thresholds
- Check input data quality

### Missing Known Events

- Decrease `min_cluster_size`
- Lower `cn_gradient_threshold`
- Check if SVs are being filtered during input parsing

### Slow Performance

- Reduce `max_path_search`
- Filter low-confidence SVs before analysis
- Consider analyzing chromosomes in parallel

---

## Version Information

This guide is applicable to OncoImplexus v0.2.1 and later.

For questions or issues, please refer to the GitHub repository or contact the maintainers.
