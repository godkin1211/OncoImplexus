# OncoImplexus Algorithm Details

**OncoImplexus** is an integrated framework for the detection and classification of complex chromosomal rearrangements (chromoanagenesis) in cancer genomes. This document details the mathematical models and algorithmic logic for its three primary modules: **Chromothripsis**, **Chromoplexy**, and **Chromoanasynthesis**.

## 1. Chromothripsis Detection

The chromothripsis detection module is adapted from the statistical framework established by Cortes-Ciriano et al. (Nature Genetics, 2020) in the *ShatterSeek* algorithm. It identifies clusters of interleaved structural variants (SVs) and evaluates them against statistical criteria defining random shattering and reassembly.

### 1.1 SV Clustering
SVs are first clustered based on genomic proximity. Let $SV_i$ and $SV_j$ be two intrachromosomal SVs with breakpoints $B_i = \{b_{i1}, b_{i2}\}$ and $B_j = \{b_{j1}, b_{j2}\}$.
A cluster $C$ is defined as a connected component in a graph where nodes are SVs and an edge exists between $SV_i$ and $SV_j$ if their spanned intervals overlap or are contained within one another.

### 1.2 Statistical Criteria

For each identified cluster, the following statistical tests are performed:

#### A. Clustering of Breakpoints (Exponential Distribution Test)
To test if breakpoints are randomly distributed (Poisson process) versus clustered, we analyze the inter-breakpoint distances.
Let $X = \{x_1, x_2, ..., x_n\}$ be the sorted genomic positions of all breakpoints in a cluster. The inter-breakpoint distances are $d_i = x_{i+1} - x_i$.
Under the null hypothesis of random breakage, the distances $d_i$ should follow an Exponential distribution:
$$ f(x; \lambda) = \lambda e^{-\lambda x} $$ 
where $\lambda$ is the rate parameter, estimated as $\hat{\lambda} = 1/\bar{d}$.
A **Kolmogorov-Smirnov (KS) test** is performed to compare the empirical distribution of $d_i$ against the theoretical exponential distribution. A low p-value allows rejection of the random breakage hypothesis, supporting clustering.

*Note: OncoImplexus (v1.0.0+) implements a critical fix in this test by correctly extracting and ordering breakpoint coordinates, resolving a coordinate-swap bug found in earlier frameworks that led to significant false-negative rates in complex cancer genomes.*

#### B. Randomness of Fragment Joins (Bernoulli/Multinomial Test)
Chromothripsis involves the random rejoining of DNA fragments. This implies that the four types of fragment joins (resulting in different SV types) should occur with equal probability.
The four join types (SV types) are:
1.  **Deletion-like (DEL):** (+, -) orientation
2.  **Duplication-like (DUP):** (-, +) orientation
3.  **Head-to-Head Inversion (h2hINV):** (+, +) orientation
4.  **Tail-to-Tail Inversion (t2tINV):** (-, -) orientation

We use a **Chi-squared Goodness-of-Fit test** (or exact Multinomial test) to test the null hypothesis:
$$ P(\text{DEL}) = P(\text{DUP}) = P(\text{h2hINV}) = P(\text{t2tINV}) = 0.25 $$ 
Significant deviation (low p-value) suggests non-random mechanisms (e.g., selection pressure), though "High Confidence" chromothripsis typically requires adhering to this random profile.

#### C. Chromosomal Enrichment
We test if the chromosome containing the cluster has significantly more SVs than expected by chance, given its length relative to the genome.
Let $N_{chr}$ be the number of SVs on the candidate chromosome, and $N_{total}$ be the total SVs in the sample. Let $L_{chr}$ and $L_{genome}$ be the effective lengths.
The expected probability is $p = L_{chr} / L_{genome}$.
A **Binomial Test** calculates the probability of observing $\ge N_{chr}$ events:
$$ P(X \ge N_{chr}) = \sum_{k=N_{chr}}^{N_{total}} \binom{N_{total}}{k} p^k (1-p)^{N_{total}-k} $$ 

#### D. Copy Number Oscillations
Chromothripsis typically produces oscillating copy number (CN) states (e.g., between 1 and 2 copies).
We calculate the maximum number of consecutive oscillating CN segments. Let $C = (c_1, c_2, ..., c_m)$ be the sequence of total copy numbers for adjacent segments in the region.
An oscillation is defined as a transition $c_i \to c_{i+1}$ where $|c_i - c_{i+1}| \le 1$.
The metric $M_{osc}$ counts the length of the longest sub-sequence satisfying oscillating conditions (2-state or 3-state).

---\n

## 2. Chromoplexy Detection (v3.0)

OncoImplexus introduces a novel **Graph-Based Chaining Algorithm** (v3.0) to detect chromoplexy, overcoming limitations in previous methods that failed to connect chains when local rearrangements or deletions occurred at breakpoints.

### 2.1 Graph Construction
We construct a graph $G = (V, E)$ where nodes $V$ represent genomic breakpoints.
The edge set $E$ is the union of three distinct edge types: $E = E_{tra} \cup E_{adj} \cup E_{del}$.

#### A. Translocation Edges ($E_{tra}$)
Represent observed inter-chromosomal structural variants.
For a translocation between $chrA:pos1$ and $chrB:pos2$:
$$ (v_{A,1}, v_{B,2}) \in E_{tra} $$ 
Weight: $w_{tra} = 1.0$

#### B. Genomic Adjacency Edges ($E_{adj}$)
Represent the physical proximity of breakpoints on the *same* chromosome, allowing the algorithm to "walk" along the chromosome to find the next translocation.
For two breakpoints $v_i, v_j$ on the same chromosome with distance $d_{ij} = |pos_i - pos_j|$:
$$ (v_i, v_j) \in E_{adj} \iff d_{ij} < D_{max\_adj} $$ 
where $D_{max\_adj}$ is the adjacency threshold (default 10Mb).
The edge weight decays exponentially with distance:
$$ w_{adj}(d_{ij}) = e^{-d_{ij} / \sigma} $$ 
where $\sigma$ is a decay scale (default 5Mb).

#### C. Deletion Bridges ($E_{del}$)
Represent inferred deletions that connect two breakpoints (e.g., a "bridge" formed by the loss of the segment between them).
A deletion bridge is identified if a CNV deletion segment $S_{del}$ spans or is adjacent to breakpoints $v_i, v_j$.
$$ (v_i, v_j) \in E_{del} $$ 
Confidence score $C_{bridge}$ is calculated based on distance, deletion size, and copy number:
$$ C_{bridge} = 0.5 \cdot e^{-d/200kb} + 0.3 \cdot S_{size} + 0.2 \cdot e^{-(CN)/1} $$ 

### 2.2 Chain Discovery (Backtracking)
We search for chains (paths) in $G$ that satisfy chromoplexy criteria using a backtracking algorithm (Depth-First Search with constraints).
A path $P = (v_1, v_2, ..., v_k)$ is a valid chromoplexy chain candidate if:
1.  **Length:** $N_{tra}(P) \ge N_{min\_tra}$ (Minimum translocations, default 3).
2.  **Complexity:** $N_{chr}(P) \ge N_{min\_chr}$ (Minimum chromosomes involved, default 3).
3.  **Graph Structure:** The path traverses a sequence of translocation edges interleaved with adjacency/deletion edges.

The search space is pruned using:
*   `max_depth`: Maximum chain length.
*   `max_neighbors`: Limit on adjacency edges per node (optimization for complex graphs).

### 2.3 Evaluation and Scoring
Each candidate chain is scored:
$$ Score_{combined} = 0.3 \cdot S_{CN} + 0.3 \cdot S_{complexity} + 0.2 \cdot F_{del} + 0.2 \cdot S_{stat} $$ 
Where:
*   $S_{CN} = e^{-\text{mean}(|CN - 2|)/2}$: Copy number stability score.
*   $S_{complexity}$: Normalized score based on number of chromosomes and translocations.
*   $F_{del}$: Fraction of edges that are deletion bridges.

---\n

## 3. Chromoanasynthesis Detection

This module detects replication-based mechanisms (Fork Stalling and Template Switching - FoSTeS / MMBIR), characterized by local clustering of SVs, copy number gradients, and linked breakpoint chains.

### 3.1 Copy Number Gradient Detection
Chromoanasynthesis often produces a stepwise increase in copy number. We detect this using a sliding window approach.
For a window of CN segments $W$, we calculate:
1.  **Spearman Correlation ($\rho$):** Correlation between segment index and Total CN.
    $$ \rho = \text{cor}(\text{rank}(index), \text{rank}(CN)) $$ 
2.  **Total Variation ($V$):** Measures local instability ("jumpiness").
    $$ V = \frac{\sum |CN_{i+1} - CN_i|}{\text{mean}(CN)} $$ 

A region is flagged if $\rho > \rho_{thresh}$ (default 0.4) or if $V$ is sufficiently high with a positive trend.

### 3.2 SV Topology Analysis
We analyze the topology of SVs within the flagged regions to distinguish FoSTeS from Breakage-Fusion-Bridge (BFB) cycles.

#### A. Linked SV Chains (FoSTeS Signature)
FoSTeS involves serial template switching. We detect this by finding chains of "linked" SVs.
Two SVs, $SV_i$ and $SV_j$, are **linked** if any of their breakpoints are within a small distance $\delta$ (default 2kb, representing the replication fork jumping distance).
$$ Linked(SV_i, SV_j) \iff \min_{b \in B_i, b' \in B_j} |b - b'| < \delta $$ 
We construct an adjacency graph of SVs and find connected components (chains).
$$ L_{chain} = \text{max}(\text{size}(\text{ConnectedComponents})) $$ 
High $L_{chain}$ ($\ge 3$ or $4$) is a strong indicator of chromoanasynthesis.

#### B. Fold-back Inversions (BFB Signature)
BFB cycles are characterized by "fold-back" inversions (very short span inversions).
$$ SV \in \text{FoldBack} \iff SV_{type} \in \{h2hINV, t2tINV\} \land |pos2 - pos1| < 20kb $$ 
High density of fold-backs penalizes the chromoanasynthesis score (suggesting BFB instead).

### 3.3 Classification
Final classification relies on a weighted complexity score:
$$ Score = 0.2 \cdot S_{grad} + 0.2 \cdot S_{var} + 0.2 \cdot S_{dup} + 0.4 \cdot S_{chain} $$ 
Events are classified as "Likely" or "Possible" based on this score and the absence of BFB features.

---\n

## 4. DNA Repair Mechanism Inference

OncoImplexus analyzes the base-level sequence features at structural variant junctions to infer the operative DNA double-strand break (DSB) repair pathways.

### 4.1 Sequence Extraction and Orientation
To detect microhomology accurately, reference sequences flanking each breakpoint are extracted and re-oriented based on the SV strand. 

Let $S_1$ and $S_2$ be the reference sequences at the two breakpoints. They are transformed into:
*   **Kept Sequence ($s_{kept}$):** The sequence retained in the rearranged genome, oriented facing the junction.
*   **Discarded Sequence ($s_{disc}$):** The sequence "lost" or "overlapped" during the rearrangement, starting exactly from the breakpoint position.

This normalization ensures that regardless of the original SV type (DEL, DUP, INV), the junction interface always follows a unified coordinate system where the junction-adjacent base of Side 1 is at the start of its discarded fragment, and Side 2 is at the end.

### 4.2 Microhomology Detection (Unified Suffix-Prefix Algorithm)
Microhomology is defined as an identical sequence present at both breakpoint interfaces that is represented only once in the resulting junction.

We identify the maximum length $k$ such that:
$$ \text{prefix}(s_{1,disc}, k) == \text{suffix}(s_{2,disc}, k) $$ 
where $k \le \text{max\_microhomology}$ (default 25 bp). 

The match must be **continuous** starting from the junction. Any mismatch terminates the search, ensuring that only biologically relevant microhomology is reported.

### 4.3 Classification Criteria
The inferred repair mechanism is determined by the length of microhomology ($L_{MH}$) and the type of SV:

| Mechanism | Criteria | Biological Significance |
|:---|:---|:---|
| **NHEJ** | $L_{MH} \in \{0, 1\}$ bp | Classic Non-Homologous End Joining; often associated with Chromothripsis. |
| **MMEJ** | $2 \le L_{MH} \le 20$ bp | Microhomology-Mediated End Joining (or Alt-EJ); suggests error-prone repair. |
| **SSA** | $L_{MH} > 20$ bp | Single-Strand Annealing; involves long homologous repeats. |
| **MMBIR** | $2 \le L_{MH} \le 5$ bp AND SV = DUP | Microhomology-Mediated Break-Induced Replication; a hallmark of Chromoanasynthesis. |

### 4.4 Insertion Detection
If no microhomology is found, the algorithm checks for non-templated insertions ($L_{ins}$) using the VCF `INSSEQ` or `ALT` fields. The presence of small insertions (typically <10 bp) further supports an NHEJ repair origin.