# Workflow Visualization Guide

This guide explains how to generate and use visual documentation for the ViroConstrictor Snakemake workflow.

## Quick Start

### Using the DAG Generation Script

```bash
# Generate all DAG visualizations with default settings
./generate_dag.sh

# Use a specific config file
./generate_dag.sh --config path/to/your/config.yaml

# Specify custom output directory
./generate_dag.sh --output docs/images/
```

### Manual DAG Generation

If you prefer manual control, you can generate DAGs directly with Snakemake:

```bash
# Detailed workflow DAG (shows all jobs)
snakemake --snakefile ViroConstrictor/workflow/main/workflow.smk \
    --configfile your_config.yaml \
    --dag | dot -Tpng -Gdpi=300 > workflow_dag.png

# Rule graph (shows rule relationships)
snakemake --snakefile ViroConstrictor/workflow/main/workflow.smk \
    --configfile your_config.yaml \
    --rulegraph | dot -Tpng -Gdpi=300 > rules_dag.png

# File dependency graph
snakemake --snakefile ViroConstrictor/workflow/main/workflow.smk \
    --configfile your_config.yaml \
    --filegraph | dot -Tpng -Gdpi=300 > files_dag.png
```

## Types of Visualizations

### 1. Detailed Workflow DAG
- **Purpose**: Shows every job that will be executed
- **Best for**: Understanding the complete execution plan
- **File size**: Large (can be very detailed)
- **Use case**: Debugging, comprehensive workflow understanding

### 2. Rule Graph
- **Purpose**: Shows relationships between rules (simplified)
- **Best for**: Understanding workflow logic and dependencies
- **File size**: Medium
- **Use case**: Documentation, presentations, architecture overview

### 3. File Dependency Graph
- **Purpose**: Shows file inputs and outputs
- **Best for**: Understanding data flow
- **File size**: Medium to large
- **Use case**: Data lineage, troubleshooting file dependencies

## Prerequisites

### Required Software

```bash
# Install Graphviz (for DOT visualization)
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install graphviz

# macOS
brew install graphviz

# Conda
conda install graphviz

# Verify installation
dot -V
```

### Snakemake Environment
Ensure you have Snakemake installed and accessible:

```bash
# Check Snakemake version
snakemake --version

# Should be >= 7.0 for best compatibility
```

## Customizing Visualizations

### Enhanced SVG Output
For higher quality, scalable images:

```bash
snakemake --dag | dot -Tsvg -Gsize="12,10!" -Gdpi=300 > workflow.svg
```

### Styled DAG with Colors
Create a custom DOT style file (`dag_style.dot`):

```dot
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname="Arial"];
    edge[fontname="Arial"];
    
    // Color scheme for different rule types
    node[fillcolor=lightblue, style="rounded,filled"] preparation;
    node[fillcolor=lightgreen, style="rounded,filled"] cleaning;
    node[fillcolor=lightyellow, style="rounded,filled"] analysis;
    node[fillcolor=lightcoral, style="rounded,filled"] results;
}
```

### Large Workflow Optimization
For very large workflows, filter to specific parts:

```bash
# Generate DAG for specific target
snakemake --snakefile ViroConstrictor/workflow/main/workflow.smk \
    --configfile your_config.yaml \
    --dag results/sample1/consensus.fasta | dot -Tpng > sample1_dag.png
```

## Integration with Documentation

### Embedding in Markdown
```markdown
![Workflow Overview](images/workflow_rules.png)

*Figure 1: ViroConstrictor workflow rule dependencies*
```

### Creating Interactive Documentation
Use tools like MkDocs with the DAG images:

```yaml
# mkdocs.yml
nav:
  - Architecture: architecture.md
  - Workflow: workflow.md

markdown_extensions:
  - attr_list  # For image sizing
```

### Automated Documentation Updates
Add DAG generation to your CI/CD pipeline:

```yaml
# .github/workflows/docs.yml
- name: Generate DAG images
  run: |
    ./generate_dag.sh --output docs/images/
    git add docs/images/
    git commit -m "Update workflow DAGs" || exit 0
```

## Troubleshooting

### Common Issues

#### "dot command not found"
```bash
# Install Graphviz
sudo apt-get install graphviz  # Linux
brew install graphviz          # macOS
```

#### "Config file errors"
The DAG generation script creates a minimal config file automatically. If you're using your own config, ensure all required fields are present.

#### "Memory errors with large DAGs"
For very large workflows:
```bash
# Increase DOT memory limit
dot -Tpng -Gmaxiter=2000 -Gnslimit=2 < workflow.dot > workflow.png
```

#### "SVG files too large"
Use PNG for large workflows:
```bash
snakemake --rulegraph | dot -Tpng -Gdpi=150 > workflow.png
```

### Performance Tips

1. **Use rule graphs** for documentation (simpler than detailed DAGs)
2. **Generate SVGs** for scalable images in documentation
3. **Filter large workflows** to specific samples or rules
4. **Cache generated images** to avoid regenerating unchanged workflows

## Best Practices

### Documentation Standards

1. **Always include a rule graph** in your architecture documentation
2. **Use descriptive filenames** (e.g., `viroconstrictor_main_workflow.png`)
3. **Update visualizations** when workflow changes significantly
4. **Include legends** explaining color coding or symbols used

### Version Control

1. **Don't commit large DAG images** to git (use `.gitignore`)
2. **Generate images in CI/CD** for documentation builds
3. **Tag major workflow changes** with corresponding DAG updates

### Maintenance

1. **Regenerate DAGs** after major workflow updates
2. **Test DAG generation** as part of your testing pipeline
3. **Document any custom styling** or generation parameters used

This visualization approach will help users understand your pipeline architecture and assist in development, debugging, and documentation maintenance.