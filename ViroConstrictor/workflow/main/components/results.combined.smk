# Combined results component for human-friendly organization
# This component creates an additional results structure organized by sample and by virus
# rather than by Virus/RefID, making it easier to navigate when working with
# many references per sample.

def modify_fasta_headers(input_files, output_file, sample_virus_refid_map):
    """
    Combine FASTA files and modify headers to include Virus and RefID information.
    
    Parameters
    ----------
    input_files : list
        List of input FASTA file paths
    output_file : str
        Path to output FASTA file
    sample_virus_refid_map : dict
        Mapping of input file paths to (virus, refid) tuples
    """
    from Bio import SeqIO
    
    with open(output_file, 'w') as out_handle:
        for infile in input_files:
            if not os.path.exists(infile) or os.path.getsize(infile) == 0:
                continue
            
            virus, refid = sample_virus_refid_map[infile]
            
            for record in SeqIO.parse(infile, "fasta"):
                # Original header format: >sampleID mincov=X
                # New format: >sampleID Virus RefID mincov=X
                header_parts = record.description.split()
                if len(header_parts) >= 2 and 'mincov=' in header_parts[-1]:
                    sample_id = header_parts[0]
                    mincov_part = header_parts[-1]
                    new_header = f"{sample_id} {virus} {refid} {mincov_part}"
                    record.id = sample_id
                    record.description = new_header
                else:
                    # Fallback: just append Virus and RefID
                    new_header = f"{record.description} {virus} {refid}"
                    record.description = new_header
                
                SeqIO.write(record, out_handle, "fasta")


# ===========================
# BY SAMPLE AGGREGATION RULES
# ===========================

# Rule to combine consensus sequences by sample with modified headers
rule combine_consensus_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{cons}{seqs}{wc.sample}.fa"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/consensus.fasta"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    run:
        # Create mapping of file to (virus, refid)
        sample_data = samples_df[samples_df["sample"] == wildcards.sample]
        file_map = {}
        for _, row in sample_data.iterrows():
            filepath = f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{cons}{seqs}{wildcards.sample}.fa"
            file_map[filepath] = (row['Virus'], row['RefID'])
        
        modify_fasta_headers(input, output[0], file_map)


# Rule to combine mutations by sample (adds Virus column)
rule combine_mutations_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{aln}{vf}{wc.sample}.tsv"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/mutations.tsv"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        # Create mapping of files to virus names
        sample_data = samples_df[samples_df["sample"] == wildcards.sample]
        file_to_virus = {}
        for _, row in sample_data.iterrows():
            filepath = f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{aln}{vf}{wildcards.sample}.tsv"
            file_to_virus[filepath] = row['Virus']
        
        # Read all input files and combine them
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                df = pd.read_csv(infile, sep='\t')
                # Add Virus column
                df['Virus'] = file_to_virus.get(infile, 'Unknown')
                dfs.append(df)
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Reorder columns to put Virus after Sample
            cols = combined_df.columns.tolist()
            if 'Sample' in cols and 'Virus' in cols:
                cols.remove('Virus')
                sample_idx = cols.index('Sample')
                cols.insert(sample_idx + 1, 'Virus')
                combined_df = combined_df[cols]
            combined_df.to_csv(output[0], sep='\t', index=False)
        else:
            # Create empty file with header
            with open(output[0], 'w') as f:
                f.write('Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n')


# Rule to combine width of coverage by sample (using correct BoC files, adds Virus and RefID columns)
rule combine_coverage_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{boc}{wc.sample}.tsv"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/Width_of_coverage.tsv"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        # Create mapping of files to virus and refid
        sample_data = samples_df[samples_df["sample"] == wildcards.sample]
        file_to_meta = {}
        for _, row in sample_data.iterrows():
            filepath = f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{boc}{wildcards.sample}.tsv"
            file_to_meta[filepath] = (row['Virus'], row['RefID'])
        
        # BoC files have NO HEADER - just tab-separated data
        # Column order: Sample_name, Width_at_mincov_1, Width_at_mincov_5, Width_at_mincov_10, Width_at_mincov_50, Width_at_mincov_100
        column_names = ['Sample_name', 'Width_at_mincov_1', 'Width_at_mincov_5', 'Width_at_mincov_10', 'Width_at_mincov_50', 'Width_at_mincov_100']
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                # Read without header and assign column names
                df = pd.read_csv(infile, sep='\t', header=None, names=column_names)
                # Add Virus and RefID columns with actual values
                virus, refid = file_to_meta.get(infile, ('Unknown', 'Unknown'))
                df['Virus'] = virus
                df['Reference_ID'] = refid
                dfs.append(df)
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Reorder columns: Sample_name, Virus, Reference_ID, then coverage columns
            combined_df = combined_df[['Sample_name', 'Virus', 'Reference_ID', 'Width_at_mincov_1', 'Width_at_mincov_5', 'Width_at_mincov_10', 'Width_at_mincov_50', 'Width_at_mincov_100']]
            combined_df.to_csv(output[0], sep='\t', index=False)
        else:
            # Create empty file with header
            with open(output[0], 'w') as f:
                f.write('Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n')


# Rule to combine amplicon coverage by sample (adds Virus and RefID columns)
rule combine_amplicon_coverage_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{prim}{wc.sample}_ampliconcoverage.csv"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/Amplicon_coverage.csv"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        # Create mapping of files to virus and refid
        sample_data = samples_df[samples_df["sample"] == wildcards.sample]
        file_to_meta = {}
        for _, row in sample_data.iterrows():
            filepath = f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{prim}{wildcards.sample}_ampliconcoverage.csv"
            file_to_meta[filepath] = (row['Virus'], row['RefID'])
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                df = pd.read_csv(infile)
                # Add Virus and RefID columns
                virus, refid = file_to_meta.get(infile, ('Unknown', 'Unknown'))
                df['Virus'] = virus
                df['Reference_ID'] = refid
                dfs.append(df)
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Move Virus and Reference_ID to front
            cols = combined_df.columns.tolist()
            if 'Virus' in cols and 'Reference_ID' in cols:
                cols.remove('Virus')
                cols.remove('Reference_ID')
                cols = ['Virus', 'Reference_ID'] + cols
                combined_df = combined_df[cols]
            combined_df.to_csv(output[0], index=False)
        else:
            # Create empty CSV
            pd.DataFrame().to_csv(output[0], index=False)


# ===========================
# BY VIRUS AGGREGATION RULES
# ===========================

# Rule to combine consensus sequences by virus with modified headers
rule combine_consensus_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{cons}{seqs}{row['sample']}.fa"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}consensus.fasta"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    run:
        # Create mapping of file to (virus, refid)
        virus_data = samples_df[samples_df["Virus"] == wildcards.Virus]
        file_map = {}
        for _, row in virus_data.iterrows():
            filepath = f"{datadir}Virus~{wildcards.Virus}/RefID~{row['RefID']}/{cons}{seqs}{row['sample']}.fa"
            file_map[filepath] = (row['Virus'], row['RefID'])
        
        modify_fasta_headers(input, output[0], file_map)


# Rule to combine mutations by virus
rule combine_mutations_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{aln}{vf}{row['sample']}.tsv"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}mutations.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        # Create mapping of files to virus names
        virus_data = samples_df[samples_df["Virus"] == wildcards.Virus]
        file_to_virus = {}
        for _, row in virus_data.iterrows():
            filepath = f"{datadir}Virus~{wildcards.Virus}/RefID~{row['RefID']}/{aln}{vf}{row['sample']}.tsv"
            file_to_virus[filepath] = row['Virus']
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    df = pd.read_csv(infile, sep='\t')
                    if not df.empty:
                        # Add Virus column
                        df['Virus'] = file_to_virus.get(infile, wildcards.Virus)
                        dfs.append(df)
                except pd.errors.EmptyDataError:
                    continue
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Reorder columns to put Virus after Sample
            cols = combined_df.columns.tolist()
            if 'Sample' in cols and 'Virus' in cols:
                cols.remove('Virus')
                sample_idx = cols.index('Sample')
                cols.insert(sample_idx + 1, 'Virus')
                combined_df = combined_df[cols]
            combined_df.to_csv(output[0], sep='\t', index=False)
        else:
            with open(output[0], 'w') as f:
                f.write('Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n')


# Rule to combine coverage by virus
rule combine_coverage_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{boc}{row['sample']}.tsv"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}Width_of_coverage.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        # Create mapping of files to virus and refid
        virus_data = samples_df[samples_df["Virus"] == wildcards.Virus]
        file_to_meta = {}
        for _, row in virus_data.iterrows():
            filepath = f"{datadir}Virus~{wildcards.Virus}/RefID~{row['RefID']}/{boc}{row['sample']}.tsv"
            file_to_meta[filepath] = (row['Virus'], row['RefID'])
        
        # BoC files have NO HEADER - just tab-separated data
        column_names = ['Sample_name', 'Width_at_mincov_1', 'Width_at_mincov_5', 'Width_at_mincov_10', 'Width_at_mincov_50', 'Width_at_mincov_100']
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    # Read without header and assign column names
                    df = pd.read_csv(infile, sep='\t', header=None, names=column_names)
                    if not df.empty:
                        # Add Virus and RefID columns with actual values
                        virus, refid = file_to_meta.get(infile, (wildcards.Virus, 'Unknown'))
                        df['Virus'] = virus
                        df['Reference_ID'] = refid
                        dfs.append(df)
                except pd.errors.EmptyDataError:
                    continue
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Reorder columns: Sample_name, Virus, Reference_ID, then coverage columns
            combined_df = combined_df[['Sample_name', 'Virus', 'Reference_ID', 'Width_at_mincov_1', 'Width_at_mincov_5', 'Width_at_mincov_10', 'Width_at_mincov_50', 'Width_at_mincov_100']]
            combined_df.to_csv(output[0], sep='\t', index=False)
        else:
            with open(output[0], 'w') as f:
                f.write('Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n')


# Rule to combine amplicon coverage by virus
rule combine_amplicon_coverage_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{prim}{row['sample']}_ampliconcoverage.csv"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}Amplicon_coverage.csv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        # Create mapping of files to virus and refid
        virus_data = samples_df[samples_df["Virus"] == wildcards.Virus]
        file_to_meta = {}
        for _, row in virus_data.iterrows():
            filepath = f"{datadir}Virus~{wildcards.Virus}/RefID~{row['RefID']}/{prim}{row['sample']}_ampliconcoverage.csv"
            file_to_meta[filepath] = (row['Virus'], row['RefID'])
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    df = pd.read_csv(infile)
                    if not df.empty:
                        # Add Virus and RefID columns
                        virus, refid = file_to_meta.get(infile, (wildcards.Virus, 'Unknown'))
                        df['Virus'] = virus
                        df['Reference_ID'] = refid
                        dfs.append(df)
                except pd.errors.EmptyDataError:
                    continue
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Move Virus and Reference_ID to front
            cols = combined_df.columns.tolist()
            if 'Virus' in cols and 'Reference_ID' in cols:
                cols.remove('Virus')
                cols.remove('Reference_ID')
                cols = ['Virus', 'Reference_ID'] + cols
                combined_df = combined_df[cols]
            combined_df.to_csv(output[0], index=False)
        else:
            pd.DataFrame().to_csv(output[0], index=False)


# =============================
# ALL SAMPLES AGGREGATION RULES
# =============================

# Global aggregation rule for consensus (from by_sample results)
rule combine_all_consensus:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/consensus.fasta",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_consensus.fasta"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    shell:
        """
        cat {input} > {output}
        """


# Global aggregation rule for mutations (from by_sample results)
rule combine_all_mutations:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/mutations.tsv",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_mutations.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    df = pd.read_csv(infile, sep='\t')
                    if not df.empty:
                        dfs.append(df)
                except pd.errors.EmptyDataError:
                    continue
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv(output[0], sep='\t', index=False)
        else:
            with open(output[0], 'w') as f:
                f.write('Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n')


# Global aggregation rule for coverage (from by_sample results)
rule combine_all_coverage:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/Width_of_coverage.tsv",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_width_of_coverage.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    # These files now have headers (created by combine_coverage_by_sample)
                    df = pd.read_csv(infile, sep='\t')
                    if not df.empty:
                        dfs.append(df)
                except pd.errors.EmptyDataError:
                    continue
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv(output[0], sep='\t', index=False)
        else:
            with open(output[0], 'w') as f:
                f.write('Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n')


# Global aggregation rule for amplicon coverage (from by_sample results)
rule combine_all_amplicon_coverage:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/Amplicon_coverage.csv",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_amplicon_coverage.csv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    run:
        import pandas as pd
        
        dfs = []
        for infile in input:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    df = pd.read_csv(infile)
                    if not df.empty:
                        dfs.append(df)
                except pd.errors.EmptyDataError:
                    continue
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv(output[0], index=False)
        else:
            pd.DataFrame().to_csv(output[0], index=False)