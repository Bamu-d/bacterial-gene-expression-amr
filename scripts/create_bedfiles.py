import os
import pandas as pd
import argparse

def create_bed_files(input_file, output_dir, chromosome):
    # Read the gene positions file
    genes = pd.read_csv(input_file, sep='\t', index_col =1)
    
    # Process forward strand genes
    forward = genes[genes['strand'] == '+'].copy()
    forward['begin'] = forward.start.sub(200)
    forward['stop'] = forward.start.add(30)
    
    # Process reverse strand genes
    reverse = genes[genes['strand'] == '-'].copy()
    reverse['begin'] = reverse.end.sub(30)
    reverse['stop'] = reverse.end.add(200)
    
    # Combine forward and reverse
    df3 = pd.concat([forward, reverse], axis=0)
    df3 = df3.drop(['strand', 'start', 'end'], axis=1)
    
    # Add chromosome information
    genes['chro'] = chromosome
    genes = genes[['chro', 'start', 'end', 'strand', 'gene_name', 'locus_tag']]
    
    # Create the final bed dataframe
    bed = genes.set_index('locus_tag').join(df3.drop(columns='gene_name').set_index('locus_tag')).reset_index()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write individual bed files
    for _, row in bed.iterrows():
        name = row['gene_name']
        locus = row['locus_tag']
        chrom = row['chro']
        start = row['begin']
        end = row['stop']
        
        with open(f'{output_dir}/{locus}.bed', 'w') as o_bed:
            o_bed.write(f'{chrom}\t{start}\t{end}\t{locus}\t{name}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create BED files from gene positions")
    parser.add_argument("input_file", help="Path to the input gene positions file")
    parser.add_argument("output_dir", help="Path to the output directory for BED files")
    parser.add_argument("--chromosome", default="NC_000913", help="Chromosome identifier (default: NC_000913)")
    args = parser.parse_args()

    create_bed_files(args.input_file, args.output_dir, args.chromosome)