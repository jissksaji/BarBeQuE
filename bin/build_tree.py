#!/usr/bin/env python3
import sys
from ete3 import NCBITaxa

def main(taxid, output_file):
    ncbi = NCBITaxa()
    
    print(f"Fetching taxonomy tree for TaxID {taxid}...")
    try:
        # Get the full taxonomy tree including intermediate nodes
        tree = ncbi.get_descendant_taxa(taxid, collapse_subspecies=False, return_tree=True, intermediate_nodes=True)
    except Exception as e:
        print(f"Error finding descendants: {e}")
        return

    print("Annotating nodes with scientific names...")
    # Gather all node taxids to translate them in bulk for performance
    all_taxids = [int(node.name) for node in tree.traverse()]
    
    # ETE3 get_taxid_translator returns a dictionary: {taxid: ['Scientific Name']}
    name_translator = ncbi.get_taxid_translator(all_taxids)
    
    # Optional: get ranks as well
    ranks = ncbi.get_rank(all_taxids)

    # Replace the node IDs with scientific names for better readability in the output tree
    for node in tree.traverse():
        tid = int(node.name)
        sci_name = name_translator.get(tid, [str(tid)])[0]
        rank = ranks.get(tid, "no rank")
        
        # Replace spaces with underscores for Newick compatibility
        safe_name = sci_name.replace(" ", "_").replace("(", "").replace(")", "").replace(":", "")
        node.name = f"{safe_name}__{rank}"
        
    print(f"Saving the tree to {output_file}...")
    # format=1 saves internal node names (useful for families/orders)
    tree.write(outfile=output_file, format=1, format_root_node=True)
    print(f"Done! Tree saved with {len(all_taxids)} nodes.")
    print("\nTip: Since the tree contains hundreds of thousands of nodes (families, species, etc.),")
    print("it is recommended to open the .nwk file using dedicated tree viewers like:")
    print(" - iTOL (Interactive Tree Of Life - web based)")
    print(" - FigTree (Desktop application)")
    print(" - Dendroscope")

if __name__ == '__main__':
    taxid = 3398
    output_nwk = "taxonomy_tree.nwk"
    
    if len(sys.argv) > 1:
        try:
            taxid = int(sys.argv[1])
        except ValueError:
            print("Please provide a valid integer taxid.")
            sys.exit(1)
            
    if len(sys.argv) > 2:
        output_nwk = sys.argv[2]
            
    main(taxid, output_nwk)
