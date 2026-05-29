from ete3 import NCBITaxa
ncbi = NCBITaxa()
taxid = ncbi.get_name_translator(["Hominidae"])["Hominidae"][0]
print("taxid", taxid)
tree = ncbi.get_descendant_taxa(taxid, collapse_subspecies=True, return_tree=True)
for n in tree.traverse():
    if n.rank == "species":
        print(n.sci_name)
