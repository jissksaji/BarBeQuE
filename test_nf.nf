def taxdump_path = "/tmp"
def possible = file("${taxdump_path}/*{accession2taxid,genbank2taxid,nucl_}*")
println possible.getClass()
println possible
