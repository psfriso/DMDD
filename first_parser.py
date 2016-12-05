import sys
import re

_test_lim = 1E10

mp_pattern = re.compile(r'^MP:[0-9]*')

class novel_gene( object ):

    def __init__(self, name, marker_id):
        self.name = name
        self.marker_id = marker_id
        self.alleles_mgi = set( )
        self.alleles_impc = set( )
        self.phenotypes = set( )
        self.new_phenotypes = set( )

    def find_alleles_in_mgi(self, l ):
        if l:
            for elem in l:
                self.alleles_mgi.add( elem )

    def find_alleles_in_impc(self, l ):
        if l :
            for elem in l:
                self.alleles_impc.add( elem )

    def collect_phenotypes( self, l ):
        if l:
            for MP_string in l:
                str_array = MP_string.split(",")
                for mp in str_array:
                    if is_mp_term(mp): self.phenotypes.add( mp )

    def new_alleles( self ):
        # alleles in impc that match (or are absent) in mgi
        not_even_in_mgi = self.alleles_impc - self.alleles_mgi

        intersection = self.alleles_impc & self.alleles_mgi

        return not_even_in_mgi | intersection # union of two sets

    def are_all_allele_new( self ):
        new = self.new_alleles( )
        old = self.known_alleles( )
        ans = ( len( self.alleles_mgi - new ) == 0 and len(new) >0 and len(old) == 0 )
        return ans

    def known_alleles( self ):
        return  self.alleles_mgi - self.alleles_impc

    def has_known_alleles( self ):
        return len( self.known_alleles() ) > 0

    def novel_phenotypes( self ):
        return self.new_phenotypes - self.phenotypes

def is_mp_term( elem ):
    return re.match(mp_pattern, elem )


def count_mp(l):
    terms = set()
    if l:
        for MP_string in l:
            str_array = MP_string.split(",")
            for mp in str_array:
                if is_mp_term(mp): terms.add( mp )
    return len( terms )



mgi = { }
with open(sys.argv[1] , 'r' ) as f:
    # MGI phenotypic allele, all avaiable data from downloaded files.
    c = 0
    for l in f:
        if c > _test_lim: break
        line_elem = l.rstrip().split("\t")
        c += 1
#        print( line_elem )
        try:
            mp_terms = line_elem[10]
        except IndexError:
            mp_terms = ""
        if line_elem[7] in mgi:
            mgi[ line_elem[7] ][0].append( line_elem[0] )
            mgi[ line_elem[7] ][1].append( line_elem[6] )
            mgi[ line_elem[7] ][2].append( mp_terms )
        else:
            mgi[ line_elem[7] ] = [ [line_elem[0]], [line_elem[6]] , [mp_terms] ]
#        print( line_elem[0] , line_elem[6], line_elem[7] , mp_terms )
f.close()


impc = {}
with open(sys.argv[2] , 'r' ) as f:
    # IMPC data, all of it ( version 5.0 beta)
    c = 0
    f.readline()
    for l in f:
        if c > _test_lim: break
        line_elem = l.rstrip().split("\t")
        c += 1
#        print( line_elem )
        mps = [ elem for elem in line_elem if is_mp_term(elem) ]
        try:
            if line_elem[1] in impc:
                impc[ line_elem[1] ][0].append( line_elem[6] )
                impc[ line_elem[1] ][1].append( line_elem[0] )
                impc[ line_elem[1] ][2].append(  ",".join(mps) )
            else:
                impc[ line_elem[1] ] = [  [line_elem[6]] , [line_elem[0]] , [ ",".join(mps) ] ]
        except IndexError:
            print( "\tSkipping\t",line_elem )

f.close()


gene_list = dict()
with open(sys.argv[3] , 'r' ) as f:
#   list of DMDD genes and maker id associated to them (one to one correspondence)
    for l in f:
        ID, marker = l.rstrip().split("\t")
        if ID in gene_list:
            if not (gene_list[ID].name == ID and gene_list[ID].marker_id == marker):
            # check
                print( "See ", gene_list[ID] )
        else:
            gene_list[ ID ] = novel_gene( ID, marker)

f.close()

for gene in gene_list:
#    print( gene_list[gene].name , gene_list[gene].marker_id )
    if gene in mgi:
        gene_list[ gene ].find_alleles_in_mgi( mgi[ gene ][0] )
        gene_list[ gene ].collect_phenotypes( mgi[ gene ][2] )
    else:
        gene_list[ gene ].find_alleles_in_mgi( None )

    if gene in impc:
        gene_list[ gene ].find_alleles_in_impc( impc[ gene ][0] )
        gene_list[ gene ].collect_phenotypes( impc[ gene ][2] )
    else:
        gene_list[ gene ].find_alleles_in_impc( None )

with open( sys.argv[4] , 'r' ) as f:
    # file 4 is DMDD's supplementary table with the new MP terms per gene
    f.readline( ) # header
    for l in f:
        aux = l.rstrip().split(",")
        gene = aux[0]
        marker = aux[1]
        mpterm = aux[4]
        try:
            assert gene_list[ gene ].marker_id == marker
            gene_list[ gene ].new_phenotypes.add( mpterm )
        except IndexError:
            print( gene, " not present in dictionary")

f.close()


# for gene in gene_list:
# #    print( gene_list[gene].name )
#     if gene_list[ gene ].are_all_allele_new():
#         print( gene_list[gene].name , gene_list[gene].marker_id , end = '\t')
#         print( "new", [ al for al in gene_list[ gene ].new_alleles() ] )
#         # for phen in gene_list[gene].novel_phenotypes():
#         #     print(phen)
#         # for phen in gene_list[gene].phenotypes:
#         #     print( phen )
#         # for phen in gene_list[gene].new_phenotypes:
#         #     print("N", phen )
#
#     if gene_list[ gene ].has_known_alleles():
#         print( gene_list[gene].name , gene_list[gene].marker_id , end = '\t')
#         print( "known" , [ al for al in gene_list[ gene ].known_alleles() ] )
#         # for phen in gene_list[gene].novel_phenoyptes():
#         #     print(phen)
#         # for phen in gene_list[gene].phenotypes:
#         #     print( phen )
#         # for phen in gene_list[gene].new_phenotypes:
#         #     print("N", phen )
#

# c = 0
# for gene in sorted(gene_list):
#     if not gene_list[ gene ].are_all_allele_new() :
#          print(gene , len(gene_list[gene].novel_phenotypes() ) )
#          c += 1
#
# print( "Number of new genes ", c , " total ")

# print( "\"","\", \"".join( gene_list["Anks6"].phenotypes | gene_list["Anks6"].new_phenotypes),"\"", sep="" )



print( "gene", "marker_id","novel", "n_mgi_alleles", "n_impc_alleles", "n_mgi_mpterms", "n_impc_mpterms",\
    "n_known_mpterms", "n_novel_mpterms", sep = "\t")

for gene in gene_list:
    n_mgi_mpterms  = count_mp( mgi[gene][2]  ) if gene in mgi else 0
    n_impc_mpterms = count_mp( impc[gene][2] ) if gene in impc else 0
    n_known_mpterms = len( gene_list[gene].phenotypes )
    n_new_phenotypes = len( gene_list[gene].new_phenotypes )
    is_novel = 'Y' if gene_list[gene].are_all_allele_new() else 'N'
    print( gene, gene_list[gene].marker_id, is_novel,len(gene_list[gene].alleles_mgi ),\
         len(gene_list[gene].alleles_impc ),n_mgi_mpterms, n_impc_mpterms, n_known_mpterms,  \
         n_new_phenotypes, sep = "\t")



sys.exit("")
