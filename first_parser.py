import sys
import re

_test_lim = 1E10

mp_pattern = re.compile(r'^MP:[0-9]*')

class novel_gene( object ):

    def __init__(self, name, marker_id):
        self.name = name
        self.marker_id = marker_id
        self.alleles_mgi = set()
        self.alleles_impc = set()

    def find_alleles_in_mgi(self, l ):
        if l:
            for elem in l:
                self.alleles_mgi.add( elem )

    def find_alleles_in_impc(self, l ):
        if l :
            for elem in l:
                self.alleles_impc.add( elem )

    def new_alleles( self ):
        # alleles in impc that match (or are absent) in mgi
        not_even_in_mgi = self.alleles_impc - self.alleles_mgi

        intersection = self.alleles_impc & self.alleles_mgi

        return not_even_in_mgi | intersection # union of two sets

    def are_all_allele_new( self ):
        new = self.new_alleles()
        ans = len( self.alleles_mgi - new ) == 0
        return ans

    def known_alleles( self ):
        return  self.alleles_mgi - self.alleles_impc

    def has_known_alleles( self ):
        return len( self.known_alleles() ) > 0



def is_mp_term( elem ):
    return re.match(mp_pattern, elem )


mgi = { }
with open(sys.argv[1] , 'r' ) as f:
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
    else:
        gene_list[ gene ].find_alleles_in_mgi( None )

    if gene in impc:
        gene_list[ gene ].find_alleles_in_impc( impc[ gene ][0] )
    else:
        gene_list[ gene ].find_alleles_in_impc( None )


c = 0
for gene in gene_list:
#    print( gene_list[gene].name )
    if gene_list[ gene ].are_all_allele_new():
        c += 1
        print( gene_list[gene].name , gene_list[gene].marker_id , end = '\t')
        print( "new", [ al for al in gene_list[ gene ].new_alleles() ] )

    if gene_list[ gene ].has_known_alleles():
        print( gene_list[gene].name , gene_list[gene].marker_id , end = '\t')
        print( "known" , [ al for al in gene_list[ gene ].known_alleles() ] )

print( "Number of new genes ", c , " total ")

# print( gene_list["Ssr2"].alleles_impc )
# print( gene_list["Ssr2"].alleles_mgi )

sys.exit("")
