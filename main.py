# Fungsi untuk mengonversi DNA ke RNA komplementer (RNA1)
def dna_to_rna_complement(dna_sequence):
    # Peta komplementer DNA -> RNA (A -> U, T -> A, C -> G, G -> C)
    complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    rna_complement = ''.join([complement[nucleotide] for nucleotide in dna_sequence])
    return rna_complement

# Fungsi untuk mengganti T dengan U dalam RNA komplementer (RNA2)
def convert_thymine_to_uracil(rna_sequence):
    return rna_sequence.replace('T', 'U')

# Fungsi untuk menerjemahkan RNA ke urutan asam amino
def translate_rna_to_protein(rna_sequence):
    # Tabel kodon RNA untuk kode genetik standar
    codon_table = {
        'AUG':'M', 'UAA':'_', 'UAG':'_', 'UGA':'_', 'UUU':'F', 'UUC':'F',
        'UUA':'L', 'UUG':'L', 'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
        'UAU':'Y', 'UAC':'Y', 'UGU':'C', 'UGC':'C', 'UGG':'W', 'CUU':'L',
        'CUC':'L', 'CUA':'L', 'CUG':'L', 'CCU':'P', 'CCC':'P', 'CCA':'P',
        'CCG':'P', 'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGU':'R',
        'CGC':'R', 'CGA':'R', 'CGG':'R', 'AUU':'I', 'AUC':'I', 'AUA':'I',
        'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'AAU':'N', 'AAC':'N',
        'AAA':'K', 'AAG':'K', 'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
        'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V', 'GCU':'A', 'GCC':'A',
        'GCA':'A', 'GCG':'A', 'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }
    
    protein_sequence = ''
    
    # Iterasi melalui urutan RNA dalam langkah 3 untuk setiap kodon
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if codon in codon_table:
            protein_sequence += codon_table[codon]
    
    return protein_sequence

# Contoh penggunaan
dna_input = input("Enter DNA sequence: ").upper()

# Transkripsi DNA ke RNA komplementer (RNA1)
rna1 = dna_to_rna_complement(dna_input)
# Ubah T menjadi U dalam RNA komplementer (RNA2)
rna2 = convert_thymine_to_uracil(rna1)
# Terjemahkan RNA2 ke urutan asam amino
protein_output = translate_rna_to_protein(rna2)

print("RNA1: ", rna1)  # RNA komplementer pertama
print("RNA2: ", rna2)  # RNA setelah mengganti T dengan U
print("Amino Acid Sequence: ", protein_output)
