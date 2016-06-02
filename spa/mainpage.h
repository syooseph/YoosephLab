/**
 * @mainpage SPA: <em>Short Protein Sequence Assembler</em>
 * This application assembles short <em>protein</em> sequences. \n
 * This application suit includes <b>three</b> programs.
 *
 * \section sec1 Post processing of gene calling program.
 * -# ORFs are predecited by FragGeneScan.
 * -# Any ORFs with INDELs are disregarded by default.
 * -# Search start/stop codons. If found, corresponding amino acid or
 * 	terminator is inserted.
 *
 * \section sec2 Preprocessing of short protein sequences.
 * -# Generate binary graph input for assembler.
 * -# Generate inverted index for assember. 
 *
 * \section sec3 Sequence Assembler.
 * -# Load sequences.
 * -# Build graph and inverted index.
 * -# Assemble sequences.
 * -# Write result.
*/

