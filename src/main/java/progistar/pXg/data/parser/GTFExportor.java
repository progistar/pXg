package progistar.pXg.data.parser;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.PBlock;
import progistar.pXg.data.XBlock;

public class GTFExportor {

	private GTFExportor() {
		
	}
	
	public static void writeGTF (PBlock pBlock, XBlock xBlock, BufferedWriter BW) throws IOException {
		String peptide = xBlock.peptideSequence;
		String nucleotide = xBlock.genomicSequence;
		String[] genomicLoci = xBlock.genomicLocus.split("\\|");
		char strand = xBlock.strand;
		String score = pBlock.score+"";
		String scanID = pBlock.getScanID();
		String rank = pBlock.rank+"";
		String mutation = xBlock.mutations;
		String events = xBlock.toEvents().get("key");
		String fastaID = xBlock.toFastaIDs().get("key");
		
		String chr = genomicLoci[0].split("\\:")[0];
		String startPos = null;
		String endPos = null;
		
		for(String genomicLocus : genomicLoci) {
			genomicLocus = genomicLocus.replace(chr+":", "");
			if(startPos == null) {
				startPos = genomicLocus.split("\\-")[0];
			}
			endPos = genomicLocus.split("\\-")[1];
		}
		
		BW.append(chr).append("\t"); // chr
		BW.append(Constants.VERSION).append("\t"); // source
		BW.append("transcript").append("\t"); // feature
		BW.append(startPos).append("\t"); // start
		BW.append(endPos).append("\t"); // end
		BW.append(score).append("\t"); // score
		BW.append(strand).append("\t"); // strand
		BW.append(".").append("\t"); // frame
		BW.append("gene_id \"").append(peptide+"("+xBlock.genomicLocus+")").append("\";"); // attributes: gene_id
		BW.append(" gene_name \"").append(peptide+"("+events+")").append("\";"); // attributes: gene_name
		BW.append(" transcript_id \"").append(scanID).append("\";"); // attributes: transcript_id
		BW.append(" event \"").append(events).append("\";"); // attributes: event
		BW.append(" mutation \"").append(mutation).append("\";"); // attributes: mutation
		BW.append(" nucleotide \"").append(nucleotide).append("\";"); // attributes: nucleotide
		BW.append(" rank \"").append(rank).append("\";"); // attributes: rank
		BW.append(" score \"").append(score).append("\";"); // attributes: score
		BW.append(" fastaID \"").append(fastaID).append("\";"); // attributes: fasta matched list
		BW.newLine();
		
		ArrayList<String> locusInfo = new ArrayList<String>();
		if(strand == '+') {
			for(int i=0; i<genomicLoci.length; i++) {
				locusInfo.add(genomicLoci[i]);
			}
		} else {
			for(int i=genomicLoci.length-1; i>=0; i--) {
				locusInfo.add(genomicLoci[i]);
			}
		}
		
		for(int i=0; i<locusInfo.size(); i++) {
			String locus = locusInfo.get(i);
			locus = locus.replace(chr+":", "");
			startPos = locus.split("\\-")[0];
			endPos = locus.split("\\-")[1];
			
			BW.append(chr).append("\t"); // chr
			BW.append(Constants.VERSION).append("\t"); // source
			BW.append("exon").append("\t"); // feature
			BW.append(startPos).append("\t"); // start
			BW.append(endPos).append("\t"); // end
			BW.append(score).append("\t"); // score
			BW.append(strand).append("\t"); // strand
			BW.append(".").append("\t"); // frame
			BW.append("gene_id \"").append(peptide+"("+xBlock.genomicLocus+")").append("\";"); // attributes: gene_id
			BW.append(" gene_name \"").append(peptide+"("+events+")").append("\";"); // attributes: gene_name
			BW.append(" transcript_id \"").append(scanID).append("\";"); // attributes: transcript_id
			BW.append(" exon_id \"").append((i+1)+"").append("\";"); // attributes: exon_id
			BW.newLine();
		}
		
	}
}
