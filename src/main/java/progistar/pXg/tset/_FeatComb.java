package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import progistar.pXg.constants.Parameters;

public class _FeatComb {

	public static String[] INIT_FEATS = {
			"SpecId",
			"Label",
			"ScanNr",
			"MainScore",	
			"Log2Reads",
			"Charge2",
			"Charge3",
			"Charge4",
			"Mass",
			"ppm",
			"SNV",
			"INDEL",
			"PC",
			"FS",
			"5`-UTR",
			"3`-UTR",
			"asRNA",
			"IGR",
			"IR",
			"ncRNA",
			"AS",
			"unknown",
			"Length",
			"Rank",
			"Peptide",
			"Proteins"};
	
	
	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/p001/S1.RAW.PEAKS.p001.pxg.pin");
		genTypeFeat(file);
	}
	
	public static void genTypeFeat (File file) throws IOException {

		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".pin", ".feat1.pin")));
		String line = null;
		String header = BR.readLine();
		
		String[] initHeaders = header.split("\t");
		
		int snvIdx = -1;
		int indelIdx = -1;
		int fsIdx = -1;
		int utr5Idx = -1;
		int utr3Idx = -1;
		int asRNAIdx = -1;
		int ncRNAIdx = -1;
		int IGRIdx = -1;
		int IRIdx = -1;
		int ASIdx = -1;
		int unknownIdx = -1;
		
		for(int i=0; i<initHeaders.length; i++) {
			if(initHeaders[i].equalsIgnoreCase("SNV")) {
				snvIdx = i;
			} else if(initHeaders[i].equalsIgnoreCase("INDEL")) {
				indelIdx = i;
			}  else if(initHeaders[i].equalsIgnoreCase("IR")) {
				IRIdx = i;
			}  else if(initHeaders[i].equalsIgnoreCase("asRNA")) {
				asRNAIdx = i;
			}  else if(initHeaders[i].equalsIgnoreCase("IGR")) {
				IGRIdx = i;
			}  else if(initHeaders[i].equalsIgnoreCase("FS")) {
				fsIdx = i;
			}  else if(initHeaders[i].equalsIgnoreCase("3`-UTR")) {
				utr3Idx = i;
			}  else if(initHeaders[i].equalsIgnoreCase("5`-UTR")) {
				utr5Idx = i;
			} else if(initHeaders[i].equalsIgnoreCase("AS")) {
				ASIdx = i;
			} else if(initHeaders[i].equalsIgnoreCase("unknown")) {
				unknownIdx = i;
			} else if(initHeaders[i].equalsIgnoreCase("ncRNA")) {
				ncRNAIdx = i;
			}
			
		}
		
		header = "SpecId\tLabel\tScanNr\tMainScore\tLog2Reads\tCharge2\tCharge3\tCharge4\tppm\tLength\tIsCanonical\tPeptide\tProteins";
		BW.append(header);
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			double isCryptic = 0;
			int snv = Integer.parseInt(fields[snvIdx]);
			int indel = Integer.parseInt(fields[indelIdx]);
			int fs = Integer.parseInt(fields[fsIdx]);
			int utr5 = Integer.parseInt(fields[utr5Idx]);
			int utr3 = Integer.parseInt(fields[utr3Idx]);
			int asRNA = Integer.parseInt(fields[asRNAIdx]);
			int IGR = Integer.parseInt(fields[IGRIdx]);
			int IR = Integer.parseInt(fields[IRIdx]);
			int ncRNA = Integer.parseInt(fields[ncRNAIdx]);
			int AS = Integer.parseInt(fields[ASIdx]);
			int unknown = Integer.parseInt(fields[unknownIdx]);
			
			isCryptic = (int) (snv * Parameters.PENALTY_MUTATION +
					indel * Parameters.PENALTY_MUTATION +
					fs * Parameters.PENALTY_FS +
					utr5 * Parameters.PENALTY_5UTR +
					utr3 * Parameters.PENALTY_3UTR +
					asRNA * Parameters.PENALTY_asRNA +
					IGR * Parameters.PENALTY_IGR +
					IR * Parameters.PENALTY_IR +
					ncRNA * Parameters.PENALTY_ncRNA +
					AS * Parameters.PENALTY_AS +
					unknown * Parameters.PENALTY_UNMAP);
			
			if(isCryptic > 0) {
				//isCryptic = 1+((double)Math.log(isCryptic+1)/Math.log(2));
				isCryptic = -1;
			} else {
				isCryptic = 1;
			}
			
			StringBuilder output = new StringBuilder();
			for(int i=0; i<10; i++) {
				if(i == 8) continue; //pass mass
				output.append(fields[i]).append("\t");
			}
			output.append(isCryptic).append("\t");
			output.append(fields[22]).append("\t");
			output.append(fields[24]).append("\t");
			output.append(fields[25]);
			BW.append(output.toString());
			BW.newLine();
		}
		
		BW.close();
		BR.close();
	}
}
