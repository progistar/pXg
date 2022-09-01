package progistar.thirdparty.msgf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

class HFasta {
	String header = null;
	String sequence = "";
}

public class HistoneID {

	public static void main(String[] args) throws IOException {
		String fileName = "/Users/gistar/projects/pXg/Sturm_JPR2021/PeptideList.tsv";
		File file = new File(fileName);
		
		BufferedReader histoneFile = new BufferedReader(new FileReader("/Users/gistar/projects/pXg/HistoneOnly.fasta"));
		String line = null;
		
		ArrayList<HFasta> histoneRecords = new ArrayList<HFasta>();
		HFasta hFasta = new HFasta();
		while((line = histoneFile.readLine()) != null) {
			if(line.startsWith(">")) {
				hFasta = new HFasta();
				histoneRecords.add(hFasta);
				hFasta.header = line;
			} else {
				hFasta.sequence += line;
			}
		}
		
		for(HFasta h : histoneRecords) {
			h.sequence = h.sequence.replace("I", "L");
		}
		histoneFile.close();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".tsv", ".histoneID.tsv")));
		
		
		String header = BR.readLine();
		
		BW.append(header).append("\tHistoneID");
		BW.newLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[1].replaceAll("[+-0123456789*]", "").replace("I", "L");
			ArrayList<String> hHeader = new ArrayList<String>();
			
			for(HFasta h : histoneRecords) {
				if(h.sequence.contains(peptide)) {
					hHeader.add(h.header);
				}
			}
			
			BW.append(line).append("\t").append(hHeader.size()+"");
			BW.newLine();
		}
		
		BR.close();
		BW.close();
	}
}
