package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.thirdparty.netMHCpan.NetMHCpanParser;
import progistar.thirdparty.netMHCpan.NetMHCpanResult;

public class Q2_2_ReadTranslationTest {

	public static void main(String[] args) throws IOException {
		
		NetMHCpanResult res = NetMHCpanParser.parseNetMHCpan("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.rmEnz.netmhcpan.xls");
		File[] files = {new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.decoy.psm"), 
				new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.target.psm")};
		File fasta = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/2204_Human_HHV4_179Contam_UniProteome.fasta");
		Hashtable[] canonicalMap = loadDB(fasta);
		
		for(File file : files) {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".psm", ".classMark.psm")));
			int idx = 0;
			if(file.getName().contains("decoy")) {
				idx = 1;
			}
			
			String line = null;
			String header = BR.readLine();
			
			BW.append("Class\t"+header+"\t"+res.getHeader());
			BW.newLine();
			StringBuilder proteins = new StringBuilder();
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String oriPeptide = fields[4].split("\\.")[1];
				String peptide = oriPeptide.replace("I", "L");
				proteins.setLength(0);
				
				for(int i=5; i<fields.length; i++) {
					if(i!=5) {
						proteins.append(";");
					}
					proteins.append(fields[i]);
				}
				
				if(canonicalMap[idx].get(peptide) == null) {
					BW.append("Noncanonical");
				}else {
					BW.append("Canonical");
				}
				
				for(int i=0; i<5; i++) {
					BW.append("\t"+fields[i]);
				}
				BW.append("\t"+proteins.toString()).append("\t"+res.getHLATyping(oriPeptide));
				
				BW.newLine();
			}
			
			
			BW.close();
			BR.close();
		}
		
	}
	
	public static Hashtable[] loadDB (File file) throws IOException {
		Hashtable[] map = new Hashtable[2];
		map[0] = new Hashtable<String, String>();
		map[1] = new Hashtable<String, String>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		ArrayList<String> targets = new ArrayList<String>();
		ArrayList<String> decoys = new ArrayList<String>();
		StringBuilder sequence = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(sequence.length() != 0) {
					targets.add(sequence.toString());
					decoys.add(sequence.reverse().toString());
				}
				sequence.setLength(0);
			} else {
				sequence.append(line.replace("I", "L"));
			}
		}
		if(sequence.length() != 0) {
			targets.add(sequence.toString());
			decoys.add(sequence.reverse().toString());
		}
		
		BR.close();
		System.out.println("Done to read database");
		// put target first
		for(int i=0; i<targets.size(); i++) {
			String record = targets.get(i);
			for(int j=0; j<8; j++) {
				int len = 8+j;
				for(int k=len; k<=record.length(); k++) {
					map[0].put(record.substring(k-len, k),"target");
				}
			}
		}
		System.out.println("Done to put targets");
		// put decoy
		for(int i=0; i<decoys.size(); i++) {
			String record = decoys.get(i);
			for(int j=0; j<8; j++) {
				int len = 8+j;
				for(int k=len; k<=record.length(); k++) {
					map[1].put(record.substring(k-len, k),"decoy");
				}
			}
		}
		System.out.println("Done to put decoys");
		
		return map;
	}
}
