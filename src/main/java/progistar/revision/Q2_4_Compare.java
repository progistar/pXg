package progistar.revision;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

class CompRecord {
	String spec;
	String line;
	String peptide;
	boolean isCanonical;

}

public class Q2_4_Compare {

	public static void main(String[] args) throws IOException {
		File cometFile = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/2.Compare/B_LCL1.comet.fdr5");
		File pXgFile = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/2.Compare/B_LCL1.pXg.res");

		BufferedReader BR = new BufferedReader(new FileReader(cometFile));
		String line = BR.readLine();

		ArrayList<CompRecord> cometRecords = new ArrayList<>();
		ArrayList<CompRecord> pXgRecords = new ArrayList<>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String spec = fields[1];
			String cls = fields[0];
			String peptide = fields[5].split("\\.")[1];
			boolean isCanonical = true;
			if(cls.startsWith("Non")) {
				isCanonical = false;
			}

			CompRecord record = new CompRecord();
			record.isCanonical = isCanonical;
			record.line = line;
			record.peptide = peptide;
			record.spec = spec;

			cometRecords.add(record);
		}

		BR.close();

		BR = new BufferedReader(new FileReader(pXgFile));
		BR.readLine();

		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String spec = fields[0].split("\\.")[0];
			String scan = fields[0].split("\\|")[1].split("\\:")[1];
			String charge = fields[0].split("\\|")[2];

			spec = spec+"."+scan+"."+scan+"."+charge;

			String cls = fields[39];
			String peptide = fields[22];
			boolean isCanonical = true;
			if(cls.startsWith("F")) {
				isCanonical = false;
			}

			CompRecord record = new CompRecord();
			record.isCanonical = isCanonical;
			record.line = line;
			record.peptide = peptide;
			record.spec = spec;

			pXgRecords.add(record);
		}

		BR.close();

		int[] overlap = new int[2];
		Hashtable<String, String> pXgNC = loadPept(pXgRecords, false);
		Hashtable<String, String> pXgC = loadPept(pXgRecords, true);
		Hashtable<String, String> cometNC = loadPept(cometRecords, false);
		Hashtable<String, String> cometC = loadPept(cometRecords, true);

		pXgC.forEach((pept, nil)->{
			if(cometC.get(pept) != null) {
				overlap[0]++;
			}
		});

		cometNC.forEach((pept, nil)->{
			if(pXgNC.get(pept) != null) {
				overlap[1]++;
			} else {
				System.out.println(nil);
			}
		});

		System.out.println("Canonical");
		System.out.println(pXgC.size()+"\t"+cometC.size()+"\t"+overlap[0]);
		System.out.println("Noncanonical");
		System.out.println(pXgNC.size()+"\t"+cometNC.size()+"\t"+overlap[1]);
	}

	public static Hashtable<String, String> loadPept(ArrayList<CompRecord> records, boolean isCanonical) {
		Hashtable<String, String> map = new Hashtable<>();

		for (CompRecord record : records) {
			if(record.isCanonical == isCanonical) {
				//map.put(record.peptide+"_"+record.spec, record.line);
				//map.put(record.spec, record.line);
				map.put(record.peptide, record.line);
			}
		}

		return map;
	}
}
