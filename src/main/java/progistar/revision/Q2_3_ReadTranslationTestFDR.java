package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

class Record implements Comparable<Record>{
	boolean isCanonical;
	boolean isTarget;
	String line;
	double score;

	@Override
	public int compareTo(Record o) {
		if(this.score > o.score) {
			return -1;
		}else if(this.score < o.score) {
			return 1;
		}
		return 0;
	}


}



public class Q2_3_ReadTranslationTestFDR {
	public static String header = null;

	public static void main(String[] args) throws IOException {
		File targetFile = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/1.SepFDR/B_LCL1.target.classMark.psm");
		File decoyFile = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/1.SepFDR/B_LCL1.decoy.classMark.psm");
		File result = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.noncanonical.fdr5");


		ArrayList<Record> targets = loadFile(targetFile, false, true);
		ArrayList<Record> decoys = loadFile(decoyFile, false, false);

		BufferedWriter BW = new BufferedWriter(new FileWriter(result));
		BW.append(header);
		BW.newLine();

		ArrayList<Record> psms = targets;
		psms.addAll(decoys);

		Collections.sort(psms);

		double targetN = 0;
		double decoyN = 0;

		int fdrIdx = 0;
		for(int i=0; i<psms.size(); i++) {
			Record record = psms.get(i);
			if(record.isTarget) {
				targetN++;
			} else {
				decoyN++;
			}

			if(targetN != 0 && decoyN/targetN < 0.05) {
				fdrIdx = i;
			}
		}

		targetN = 0;
		decoyN = 0;

		for(int i=0; i<=fdrIdx; i++) {
			Record record = psms.get(i);

			if(record.isTarget) {
				BW.append(record.line);
				BW.newLine();
				targetN++;
			} else {
				decoyN++;
			}
		}

		System.out.println(targetN+"/"+decoyN);
		BW.close();
	}

	public static ArrayList<Record> loadFile (File file, boolean hasCanonical, boolean isTarget) throws IOException {
		ArrayList<Record> records = new ArrayList<>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;

		header = BR.readLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String canonical = fields[0];
			String mhc = fields[13];
			boolean isCanonical = true;
			double score = Double.parseDouble(fields[2]);
			if(canonical.startsWith("Non")) {
				isCanonical = false;
			}

			if(mhc.equalsIgnoreCase("NB") || (isCanonical != hasCanonical)) {
				continue;
			}

			Record record = new Record();
			record.isCanonical = isCanonical;
			record.line = line;
			record.score = score;
			record.isTarget = isTarget;

			records.add(record);
		}
		BR.close();
		return records;
	}
}
