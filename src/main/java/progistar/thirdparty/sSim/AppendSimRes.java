package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class AppendSimRes {

	public static void main(String[] args) throws IOException {
		File simResFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/ProteomToolsHLA/SimRes.tsv");
		File ncPSMFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/ProteomToolsHLA/UniqueNoncanonicalPSMs.tsv");

		BufferedWriter BW = new BufferedWriter(new FileWriter("simRes.annotation.tsv"));
		BufferedReader BR = new BufferedReader(new FileReader(ncPSMFile));
		String line = null;

		Hashtable[] reads = new Hashtable[4];
		Hashtable[] psms = new Hashtable[4];

		for(int i=0; i<reads.length; i++) {
			reads[i] = new Hashtable<String, String>();
			psms[i] = new Hashtable<String, Integer>();
		}


		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[20];
			String charge = fields[10];
			String read = fields[35];
			String subject = fields[46];
			int subjectIdx = Integer.parseInt(subject.split("\\s")[1]) - 1;

			String key = peptide+"_"+charge;
			reads[subjectIdx].put(key, read);
			Integer count = 0;
			if(psms[subjectIdx].get(key) != null) {
				count = (Integer) psms[subjectIdx].get(key);
			}
			count++;
			psms[subjectIdx].put(key, count);
		}

		BR.close();

		BR = new BufferedReader(new FileReader(simResFile));
		String header = BR.readLine() +"\tReads for subject 1\tReads for subject 2\tReads for subject 3\tReads for subject 4"
				+ "\tPSMs for subject 1\tPSMs for subject 2\tPSMs for subject 3\tPSMs for subject 4";

		BW.append(header);
		BW.newLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			BW.append(line);
			String key = fields[0];
			// write reads
			for (Hashtable read : reads) {
				BW.append("\t");
				if(read.get(key) == null) {
					BW.append("0");
				} else {
					BW.append((String) read.get(key));
				}
			}

			// write PSMs
			for(int i=0; i<reads.length; i++) {
				BW.append("\t");
				if(psms[i].get(key) == null) {
					BW.append("0");
				} else {
					BW.append((psms[i].get(key))+"");
				}
			}
			BW.newLine();
		}


		BR.close();
		BW.close();
	}
}
