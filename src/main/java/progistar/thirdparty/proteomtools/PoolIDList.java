package progistar.thirdparty.proteomtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class PoolIDList {

	public static void main(String[] args) throws IOException {
		File expFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/ProteomToolsHLA/ProteomToolsHLAList.tsv");
		File idFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/ProteomToolsHLA/ProteomToolsIDList.tsv");
		File pxgFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/ProteomToolsHLA/UniqueNoncanonicalPSMs.tsv");
		Hashtable<String, String> experiementalMapper = loadMapper(expFile, 1);
		Hashtable<String, String> idMapper = loadMapper(idFile, 0);
		
		BufferedReader BR = new BufferedReader(new FileReader(pxgFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pxgFile.getAbsolutePath().replace(".tsv", ".ptools.tsv")));
		String line = null;
		
		int peptideIdx = 20;
		int subjectIdx = 46;
		int alcIdx = 7;
		ArrayList<String> records = new ArrayList<String>();
		Hashtable<String, ArrayList<String>> subjectList = new Hashtable<String, ArrayList<String>>();
		
		while((line = BR.readLine()) != null) {
			if(records.size() != 0) {
				String[] fields = line.split("\t");
				String peptide = fields[peptideIdx];
				String subject = fields[subjectIdx];
				
				ArrayList<String> subjects = subjectList.get(peptide);
				
				if(subjects == null) {
					subjects = new ArrayList<String>();
					subjectList.put(peptide, subjects);
				}
				
				subjects.add(subject+"_"+fields[alcIdx]);
				
			}
			records.add(line);
		}
		
		// subject aggregation
		Hashtable<String, String> isDuplicated = new Hashtable<String, String>();
		BW.append(records.get(0)+"\tSubject 1\tSubject 2\tSubject 3\tSubject 4");
		BW.newLine();
		
		String[] subjectNames = {"Subject 1", "Subject 2", "Subject 3", "Subject 4"};
		for(int i=1; i<records.size(); i++) {
			line = records.get(i);
			String[] fields = line.split("\t");
			String peptide = fields[peptideIdx];
			
			if(isDuplicated.get(peptide) != null) {
				continue;
			}
			isDuplicated.put(peptide, "");
			BW.append(line);
			
			ArrayList<String> subjects = subjectList.get(peptide);
			for(int j=0; j<subjectNames.length; j++) {
				StringBuilder list = new StringBuilder();
				for(String subject : subjects) {
					String[] subjectField = subject.split("\\_");
					if(subjectField[0].equalsIgnoreCase(subjectNames[j])) {
						if(list.length() != 0) {
							list.append(";");
						}
						list.append(subjectField[1]);
					}
				}
				if(list.length() == 0) {
					list.append("NA");
				}
				BW.append("\t").append(list.toString());
			}
			
			String exp = experiementalMapper.get(peptide);
			String id = null;
			if(exp != null) {
				id = idMapper.get(peptide);
			} else {
				exp = "NA";
			}
			
			BW.append("\t").append(exp);
			if(id != null) {
				BW.append("\t").append(id);
			}
			BW.newLine();
		}
		
		
		BW.close();
		BR.close();
	}
	
	public static Hashtable<String, String> loadMapper (File file, int keyIdx) throws IOException {
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[keyIdx];
			
			if(mapper.get(key) == null) {
				mapper.put(key, line);
			} else {
				System.out.println(key +" is duplicated !");
			}
		}
		
		BR.close();
		
		return mapper;
	}
}
